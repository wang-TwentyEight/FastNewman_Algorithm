#include <iostream>
#include <fstream>
#include <string>
#include "stdlib.h"
#include "time.h"
#include "math.h"

#include "maxheap.h"
#include "vektor.h"
//#define MAX_VERTEX_NUM 34

using namespace std;

// ------------------------------------------------------------------------------------
// Edge object - defined by a pair of vertex indices and *edge pointer to next in linked-list
class edge {
public:
	int     so;					// originating node
	int     si;					// terminating node
	edge    *next;					// pointer for linked list of edges
	
	edge();						// default constructor
	~edge();						// default destructor
};
edge::edge()  { so = 0; si = 0; next = NULL; }
edge::~edge() {}

// ------------------------------------------------------------------------------------
// Nodenub object - defined by a *node pointer and *node pointer 
struct nodenub {
	tuple	*heap_ptr;			// pointer to node(max,i,j) in max-heap of row maxes
	vektor    *v;					// pointer stored vector for (i,j)
};

// ------------------------------------------------------------------------------------
// tuple object - defined by an real value and (row,col) indices
#if !defined(TUPLE_INCLUDED)
#define TUPLE_INCLUDED
struct tuple {
	double    m;					// stored value
	int		i;					// row index
	int		j;					// column index
	int		k;					// heap index
};
#endif

// ordered pair structures (handy in the program)
struct apair { int x; int y; };
#if !defined(DPAIR_INCLUDED)
#define DPAIR_INCLUDED
class dpair {
public:
	int x; double y; dpair *next;
	dpair(); ~dpair();
};
dpair::dpair()  { x = 0; y = 0.0; next = NULL; }
dpair::~dpair() {}
#endif

// ------------------------------------------------------------------------------------
// List object - simple linked list of integers
class list {
public:
	int		index;				// node index
	list		*next;				// pointer to next element in linked list
	list();   ~list();
};
list::list()  { index= 0; next = NULL; }
list::~list() {}

// ------------------------------------------------------------------------------------
// Community stub object - stub for a community list
class stub {
public:
	bool		valid;				// is this community valid?
	int		size;				// size of community
	list		*members;				// pointer to list of community members
	list		*last;				// pointer to end of list
	stub();   ~stub();
};
stub::stub()  { valid = false; size = 0; members = NULL; last = NULL; }
stub::~stub() {
	list *current;
	if (members != NULL) {
		current = members;
		while (current != NULL) { members = current->next; delete current; current = members; }
	}
}
// ------------------------------------------------------------------------------------
// 函数声明 --------------------------------------------------------------

void buildDeltaQMatrix();
//void buildFilenames();
void dqSupport();
void groupListsSetup();
void groupListsStats();
void groupListsUpdate(const int x, const int y);
void mergeCommunities(int i, int j);
void readAndwrite();
void readInputFile();

// ------------------------------------------------------------------------------------
// 程序参数 -----------------------------------------------------------------

struct netparameters {
	int			n;				// number of nodes in network
	int			m;				// number of edges in network
	int			maxid;			// maximum node id
	int			minid;			// minimum node id
}; netparameters    gparm;

struct groupstats {
	int			numgroups;		// number of groups
	double		meansize;			// mean size of groups
	int			maxsize;			// size of largest group
	int			minsize;			// size of smallest group
	double		*sizehist;		// distribution of sizes
}; groupstats		gstats;

// ------------------------------------------------------------------------------------
// ----------------------------------- 全局变量 -------------------------------


edge		*e;				// initial adjacency matrix (sparse)
edge		*elist;			// list of edges for building adjacency matrix
nodenub   *dq;				// dQ matrix
maxheap   *h;				// heap of values from max_i{dQ_ij}
double    *Q;				// Q(t)
dpair     Qmax;			// maximum Q value and the corresponding time t
double    *a;				// A_i
apair	*joins;			// list of joins
stub		*c;				// link-lists for communities

enum {NONE};

int    supportTot;
double supportAve;





void main(){
	
	int i;
	edge   *current;
	readInputFile();
	/*for(i =0;i<gparm.n;i++){
		current  = &e[i];
		while(current!=NULL){

		cout <<" e: ["<< i <<"] :"<< current->so<<"\t"<< current->si<<"\t"<<endl;
		current=current->next;
		}
	}*/
				// Allocate data structures for main loop
	a     = new double [gparm.maxid];
	Q     = new double [gparm.n+1];
	joins = new apair  [gparm.n+1];
	for (int i=0; i<gparm.maxid; i++) { a[i] = 0.0; }
	for (int i=0; i<gparm.n+1;   i++) { Q[i] = 0.0; joins[i].x = 0; joins[i].y = 0; }
	int t = 1;                       //*****************
	Qmax.y = -4294967296.0;  Qmax.x = 0;
	groupListsSetup(); 		// will need to track agglomerations

	cout << "now building initial dQ[]" << endl;
	buildDeltaQMatrix();							// builds dQ[] and h
	/*dpair  *list,*current1;
	for(i =0;i<gparm.n;i++){
		list     = dq[i].v ->returnTreeAsList();
		current1 = list;
		while(current1!=NULL){
			cout<<"dq["<<i<<"]‘s dpairs"<<"X:"<<current1->x<<"Y:"<<current1->y<<endl;
			current1=current1->next;
		}
	
	}*/
	// initialize f_joins, f_support files
	ofstream fjoins("f_joints", ios::in);
	fjoins << -1 << "\t" << -1 << "\t" << Q[0] << "\t0\n";
	fjoins.close();
	dqSupport();
	cout<< 0 << "\t" << supportTot << "\t" << supportAve << "\t" << 0 << "\t->\t" << 0 << "\n";



	// ----------------------------------------------------------------------
	// Start FastCommunity algorithm
	cout << "starting algorithm now." << endl;
	tuple  dQmax, dQnew;
	int isupport, jsupport;
	while (h->heapSize() > 2){
		// ---------------------------------
		// Find largest dQ
		 h->printHeapTop10(); cout << endl; 
		dQmax = h->popMaximum();					// select maximum dQ_ij // convention: insert i into j
		if (dQmax.m < -4000000000.0) { break; }		// no more joins possible
		cout << "Q["<<t-1<<"] = "<<Q[t-1];
		
		// ---------------------------------
		// Merge the chosen communities
		cout << "\tdQ = " << dQmax.m << "\t  |H| = " << h->heapSize() << "\n";
		if (dq[dQmax.i].v == NULL || dq[dQmax.j].v == NULL) {
			cout << "WARNING: invalid join (" << dQmax.i << " " << dQmax.j << ") found at top of heap\n";// cin >> pauseme;
		}
		isupport = dq[dQmax.i].v->returnNodecount();
		jsupport = dq[dQmax.j].v->returnNodecount();
		if (isupport < jsupport) {
			cout << "  join: " << dQmax.i << " -> " << dQmax.j << "\t";
			cout << "(" << isupport << " -> " << jsupport << ")\n";
			mergeCommunities(dQmax.i, dQmax.j);	// merge community i into community j
			joins[t].x = dQmax.i;				// record merge of i(x) into j(y)
			joins[t].y = dQmax.j;				// 
		} else {								// 
			cout << "  join: " << dQmax.i << " <- " << dQmax.j << "\t";
			cout << "(" << isupport << " <- " << jsupport << ")\n";
			dq[dQmax.i].heap_ptr = dq[dQmax.j].heap_ptr; // take community j's heap pointer
			dq[dQmax.i].heap_ptr->i = dQmax.i;			//   mark it as i's
			dq[dQmax.i].heap_ptr->j = dQmax.j;			//   mark it as i's
			mergeCommunities(dQmax.j, dQmax.i);	// merge community j into community i
			joins[t].x = dQmax.j;				// record merge of j(x) into i(y)
			joins[t].y = dQmax.i;				// 
		}									// 
		Q[t] = dQmax.m + Q[t-1];	// record Q(t)


		// ---------------------------------
		// Record the support data to file
		ofstream fjoins("f_joints", ios::app);
		fjoins << joins[t].x-1 << "\t" << joins[t].y-1 << "\t";	// convert to external format
		if ((Q[t] > 0.0 && Q[t] < 0.0000000000001) || (Q[t] < 0.0 && Q[t] > -0.0000000000001))
			{ fjoins << 0.0; } else { fjoins << Q[t]; }
		fjoins << "\t" << t << "\n";
		fjoins.close();
		if (Q[t] > Qmax.y) { Qmax.y = Q[t]; Qmax.x = t; }
		
		t++;			// increment time
	
	}
							
	
	return;
}

// ------------------------------------------------------------------------------------
void buildDeltaQMatrix(){
	// First we compute e_{i,j}, and the compute+store the a_{i} values. These will be used
	// shortly when we compute each dQ_{i,j}.
	edge    *current;
	double  eij = (double)(0.5/gparm.m);
	for(int i=1;i<gparm.maxid;i++){
		a[i]=0.0;
		if(e[i].so !=0){
			current =&e[i];
			a[i]    += eij;
			while(current->next !=NULL){
				a[i]   += eij;
				current =current->next;
			}
			Q[0] += -1.0*a[i]*a[i];		
		}	
	}
	//现在建立空的稀疏矩阵dq[]
	dq = new nodenub [gparm.maxid];						// initialize dq matrix
	for (int i=0; i<gparm.maxid; i++) {					// 
		dq[i].heap_ptr = NULL;							// no pointer in the heap at first
		if (e[i].so != 0) { dq[i].v = new vektor(2+(int)floor(gparm.m*a[i])); }
		else {			dq[i].v = NULL; }
	}
	h = new maxheap(gparm.n);						// allocate max-heap of size = number of nodes

	//b把dQ_{i,j}插入到dq[i]之中  
	double    dQ;
	tuple	dQmax;										// for heaping the row maxes
	//tuple*    itemaddress;									// stores address of item in maxheap

	for (int i=1; i<gparm.maxid; i++) {
		if (e[i].so != 0) {
			current = &e[i];								// grab first edge
			dQ      = 2.0*(eij-(a[current->so]*a[current->si]));   // compute its dQ
			dQmax.m = dQ;									// assume it is maximum so far
			dQmax.i = current->so;							// store its (row,col)
			dQmax.j = current->si;							// 
			dq[i].v->insertItem(current->si, dQ);				// insert its dQ
			while (current->next != NULL) {					// 
				current = current->next;						// step to next edge
				dQ = 2.0*(eij-(a[current->so]*a[current->si]));	// compute new dQ
				if (dQ > dQmax.m) {							// if dQ larger than current max
					dQmax.m = dQ;							//    replace it as maximum so far
					dQmax.j = current->si;					//    and store its (col)
				}
				dq[i].v->insertItem(current->si, dQ);			// insert it into vector[i]
			}
			dq[i].heap_ptr = h->insertItem(dQmax);				// store the pointer to its loc in heap
		}
	}
	delete [] elist;
	delete [] e;
	return;

}
	
// ------------------------------------------------------------------------------------
void dqSupport(){
	int    total = 0;
	int    count = 0;
	for (int i=0; i<gparm.maxid; i++) {
		if (dq[i].heap_ptr != NULL) { total += dq[i].v->returnNodecount(); count++; }
	}
	supportTot = total;
	supportAve = total/(double)count;
	return;
}

// ------------------------------------------------------------------------------------
void groupListsSetup(){
	list   *newList;
	c   = new  stub [gparm.maxid];
	for(int i=0; i<gparm.maxid; i++){
		if(e[i].so !=0){
			newList	= new list;
			newList->index = i;
			c[i].members   = newList;
			c[i].size      = 1;
			c[i].valid	   = true;
		}
	}
	return;
}


// ------------------------------------------------------------------------------------
void mergeCommunities(int i, int j){
	dpair *list, *current, *temp;
	tuple newMax;
	int t = 1;
	list    = dq[i].v->returnTreeAsList();			// get a list of items in dq[i].v
	current = list;
	cout << "stepping through the "<<dq[i].v->returnNodecount() << " elements of community " << i << endl;
	while (current!=NULL){
		if (current->x != j){
			if (dq[j].v->findItem(current->x)){
				dq[current->x].v->insertItem(j,current->y);
				dq[current->x].v->deleteItem(i);
				newMax = dq[current->x].v->returnMaxStored();
				h->updateItem(dq[current->x].heap_ptr, newMax);
				dq[j].v->insertItem(current->x,current->y);
			}else{
				double axaj = -2.0*a[current->x]*a[j];
				dq[current->x].v->insertItem(j,current->y + axaj);
				dq[current->x].v->deleteItem(i);
				newMax = dq[current->x].v->returnMaxStored();
				h->updateItem(dq[current->x].heap_ptr, newMax);
				dq[j].v->insertItem(current->x,current->y + axaj);
			}
		
		}
		temp    = current;
		current = current->next;						// move to next element
		delete temp;
		temp = NULL;
		t++;
	}
	dq[j].v->printTree();
	dq[j].v->deleteItem(i);
	newMax = dq[j].v->returnMaxStored();
	h->updateItem(dq[j].heap_ptr, newMax);

	// Again, the first thing we do is get a list of the elements of [j], so that we may
	// step through them and determine if that element constitutes an ijx-chain which
	// would require some action on our part.
	list = dq[j].v->returnTreeAsList();			// get a list of items in dq[j].v
	current = list;							// store ptr to head of list
	t       = 1;

	while (current != NULL) {
		if ((current->x != i) && (!dq[i].v->findItem(current->x))){
			double axai = -2.0*a[current->x]*a[i];
			dq[current->x].v->insertItem(j, axai);
			newMax = dq[current->x].v->returnMaxStored();
			h->updateItem(dq[current->x].heap_ptr, newMax);
			dq[j].v->insertItem(current->x, axai);			// (step 4)
			newMax = dq[j].v->returnMaxStored();			// (step 6)
			h->updateItem(dq[j].heap_ptr, newMax);
		}
		temp    = current;
		current = current->next;						// move to next element
		delete temp;
		temp = NULL;
		t++;	
	}
	a[j] += a[i];								// (step 7)
	a[i] = 0.0;
	delete dq[i].v;							// (step 8)
	dq[i].v        = NULL;						// (step 8)
	dq[i].heap_ptr = NULL;						//
	
	return;


}
// ------------------------------------------------------------------------------------
void readInputFile()
{
	int numnodes=0;
	int numlinks=0;
	int s,t,f;
	edge **last;
	edge *newedge;
	edge *current;
	bool existsFlag;
	ifstream fscan("save.txt", ios::in);
	while (fscan >> s >> f) {					// read friendship pair (s,f)
		numlinks++;							// count number of edges
		if (f < s) { t = s; s = f; f = t; }		// guarantee s < f
		if (f > numnodes) { numnodes = f; }		// track largest node index
	}
		fscan.close();
	cout << "  edgecount: ["<<numlinks<<"] total (first pass)"<<endl;
	cout << "  numnodes:  ["<<numnodes<<"] total (first pass)"<<endl;
	
	gparm.maxid = numnodes+2;					// store maximum index
	elist = new edge [2*numlinks];				// create requisite number of edges
	int ecounter = 0;							// index of next edge of elist to be used

	// Now that we know numnodes, we can allocate the space for the sparse matrix, and
	// then reparse the file, adding edges as necessary.
	cout << " allocating space for network." << endl;
	e        = new  edge [gparm.maxid];			// (unordered) sparse adjacency matrix
	last     = new edge* [gparm.maxid];			// list of pointers to the last edge in each row
	numnodes = 0;								// numnodes now counts number of actual used node ids
	numlinks = 0;	// numlinks now counts number of bi-directional edges created
	ifstream fin("save.txt",ios::in);
	while(fin >> s >> f){
		s++;f++;
			if (f < s) { t=s; s=f;f=t; }
					numlinks++;							// increment link count (preemptive)
				if (e[s].so == 0) {						// if first edge with s, add s and (s,f)
					e[s].so = s;						// 
					e[s].si = f;						// 
					last[s] = &e[s];					//    point last[s] at self
					numnodes++;						//    increment node count
				} else {								//    try to add (s,f) to s-edgelist
					current = &e[s];					// 
					existsFlag = false;					// 
					while (current != NULL) {			// check if (s,f) already in edgelist
						if (current->si==f) {			// 
							existsFlag = true;			//    link already exists
							numlinks--;				//    adjust link-count downward
							break;					// 
						}							// 
						current = current->next;			//    look at next edge
					}								// 
					if (!existsFlag) {					// if not already exists, append it
						newedge = &elist[ecounter++];		//    grab next-free-edge
						newedge -> so = s;				// 
						newedge -> si = f;				// 
						last[s] -> next = newedge;		//    append newedge to [s]'s list
						last[s]         = newedge;		//    point last[s] to newedge
					}								// 
				}									// 
		
				if (e[f].so == 0) {						// if first edge with f, add f and (f,s)
					e[f].so = f;						// 
					e[f].si = s;						// 
					last[f] = &e[f];					//    point last[s] at self
					numnodes++;						//    increment node count
				} else {								// try to add (f,s) to f-edgelist
					if (!existsFlag) {					//    if (s,f) wasn't in s-edgelist, then
						newedge = &elist[ecounter++];		//       (f,s) not in f-edgelist
						newedge -> so = f;				// 
						newedge -> si = s;				// 
						last[f] -> next = newedge;		//    append newedge to [f]'s list
						last[f]		 = newedge;		//    point last[f] to newedge
					}								// 
				}									
				existsFlag = false;						// reset existsFlag


	}
	fin.close();
	gparm.m = numlinks;
	gparm.n = numnodes;
	return ;
}

// ------------------------------------------------------------------------------------

