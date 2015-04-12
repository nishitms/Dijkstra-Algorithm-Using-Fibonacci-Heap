// Dijkstra's Algorithm Using Fibonacci heap and Simple scheme

//--------------------------------------HEADER FILES--------------------------------------------------------------------------

#include <iostream>                             // For I/O operations
#include <string>                               // For String operations
#include <stdlib.h>                             // For string to Int conversion
#include <stdio.h>                              // For operations on files and characters
#include <fstream>                              // Taking input from file
#include <math.h>                               // Using the log function
#include <time.h>                               // For clock

using namespace std;

//-------------------------------------GLOBAL VARIABLES AND MACROS--------------------------------------------------------------

#define INFINITY 99999999999                    // Gives INFINITY the following number, used for solving dijkstra's algo
long long count=0;                              // Counts the total number of edges added to the graph
long long ignored_count = 0;                    // Counts the number of edges not added and ignored because they repeated
static int count_call =0;                       // Counter used for checking the number of times a function is called
int * check_exist;                              // Array (Dynamic) used for Depth First Search
static int check_count = 0;                     // Used for checking if count is equal to numer of existing nodes(connected)
long ** short_array ;                           // Used mainly in simple scheme , for storing the source, edge , flag and cost
int num_heap_elements= 0;                       // Used to count the number of elements inserted in heap
int * check;                                    // Random array used for checking purpose

//-----------------------------------------Function Prototypes---------------------------------------------------------------

/*
This function is used for performing depth first search and flags the nodes that exist and keeps a track of unconnected nodes
Performs the search in the list from the source element
Takes the graph as a parameter.
*/
	void depth_first (struct graph * mgraph ,int * check_exist , int n);

/*
This function is used to initialize the array in which the cost of nodes is going to be stored
Used in both simple dijkstra and fibonacci dijkstra.
Takes a dynamic array and initializes for the total number of vertices
*/
	void initialize_array ( long ** & short_array, int num_vertices , int source );

/*
As the name suggests, this function is used for implementation of dijkstra's algorithm, the destimnation parameter was
just for some testing
*/
	void simple_dijkstra ( struct graph * mgraph , int num_vertices , int source , int dest);

/*
This function is used in simple dijkstra algorithm, it is used for finding the minimum source,
it then traverses through and returns the minimum source node
*/
	int find_min_array ( struct graph * mgraph, long ** & short_array, int num_vertices , int source );

/*
This function is used for updating the cost of each node , as the path progresses
*/
	void update_array ( long ** & short_array , int node_num , int node_cost, int num_vertices);

/* This function was used for checking whether we are getting the right values or not */

	int check_short_array(long ** & short_array,int num_vertices);

/*This is the most important function , what this function does is that whenever a remove min is performed, the whole heap needs to be
    consolidated, i.e the trees with same degree's get joined/consolidated */

	void consolidate(int num_vertices);

/*This function is important , as it removes the minimum element from the heap, uses the pointer to root node (minimum)*/

	void remove_min(int num_vertices);


//-----------------------------------------------------STRUCTURES AND SOME FUNCTIONS----------------------------------------------------

/* This structure was mainly used for storing various values into the fibo node, like
    the cost , the source, the edge value, the parent, child , left, right nodes etc
*/
typedef struct fibo_node
{
        int cost_f,edge_val, edge_source, degree , flag;
        fibo_node* parent , * child ,* left,* right;
        bool mark;
}node;
node *root_node=NULL;
node *n[5000];

// A structure to contain all the Vertex nodes and their cost corresponding to a vertex
struct edge_node{
        int edge_val;                             // Will store the vertex number to which the head vertex is connected
        long long edge_cost;                      // Will contain the cost from head vertex to this vertex
        struct edge_node * next_edge;             // Will point to next vertex to which the head vertex is connected
};

struct head_vertex{

        struct edge_node * head;            // The head vertex of each list
};

// Make a new node, stores the values and then returns the node
struct edge_node * new_edge_node ( int next_edge_val , long long next_edge_cost){

        struct edge_node * new_node = (struct edge_node * ) malloc(sizeof(struct edge_node));
        new_node->edge_val = next_edge_val;
        new_node->edge_cost = next_edge_cost ;
        new_node->next_edge = NULL;
        count ++ ;
        return new_node;
};

/*
 Now I have made a list i.e. 0 -> | 1 3 | -> |2 1 | , where the head vertex is 0 .
 Now the next step is to make an array of this list's, so we will have various head nodes
 and their corresponding edges with costs. This array of lists is required to make a graph.
*/

struct graph{
        int num_vertices;
        struct head_vertex * array ;
};

// Makes the graph
struct graph * makegraph ( int vertices){
        struct graph * mgraph = (struct graph *) malloc(sizeof(struct graph));
        mgraph->num_vertices = vertices ;

        mgraph->array = (struct head_vertex * ) malloc(vertices * sizeof(struct head_vertex));

        for (int i= 0 ; i< vertices ; i++)
            mgraph->array[i].head = NULL ;
        return mgraph;
};

/*
This function is used for adding an edge to the array list, it takes the graph as a parameter,
and all the required information for an edge, i.e the source to edge, the destination, the cost */

/* More Importantly, this function ignores the edges that have been previously added */

void addEdge ( struct graph * mgraph , int source , int next_edge_val , long long next_edge_cost)
{
        int node_count = 0 , check = 0 , * check_arr;

        struct edge_node * new_node = new_edge_node(next_edge_val , next_edge_cost);
        struct edge_node * check_no = mgraph->array[source].head;
        while(check_no != NULL)
        {
            node_count ++ ;
        //   cout<< "Source :"<< source <<":"<<check_no->edge_val <<"\n";
            check_no = check_no -> next_edge ;
        }

        check_no = mgraph->array[source].head;
        //cout<<" \n Running for loop : " << node_count;

        for (int i=0 ; i<node_count ; i++)
        {
            if ( (check_no->edge_val == next_edge_val) || (source == next_edge_val ))
                    {
                      //  cout<<"\n Source "<<source << "Val "<< check_no->edge_val << "=" << next_edge_val;
                        check = 1;
                        break;
                    }
                    else
                    {
                        check = 0;
                        check_no = check_no -> next_edge ;
                    }

        }


        if ( check == 0)
        {
            new_node->next_edge = mgraph->array[source].head;
            mgraph->array[source].head = new_node;

            new_node= new_edge_node(source , next_edge_cost);
            new_node->next_edge = mgraph->array[next_edge_val].head;
            mgraph->array[next_edge_val].head = new_node;
            cout<<"\n Adding Edge ";
        }
        else
        {
          //  cout<<"\n This node already exists, cannot be added , Sorry";             //Ignores the edge
            ignored_count ++ ;
            count--;
        }

}

// This function is used for showing the graph formed

void showGraph ( struct graph * mgraph)
{
        cout<<"\n       Vertices                Edge|Cost   \n";
        for(int i = 0; i< (mgraph->num_vertices) ; i++)
        {
            struct edge_node * traverse = mgraph->array[i].head;

            cout<<"\n Head vertex " << i << " is connected to :";
            while(traverse)
            {
                cout << " ->" << traverse->edge_val << "|" << traverse->edge_cost ;
                traverse = traverse -> next_edge ;
            }
            cout<<"\n";
        }
        cout<<"\n Total number of edges = "<<count;
}

void initialize (int vertices)
{
  //      min_element= new fnode;
  //      min_element=NULL;
    for (int i=0;i<vertices;i++)
    {
        n[i]= new node;
//        if(n[i]->parent != NULL)
//        cout<<" i : "<<i<< "parent :"<< n[i]->parent<<"\n";
    }
}

void insert(fibo_node *n1)
{
	if(root_node==NULL)                     //This would then be the first node entering the Fibonacci heap_node
	{
		root_node=n1;                       //This would become the root node
		n1->right=n1;                       //The right pointer will point towards itself
		n1->left=n1;                        //The left pointer will point towards itself
	}
	else
	{
		node *n2;
		n2=root_node;
		while(n2->right!=root_node)         // Reach the Last node
			n2=n2->right;                   // Increase the pointer till it reaches the last node
		n2->right=n1;                       //Now we have reached the last node , which will then point to the new node
		n1->left=n2;                        // The new node will point back to the last node
		n1->right=root_node;                // Last node will point to the root node
		root_node->left=n1;                 //  Root node will point to the last node
		if((n1->cost_f) <(root_node->cost_f ))
			root_node=n1;                   //If the new node has a less cost then it would become the root_node
	}
	num_heap_elements++;
}
/*This function is used for reducing the cost and also performing the cut and cascading cut (which is the most important operation )
of a fibonacci heap, this function decreases the key of a node and arranges itself inthe heap, either by cutting from the tree , or
remaining there if it obeys the heap property, now cascading cut is the operation that only one child can be cut from a node that
is not a root node*/
void update_Min()
{
    node * traversal;
    traversal = root_node -> right;
    int nod = traversal->cost_f, flag = 1;

    while ((traversal->cost_f != nod) || (flag)) {
        flag = 0;

        if (root_node->cost_f >= traversal->cost_f)
            root_node = traversal;

        traversal = traversal->right;
    }
}
void find_min(node *n1)
{
    node *t, *min_node;
    t=n1;
    min_node=t;
    do
    {
        if(t->cost_f < min_node->cost_f)
            min_node=t;
        t=t->right;
    }while(t!=n1);
    root_node = min_node;
    cout<<"----"<<root_node->cost_f;
}
void decrease_key(node *n1)
{
    node *parent;
    if (n1->parent != NULL)
        parent = n1->parent;
    else // inside else means n1 has no parent. It is the root of one of the trees in the fheap
    {
        update_Min();
        return;
    }
    if (n1->cost_f >= parent->cost_f) // no need for cascade cut. simply return
        return;

    //Casade Cut

    if (parent->parent == NULL) // that means 'parent' node is the root of that particular heap
    { // that means simply cut off fptr and add it to the top list


        n1->left->right = n1->right; // removing fptr from the sibling list
        n1->right->left = n1->left;

        if (n1->right != n1) // change the child pointer of the parent
            n1->parent->child = n1->right;

        else // that means there is no sibling of fptr. Make the child pointer of the parent null
            n1->parent->child = NULL;

        n1->right = n1->right; // add fptr to the top level circular list.
        root_node->right->left = n1;
        n1->left = root_node;
        root_node->right = n1;

        update_Min();
        return;
    }

    while ((parent->mark)&&(parent != root_node)) {

        node* gParent;

        if (parent->parent != NULL)
            gParent = parent->parent;

        parent->parent = NULL;

        if (gParent->child == parent)
            gParent->child = parent->right;

        parent->left->right = parent->right;
        parent->right->left = parent->left;

        parent->right = root_node->right;
        root_node->right->left = parent;
        parent->left = root_node;
        root_node->right = parent;
        parent->mark = 1;

        parent = gParent;
    }

    parent->mark =0;
    update_Min();

}
/*
void decrease_key(node *n1)
{
	node *parent=n1->parent;
	node *n2,*n3;
	if(parent!=NULL && n1->cost_f < parent->cost_f)
	{
		if(parent->child == n1)
		{
			if(n1->right==n1 && n1->right==n1)
				parent->child=NULL;
			else
				parent->child=n1->right;
		}
		else
		{
			n2=n1->left;
			n3=n1->right;
			n2->right=n3;
			n3->left=n2;
		}
		n2=root_node;
		while(n2->right!=root_node)
		{
			n2=n2->right;
		}
		n1->right=root_node;
		n1->left=n2;
		root_node->left=n1;
		n2->right=n1;
		n1->parent=NULL;
		n1->mark=0;
		if(parent->mark==0)
			parent->mark=1;
		else
			decrease_key(parent);
	}
}
*/
/*This is the function that implements dijkstra's algorithm using a fibonacci heap */
void dijkstra_fh(struct graph* fgraph, long ** short_array, int source , int num_vertices )
{
    int min_node, min_cost , cnt=1;
    initialize(num_vertices);
  //  node * n[num_vertices];


//Initializes the array in which we will keep the flags if the node has been visited or not
    initialize_array( short_array,num_vertices ,source );

    //Initialize and create all the new nodes
    for(int i=0;i<num_vertices;i++)
    {
        n[i]=new node;
        n[i]->cost_f=short_array[i][1];
        n[i]->child=NULL;
        n[i]->left=NULL;
        n[i]->right=NULL;
        n[i]->parent=NULL;
        n[i]->degree=0;
        n[i]->mark=0;
        n[i]->edge_source = i;

    }
   /* for(int i=0;i<num_vertices;i++)
    {
       //cout<<"\n i "<<i;
       insert(n[i]);
       //cout<<"\n val"<<n[i]->cost_f;
    }*/


    // Points to the source node and checks the flag and then traverses the nodes with which source is connected
    struct edge_node * traverse = fgraph->array[source].head;
    while(traverse)
    {
        if(!short_array[traverse->edge_val][2])
            {
                if  (traverse->edge_cost < (n[traverse->edge_val]->cost_f))
                {
                    n[traverse->edge_val]->cost_f = traverse->edge_cost;
                }
                traverse=traverse->next_edge;
            }
    }
     for(int i=0;i<num_vertices;i++)
    {
           if(!short_array[i][2])
           insert(n[i]);
    }

    for(int i=0;i<num_vertices-1;i++)
    {
        if(root_node!=NULL)
        {
            min_node=root_node->edge_source;
            min_cost=root_node->cost_f;
            short_array[min_node][2]=1;
		}
		else
            cout<<"\n Empty";

		    remove_min(num_vertices);
        cout<<"----------------------------------"<< i ;

        traverse = fgraph->array[min_node].head;
        while(traverse)
        {
            if(((min_cost + traverse->edge_cost) < n[traverse->edge_val]->cost_f) && short_array[traverse->edge_val][2]==0)
            {
				n[traverse->edge_val]->cost_f=min_cost + traverse->edge_cost;
				//decrease_key(n[traverse->edge_val]);
			}
			traverse=traverse->next_edge;
		}
    }

		for(int i=0;i<num_vertices;i++)
            cout<<"\n"<<n[i]->cost_f<<"/fibo/"<<" Cost from " << source <<" to "<<i<<"\n";      //Displays the output that you see
   //          n[i]->cost_f=check[i];
}



//-----------------------------------------------------MAIN------------------------------------------------------------------------------

int main ( int argcount , char * argvect[])     // argcount is the number of arguments passed in command line and argvect contains the values
{
     int ** store_values = NULL;
     string access, ran_mode = "-r", usr_mode_simple = "-s" , usr_mode_fheap= "-f";
     int n,x ,edg , dest;                              // edg is of type int edge
     float d;
     long long edge;                            // n = num of vertices , d = % edge density , x= source node number , edge = no. of edges
     clock_t start , simple_scheme , fheap_scheme;
     access = argvect[1];                       // Stores the second argument
    // cout<<"\n First : "<< argvect[0]<<"\n"<<argcount;
    // cout<<" \n Access: "<< access;

    /* For Random mode */
        if ( argcount == 5 && (access.compare(ran_mode)== 0))
        {

            n= (int) strtol (argvect[2] , NULL, 10 );
            d= (float) strtof (argvect[3] , NULL );
            x= (int) strtol (argvect[4] , NULL, 10 );
            edge =ceil ((n * (n-1)* d)/200 );
            cout<<"\n Random Mode";
            cout<<"\n n = "<<n;
            cout<<"\n d = "<<d;
            cout<<"\n x = "<<x;
            cout<<"\n edge = "<<edge;


            struct graph * fgraph = makegraph(n);

            cout<<"\n Array List";
            int cost_v,source_v,dest_v;
            for (long long i=0 ; i < edge ; i ++)
            {
                            cost_v= (rand() % 1000 ) + 1;
                            source_v= rand() % n ;
                            dest_v= rand() % n ;
                            addEdge(fgraph,source_v,dest_v,cost_v);
            }
            showGraph(fgraph);
            depth_first (fgraph ,check_exist ,n);
            cout<<"\n Graph Generated";

            start= clock();
            simple_dijkstra (fgraph , n , x , dest);
            simple_scheme = clock() - start;


            start=clock();
          //  dijkstra_fh(fgraph,short_array,check,x,n);
            fheap_scheme = clock() - start;
            cout<<"--------------------------------------------------------------------";


            cout<<" \n Max Number of edges possible  = "<<edge;
            cout<<" \n Number of non repeated edges  ="<< count/2;
            cout<<" \n Number repeated edges ignored ="<< ignored_count ;

            cout<<"\n Time : " << simple_scheme ;
           // cout<<"\n Time : " << fheap_scheme;

        }
        else
        {
           // For user input mode - simple scheme */
            if ( argcount == 3 && (access.compare(usr_mode_simple)==0))
            {
                cout<<"\n Simple Heap File Mode ";

                 ifstream  my_file;
                 my_file.open(argvect[2]);
                 if ( my_file != 0)
                 {
                    my_file>>x;     // x has value
                    my_file>>n;     // n has vertices
                    my_file>>edg;  // edge has number of edges
                    edge =edg;
                    cout<<"\n Source : " << x;
		            cout<<"\n Vertices :"<< n ;
		            cout<<"\n Edge : "<<edge;

                    struct graph * fgraph = makegraph(n);
                    int cost_v,source_v,dest_v;
                    for (long long i=0 ; i < edge ; i ++)
                    {

                            my_file>>source_v;
                            my_file>>dest_v;
                            my_file>>cost_v;
                            addEdge(fgraph,source_v,dest_v,cost_v);
                    }
                    showGraph(fgraph);
                    depth_first(fgraph ,check_exist , n);
                    dest= n - 1;
                    start= clock();
                    simple_dijkstra (fgraph , n , x , dest);
                    simple_scheme = clock() - start;
                    //cout<<"\n Time Simple: " << simple_scheme ;
                    my_file.close();
                }
                else
                    cout<<"\n Sorry ! Unable to open file \n ";

            }
            else
            {
                //For user input mode - fheap scheme
                if ( argcount == 3 && (access.compare(usr_mode_fheap)==0))
                {
                 cout<<"\n Fheap mode " ;

                 ifstream  my_file;
                 my_file.open(argvect[2]);
                 if ( my_file != 0)
                 {

                    my_file>>x;     // x has value
                    my_file>>n;     // n has vertices
                    my_file>>edg;  // edge has number of edges
                    edge =edg ;
                    cout<<"\n Source : " << x;
		            cout<<"\n Dest :"<< n ;
		            cout<<"\n Edge : "<<edge;

                    struct graph * fgraph = makegraph(n);
                    int cost_v,source_v,dest_v;
                    for (long long i=0 ; i < edge ; i ++)
                    {

                            my_file>>source_v;
                            my_file>>dest_v;
                            my_file>>cost_v;
                            addEdge(fgraph,source_v,dest_v,cost_v);
                    }
                    showGraph(fgraph);
                    depth_first (fgraph ,check_exist ,n);
                    start= clock();
                    dijkstra_fh(fgraph,short_array,x,n);
                    fheap_scheme = clock() - start;
                    //cout<<"\n Time Fibo: " << fheap_scheme ;
                    my_file.close();
                 }
                    else
                        cout<<"\n Sorry ! Unable to open file ";


                }
                else
                {
                    cout<<"\n Invalid Choice ";
                    cout<<"\n usage : program name <-r> n d x";
                    cout<<"\n usage : program name <-s> file-name ";
                    cout<<"\n usage : program name <-f> file-name ";
                }
            }
        }
        return 0;
}

//---------------------------------------------------------------FUNCTION BODY----------------------------------------------------------

void depth_first (struct graph * mgraph ,int * check_exist , int n)
{
    check_exist = new int [n];
    for(int i=0; i < n; i++)
        check_exist[i]=0;
        int i=0;

    while(i<n)
    {
        struct edge_node * check_no = mgraph->array[i].head;
        if(check_no != NULL)
           check_exist[i]=1;
        i++;
    }
    for (int i=0; i<n ;i++)
    {
        if(check_exist [i] == 1 )
        {
            check_count++;
        }

        else
            continue;
    }

}
/*Performs simple mode dijkstra's algorithm */

 void simple_dijkstra ( struct graph * mgraph , int num_vertices , int source , int dest)
{
    long min_cost;
    int min_node;
    initialize_array ( short_array,num_vertices ,source );
    /*short_array = new long * [num_vertices];
    for (int i=0;i<num_vertices;i++)
        short_array [i] = new long [3];

    for(int i=0;i<num_vertices;i++)
	{
		short_array[i][0]=i;            //Edge
		short_array[i][1]=INFINITY;     //Cost
		short_array[i][2]=0;            //Status
	}

	short_array[source][1]=0;           //Cost
	short_array[source][2]=1;           //Status */

    struct edge_node * traverse = mgraph->array[source].head;
   // traverse=traverse->next_edge;
	while(traverse!=NULL)
	{
		if(!short_array[traverse->edge_val][2])
		{
			if(traverse->edge_cost < short_array[traverse->edge_val][1])
			{
				short_array[traverse->edge_val][1]= traverse->edge_cost;
			}
		}
		traverse=traverse->next_edge;
	}

    for(int i=0;i<num_vertices-1;i++)
	{
		min_node=-1;
		min_cost=INFINITY;

		for(int j=0;j<num_vertices;j++)
		{
			if(short_array[j][1]<=min_cost && short_array[j][2] == 0)
			{
				min_cost=short_array[j][1];
				min_node=j;
			}
		}


        short_array[min_node][2]=1;           //Status

        traverse = mgraph->array[min_node].head;
     //   traverse=traverse->next_edge;

		while(traverse)
		{
			if((min_cost+traverse->edge_cost)< short_array[traverse->edge_val][1] && short_array[traverse->edge_val][2] == 0 )
			{
				short_array[traverse->edge_val][1]=min_cost+traverse->edge_cost;
			}
			traverse=traverse->next_edge;
		}

	}
   /* int main_source=source;
    initialize_array ( short_array,num_vertices ,source );
  //source = find_min_array ( short_array,num_vertices ,source );
    for(int i=0; i< num_vertices ; i++)
    {
        struct edge_node * traverse = mgraph->array[source].head;
        while(traverse)
        {
            if  (( ( traverse->edge_cost ) + short_array[source][1] )  < (short_array[traverse->edge_val][1]) && (short_array[traverse->edge_val][2]==0) )
                update_array (short_array , traverse->edge_val , (( traverse->edge_cost ) + short_array[source][1] ), num_vertices);
            traverse=traverse->next_edge;
          //  new_source = find_min_array ( short_array,num_vertices ,source );
        }
        source = find_min_array ( mgraph ,short_array, num_vertices ,source );

       // if(source != dest)
       //     cout<<"source "<< source;
        //  cout<<"\n count"<< count_call;
      //if(j == (count/2)+1 )
       //     break;
     //   else
     //       j++;
          if(check_short_array(short_array,num_vertices) != 1 )
            break;
          else
            continue;

          cout<<"\n Depth First Performed, Existing elements found";*/


    for(int i=0; i< num_vertices ; i++)
    {
     //   check[i]=short_array[i][i];
    //cout<< "\n Source : "<< find_min_array ( short_array,num_vertices ,source );
    //for(int i = 0 ; i <num_vertices ; i++)
    // cout << "\n S :" << short_array[i][0]<< "\n C :" << short_array[i][1]<< "\n F :" << short_array[i][2]<<"\n";
        cout<<"\n"<<short_array[i][1]<< " // "<< "  Cost from " <<source<< " to " <<short_array[i][0]<<"\n";

    }
}

 void initialize_array ( long ** & short_array, int num_vertices , int source )
{
            short_array = new long * [num_vertices];
            for (int i=0;i<num_vertices;i++)
                short_array [i] = new long [3];


            for (int i =0 ; i<num_vertices ; i++)
            {
                if(i==source)
                {
                    short_array [i][0] = i ;
                    short_array [i][1] = 0;
                    short_array [i][2] = 1;
                }
                else
                {
                        short_array [i][0] = i ;
                        short_array [i][1] = INFINITY;
                        short_array [i][2] = 0;
                }

            }
            check =new int [num_vertices];

}
int find_min_array ( struct graph * mgraph,long ** & short_array, int num_vertices , int source )
{
            int source_val;

            static long min_cost, min_cost_source;

            if(count_call == 0)
            {
                    source_val= source;
         //           cout<< "\n"<<min_cost<< " // " << " Cost from " << source_val << " to " << source_val ;
                    count_call++;
            }
            else
                 count_call++;

      //      for(int i = 0 ; i <num_vertices ; i++)
      //          cout << "\n S :" << short_array[i][0]<< "\n C :" << short_array[i][1]<< "\n F :" << short_array[i][2]<<"\n";

 struct edge_node * traverse = mgraph->array[source].head;



            for (int i =0 ; i<num_vertices -1 ; i++)
            {
                if ( (short_array [i][2] == 0 ))
                {
                    while (traverse)
                                {
                                    if( i == (traverse->edge_val))
                                    {
                                     //  cout<<"\n Foun Min "<<i;
                                       min_cost=short_array [i][1];
                                       break;
                                    }
                                    else
                                    {
                                        traverse = traverse->next_edge;
                                    }

                                }
                }
                else
                    continue;
            }


            traverse = mgraph->array[source].head;
      //  cout<<"\n Source"<<source;
            for (int i =0 ; i<num_vertices ; i++)
            {
              //  if ( (short_array [i][2] == 1 ))
               // {
//cout<<"\n IN";
// cout<< short_array[i][0] <<"\n";
                            int cnt=0,temp;
                               while (traverse)
                                {
                                 if(cnt==0 && (short_array [traverse->edge_val][2]==0) )
                                 {

                                        min_cost = short_array [traverse->edge_val][1];
                                        min_cost_source = traverse->edge_val;
                                        cnt++;
                                 }
                                 else
                                 {
                                     temp=short_array [traverse->edge_val][1] ;
                                     if(temp<min_cost && (short_array [traverse->edge_val][2]==0))
                                     {
                                        min_cost = short_array [traverse->edge_val][1];
                                        min_cost_source = traverse->edge_val;
                                     }
                                 }

                                  //  }
                                 //   else
                                        traverse = traverse->next_edge;
                                }
                }



            short_array[ min_cost_source ] [2] = 1;
            source=min_cost_source;
            return min_cost_source;
}
void update_array ( long ** & short_array , int node_num , int node_cost, int num_vertices)
{

	    short_array [node_num][1] = node_cost;
}
int check_short_array(long ** & short_array, int num_vertices)
{
            for (int i =0 ; i<num_vertices ; i++)
            {
                if ( (short_array [i][2] == 1 ))
                    continue;
                else
                    return 1;
            }
    return 0;
}

void consolidate(int num_vertices)
{
    int d = (log(num_heap_elements)/(log(2)));
    //int n=num_vertices;
	node *f_line[num_vertices],*n1,*n2;
//	for(int i=0;i<(log(num_heap_elements)/(log(2)));i++)
	for(int i=0;i<d;i++)
	{
		f_line[i]=NULL;
	}
	n1=root_node;
	node *n3;
	node *n4;
	node *n5;
	node *n6;
	for(int i=0;i<d;i++)
	{
		if(f_line[n1->degree]==NULL)
        {
            cout<<"\n //1 ";
			f_line[n1->degree]=n1;
			n1=n1->right;
		}
		else
		{
			n2=f_line[n1->degree];
			if(n1==n2)
			{
			    		cout<<"\n //2 ";
				n1=n1->right;
				continue;
			}
			if(n1->cost_f <= n2->cost_f)
			{
				if(n2==root_node)
					root_node=n1;
				n5=n2->left;
				n6=n2->right;
				n5->right=n6;
				n6->left=n5;
				if(n1->child!=NULL)
				{
					n3=n1->child;
					n4=n3->right;
					n3->right=n2;
					n2->left=n3;
					n2->right=n4;
					n4->left=n2;
					n2->parent=n1;
							cout<<"\n //3 ";
				}
				else
				{
					n1->child=n2;
					n2->right=n2;
					n2->left=n2;
					n2->parent=n1;
				}
				f_line[n1->degree]=NULL;
				n1->degree+=1;
			}
			else
			{
			    		cout<<"\n //4 ";
				n5=n1;
				n1=n2;
				n2=n5;
				n5=n2->left;
				n6=n2->right;
				n5->right=n6;
				n6->left=n5;
				if(n1->child!=NULL)
				{
					n3=n1->child;
					n4=n3->right;
					n3->right=n2;
					n2->left=n3;
					n2->right=n4;
					n4->left=n2;
					n2->parent=n1;
				}
				else
				{
					n1->child=n2;
					n2->right=n2;
					n2->left=n2;
					n2->parent=n1;
				}
				f_line[n1->degree]=NULL;
				n1->degree+=1;
						cout<<"\n //5 ";
			}
		}
	}
}

void remove_min(int num_vertices)
{
	node *n1,*n2,*n3,*n4, *n5;
	if(root_node->right == root_node && root_node->left == root_node)
	{
		if(root_node->child==NULL)
			root_node=NULL;
		else
		{
			n1=root_node->child;
			n2=n1;
			do
			{
			    n2->parent=NULL;
				n2=n2->right;
				num_heap_elements++;
			}while(n2!=n1);

			update_Min();
			find_min(n1);
			//root_node=n2;

			consolidate(num_vertices);
		}
	}
	else
	{
	    cout<<"\n Here ";
	    cout<<root_node->parent;
		n1=root_node->left;
		n2=root_node->right;
		if(root_node->child==NULL)
		{
		   // cout<<"\n*";
			n1->right=n2;
			n2->left=n1;
            update_Min();
			find_min(n2);
		}
		else
		{
		    cout<<"\n"<<root_node->parent;
		    cout<<"\n $ ";
			n3=root_node->child;
			n4=n3;
			do
            {
                cout<<"\n }";
                n4->parent= NULL;
                n4=n4->right;
            }while(n4->right!=n3);
            n5=n3->left;
			n1->right=n3;
			n3->left=n1;
			n2->left=n5;
			n5->right=n2;
			cout<<":";
			find_min(n3);
			update_Min();
		}
            consolidate(num_vertices);
	}
}





















