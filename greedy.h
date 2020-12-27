#include<bits/stdc++.h>
using namespace std;
namespace planner_info
{


int obs_ts=INT_MAX;
long long plan_till=100;
//long long horizon_len=60;
int cy_covered=-1;
int sz=10,qtot=3;
vector<int> dest;
int time=0;
unsigned long long counter=0;
bool add_path=false;

struct point
{
    int x,y;
    point(int a,int b)
    {
        x=a;
        y=b;
    }
    point()
    {

    }
};


int nrow,ncol;
int dir=5,total_system_states=0,comp_red=0,tot_comp=0,grid_sz;
int nb[][2]={{0,0},{0,-1},{-1,0},{1,0},{0,1}};

unordered_map<long long int,int> grid;
vector< vector< vector< vector<int> > > > trans(qtot,vector< vector< vector<int> > >(qtot) );
vector< vector< vector< vector<int> > > > negtrans(qtot,vector< vector< vector<int> > >(qtot) );
vector<int> neg_trans_to_neighbour;
vector<vector<int> > pos_system_state;
vector<unordered_map<int,float> > adj;
vector<unordered_map<int,int> > updated;

unordered_map<string,vector<int> > prop;
vector<vector<int> > qsystate;
vector< vector< vector<int> > > qsytrans(qtot,vector< vector<int> >(qtot));
vector<vector< vector<int> > > prop_sys_states(1000);
//int k=0;

unsigned long long int key(vector<int> &v)
{
    long long int val;
    if(v.size()==3) 
    val = v[0]*sz*sz+v[1]*sz+v[2];
    else
    val = v[0]*sz+v[1];

    return val;
}


long long int key_r(vector<int> &v)
{
    long long int val;
    if(v.size()==3) 
    val = v[2]*sz*sz+v[1]*sz+v[0];
    else
    val = v[1]*sz+v[0];

    return val;
}


long long int key_v(int a,int b,int c)
{
 return(a*sz*sz+b*sz+c);
}


string conv_vec_to_string(vector<int> &vec)
{
    string s="";
    for(int i=0;i<vec.size();i++)
    {
        stringstream ss;
        ss << vec[i];
        s= s+ss.str();
        s+=",";
    }
    return s;
}


void printvec(vector<int> &vec)
{
    for(int i=0;i<vec.size();i++)
        cout<<vec[i];
}

bool valid(int y,int x)
{
    if(x>=0 && x<ncol && y>=0 && y<nrow && grid[y*sz+x]!=-1)
      return 1;
      return 0;
}


struct cycle_node
{
  int x,y,state;
  vector<cycle_node*> non_adj;
  float f;
  cycle_node* par;
  bool operator==(const cycle_node* a)const
  {
    if(this->x==a->x && this->y==a->y && this->state==a->state)
      return true;
      return false;
  }

  bool operator==(const cycle_node& a)const
  {
    if(x==a.x && y==a.y && state==a.state)return true;
    return false;
  }

  bool operator!=(cycle_node* a)
  {
    if(this->x!=a->x || this->y!=a->y || this->state!=a->state)return true;
    return false;
  }

   bool operator<(cycle_node* a)
  {
        //cout<<"comparing "<<a->coord.x<<","<<a->coord.y<<" and "<<b->coord.x<<","<<b->coord.y<<endl;
        if(f < a->f){return true;}
        else if(f > a->f){return false;}
        else if(x==a->x && y==a->y && state==a->state)return(false);
        else return(true);
  }

};



cycle_node* new_cycle_node(int x,int y,int s)
{
    cycle_node* nd = new cycle_node;
    nd->x=x;
    nd->y=y;
    nd->state = s;
    nd->f=FLT_MAX;
    nd->par=NULL;
    return nd;
}


string state_to_str(cycle_node* nd)
{
string str=to_string(nd->x)+","+to_string(nd->y)+","+to_string(nd->state);
return(str);
}


struct key_comparator_cycle
{
    bool operator()(cycle_node* a , cycle_node* b)
    {
        //cout<<"comparing "<<a->coord.x<<","<<a->coord.y<<" and "<<b->coord.x<<","<<b->coord.y<<endl;
        if(a->f < b->f){return true;}
        else if(a->f > b->f){return false;}
        else if(a->x==b->x && a->y==b->y && a->state==b->state)return(false);
        else return(true);
    }
};

typedef map<cycle_node*,long long int,key_comparator_cycle > cycle_oq;
cycle_oq cycle_queue;

struct node
{
   int x,y,state;
   int search;
   float g,h,f;
   node* par;
   node* next;


 bool operator==(const node* a)const
 {
    if(x==a->x && y==a->y && state==a->state){return true;}
    return false;
 }

 bool operator==(const node& a)const
 {
   if(x==a.x && y==a.y && state==a.state)return true;
   return false;
 }

 bool operator!=(node* a)
 {
   if(x!=a->x || y!=a->y || state!=a->state)return true;
   return false;
 }

 bool operator<(const node* a)const
 {
   return(f < a->f);
 }

};
node* start;
node* goal;


float get_initial_heuristic(node* a,node* b)
{
 return(abs(a->x-b->x)+abs(a->y-b->y));
}


int step_heuristic(node* a,node* b)
{
 float max,min;
 max=(abs(a->x-b->x)>abs(a->y-b->y))?abs(a->x-b->x):abs(a->y-b->y);
 min=(abs(a->x-b->x)<abs(a->y-b->y))?abs(a->x-b->x):abs(a->y-b->y);
 
 return(min+(max-min));
}

node* new_node(int x,int y,int state)
{

   node* temp=new node;
   temp->x=x;
   temp->y=y;
   temp->state=state;
   temp->search=0;
   temp->h=FLT_MAX;
   temp->next=NULL; 
   temp->g=FLT_MAX;
   temp->f=FLT_MAX;

  return(temp);
}




struct node_comparator
{
    bool operator()(const node* a ,const node* b)const
    {
        //cout<<"comparing "<<a->x<<","<<a->y<<" and "<<b->x<<","<<b->y;
        if(a->f < b->f)return(true);
        if(a->f > b->f)return(false);
        if(a->f==b->f && a->h<b->h)return true;
        if(a->f==b->f && a->h>b->h)return false;;        
        if(a->x==b->x && a->y==b->y && a->state==b->state)return(false); 
        return true;    
   
    }
};


struct node_pq_comparator
{
    bool operator()(const node* a ,const node* b)const
    {
        //cout<<"comparing "<<a->x<<","<<a->y<<" and "<<b->x<<","<<b->y;
        if(a->f > b->f)return(true);
        if(a->f < b->f)return(false);
        if(a->x==b->x && a->y==b->y && a->state==b->state)return(false); 

   
    }
};

struct h_comp
{
    bool operator()(node* a , node* b)
    {
        //cout<<"comparing "<<a->x<<","<<a->y<<" and "<<b->x<<","<<b->y<<endl;

        if(a->h < b->h)return(true);
        if(a->h > b->h)return(false);
        if(a->x==b->x && a->y==b->y && a->state==b->state)return(false); 
   
    }
};


typedef map<node*,long long int,node_comparator> node_oq;
typedef map<node*,long long int,h_comp > temp_oq;

unordered_map<unsigned long long,node*> pv_created;
unordered_map<node*,unsigned long long> dy_obs;
list<node*> dec_cost;
list<node*>closed;

vector<node*> neigh;

float cal_heuristic_cost(point a,point b)
{
    float mx = max(abs(a.y-b.y),abs(a.x-b.x));
    float mn = min(abs(a.y-b.y),abs(a.x-b.x));
    return mn*1.5+(mx-mn);

}


struct dstar_inf
{
  int x,y,state;
  bool operator==(const dstar_inf& a)const
  {
    if(x==a.x && y==a.y && state==a.state)return true;
    return false;
  }
};



//Planner object

struct planner_object
{

unsigned long long plan_dur;
vector<pair<unsigned long long,cycle_node*>> change;
vector<pair<unsigned long long,cycle_node*>> path;

};


cycle_node* s_start;
cycle_node* s_goal;
cycle_node* s_last;
cycle_node* new_start;
int km=0;

unordered_map<long long int,cycle_node*>cy_created;
unordered_map<long long int,bool> cy_visited;
unordered_map<long long int,int> obs1;
unordered_map<long long int,int> obs2;


dstar_inf temp1,temp2;
vector<dstar_inf>dest_coords;
vector<cycle_node*>dest_state;
vector<pair<int,cycle_node*>> suf_cycle;
vector<pair<int,cycle_node*>> pref_path;

class MyHashFunction { 
public:  
    size_t operator()(const dstar_inf& a) const
    { 
        return (a.y*sz*sz + a.x*sz + a.state);
    } 
};


struct plan_node
{
  
  unsigned long long from;
  float cy_len;
  vector<cycle_node*>plan;

};


plan_node* new_plan_node(unsigned long long from,float cy_len)
{
 plan_node* temp_node=new plan_node;
 temp_node->from=from;
 temp_node->cy_len=cy_len;
 return(temp_node);

}


struct comp_by_time
{
  bool operator()(const std::string& a, const std::string& b) const
  {
          unsigned a_str = a.find_last_of("_");
          unsigned a_time=stoi(a.substr(a_str+1,a.size()-1));
          
          unsigned b_str = b.find_last_of("_");
          unsigned b_time=stoi(b.substr(b_str+1,b.size()-1));  
          
          if(a_time<=b_time)
          return(true);
          else
          return(false);
  
  }

};


unordered_map<dstar_inf,cycle_node*,MyHashFunction>graph_root;
map<cycle_node*,pair<float,vector<pair<int,cycle_node*>>>> all_suf_cycle;
vector<unsigned long long> num_cy;

cycle_node* current_best;
cycle_node* local_best;
float best_cost=FLT_MAX;
bool is_pref=false;


long long int global_ts=0;
map<unsigned long long,vector<pair<int,int> > >dynamic_obstacles;
map<pair<int,int>,long long> grid_dy;
unordered_map<cycle_node*,vector<plan_node>> all_cy_cost;
vector<pair<float,cycle_node*>>suf_dy_path;

//vector<node*,vector<pair<int,vector< pair<unsigned long long,node*>>>>>plan;

map<pair<unsigned long long,unsigned long long>,float>static_edge_len;
map<pair<unsigned long long,unsigned long long>,vector<pair<unsigned long long,float>>>dyn_edge_len; 

map< int,cycle_node*> pwh_plan;
map< int,cycle_node*> temp_plan;
map< int,cycle_node*> pwh_pref;
map< int,cycle_node*> pwh_suf;
vector<cycle_node*> edge_dy_path;

//vector<pair<unsigned int,cycle_node*>> pwh_suffix;
cycle_node* prev_f=NULL;
cycle_node* cur_f=NULL;
int cycle_count=0;



/* GUROBI Data Structures */
unordered_map<string,cycle_node*> str_to_cy;
map<int,unordered_set<string>> obj;
map<string,vector<string>,comp_by_time> cons;
unordered_set<string> last_cons;
unordered_set<string> gen_cons;
map<int,unordered_set<string>> all_var;
unordered_map<string,unordered_set<string>> int_pos;
unordered_map<string,unordered_set<string>> int_last_pos;
unordered_set<string> ind_var;
map<pair<int,int>,int> cmap;

}
