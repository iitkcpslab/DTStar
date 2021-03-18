#include"greedy.h"
#include<unistd.h>
#define mx_steps 50000
#define INF 100000 
#define ull unsigned long long

using namespace std;
using namespace planner_info;

fstream reader2;
cycle_node* origin=new_cycle_node(0,0,0);
cycle_node* new_origin;


pair<int,int> create_automata_trans_table(char ltl_query[])
{
	
    string table_text="";
    const int BUFSIZE=128;
    char buf[BUFSIZE];
    FILE *fp;
    ofstream outfile;
    
    outfile.open("trans.dat");
    
    int automata_transition_cnt=0;
    pair<int,int> automata_info;

    //writing the automaton transition table returned by LTL2BA converter to trans.dat
    if ((fp = popen(ltl_query, "r")) == NULL)
    {
        printf("Error opening pipe!\n");
        return automata_info;
    }
    
    //copying the transition table to s
    while (fgets(buf, BUFSIZE, fp) != NULL) 
    {
        table_text = table_text+buf;
    }

    if(pclose(fp))  
    {
        printf("Command not found or exited with error status\n");
        return automata_info;
    }

    //maps each automaton state to a  number
    unordered_map<string,int> mp;
    unordered_map<string,int> state;
    //stores all automaton states
    vector<string> state_list;

    string tmp;
    int i=0;
    while(table_text[i]!='T' && table_text[i]!='a')
        i++;

    int state_no=0;
    
    vector<string> table;

    while(i<table_text.size())
    {
       tmp="";
       while(i<table_text.size() && table_text[i]!=':')
        {
           tmp+=table_text[i];
           i++;
        }
        if(tmp[0]=='}')
            break;
       if(state.find(tmp)==state.end())
        {
          state[tmp] = state_no;
          state_list.push_back(tmp);
          state_no++;
        }

        while(i<table_text.size() && tmp!="fi;")
        {
        tmp = "";
            while(i<table_text.size() && table_text[i]!='\t' && table_text[i]!=' ' && table_text[i]!='\n')
            {
                if(table_text[i]!='(' && table_text[i]!=')')
                tmp+=table_text[i];
                i++;
            }
            if(tmp!=" " && tmp!="" && tmp!="if" && tmp!="::" && tmp!="->" && tmp!="goto" &&  tmp!="fi;")
            table.push_back(tmp);
            i++;
        }

    }

    for(i=0;i<state_list.size();i++)
    {
        mp[state_list[i]]=i;
        if(state_list[i][0]=='a' && state_list[i][1]=='c' && state_list[i][2]=='c')
            dest.push_back(i);
    }

   vector<vector<int> > transition_condn_list;
   vector<int> transition_condn;
   int x,y;
   int automaton_state=0;
   i=1;

   while(i<table.size())
   {
        if(table[i][0]==':')
        {
            automaton_state++;
        }
        else if(table[i]=="||")
        {
           transition_condn_list.push_back(transition_condn);
           transition_condn = vector<int>(0);
       }
       else if(table[i]=="&&")
       {
       }
       else if(mp.find(table[i])!=mp.end())
       {
           transition_condn_list.push_back(transition_condn);
           for(x=0;x<transition_condn_list.size();x++)
           {
               automata_transition_cnt++;
               outfile<<automaton_state<<" "<<mp[table[i]]<<" "<<transition_condn_list[x].size()<<" ";
               for(y=0;y<transition_condn_list[x].size();y++)
               {
                   outfile<<transition_condn_list[x][y]<<" ";
               }
               outfile<<"\n";
           }
           transition_condn_list = vector<vector<int> >(0);
           transition_condn = vector<int>(0);
       }
       else if(table[i]=="1")
       {
           transition_condn.push_back(0);
       }
       else if(table[i][0]=='!')
       {
           int v=0;
           for(x=2;x<table[i].size();x++)
               v = v*10+(table[i][x]-'0');
           transition_condn.push_back(-1*v);
       }
       else if(table[i][0]=='p')
       {
           int v=0;
           for(x=1;x<table[i].size();x++)
               v = v*10+(table[i][x]-'0');
           transition_condn.push_back(v);
       }

       else if(table[i]=="skip")      
        {automata_transition_cnt++; outfile<<automaton_state<<" "<<automaton_state<<" 1 0\n";}

       i++;
    }
    outfile.close();

    automata_info.first = state_list.size();
    automata_info.second = automata_transition_cnt;
    return automata_info;
}

void initializegrid(char *grid_info_file)
{
    string s;
    int i,j,k,no_of_trans,literal,num_obstacles;
    char buff[1000];
    FILE *ft;
    ft = fopen("query.dat", "r");
    /*reading LTL query*/
    fgets(buff, 1000, (FILE*)ft);

    pair<int,int> automata_info;
    /**reading the automata***/
    automata_info = create_automata_trans_table(buff);
    no_of_trans = automata_info.second;
    qtot = automata_info.first;

    /**vector for storing the automata transition table**/
    trans = vector< vector< vector< vector<int> > > >(qtot,vector< vector< vector<int> > >(qtot) );
    /**storing automata transitions on conjunction of negative literals**/
    negtrans = vector< vector< vector< vector<int> > > >(qtot,vector< vector< vector<int> > >(qtot) );
    neg_trans_to_neighbour = vector<int>(qtot);
    ifstream ifile;
    ifile.open("trans.dat");

    vector<int> automata_states;
    vector<int> transition_condn;
    
    for(i=0;i<no_of_trans;i++)
    {
        transition_condn = vector<int>(0);
        automata_states = vector<int>(2);
        ifile>>automata_states[0];
        ifile>>automata_states[1];
        int trans_condn_len=0,neg_literals=0,strict_neg_literals=0;
        ifile>>trans_condn_len;
        for(j=0;j<trans_condn_len;j++)
        {
            ifile>>literal;
            transition_condn.push_back(literal);
            if(literal <= 0)
                neg_literals++;
            if(literal < 0)
                strict_neg_literals++;
        }
        trans[automata_states[0]][automata_states[1]].push_back(transition_condn);
        /**if transition condition is conjunction of negative literals**/ 
        if(neg_literals==trans_condn_len)
        {
            if(automata_states[0]!=automata_states[1] && strict_neg_literals==trans_condn_len)
            neg_trans_to_neighbour[automata_states[0]]=1;
            negtrans[automata_states[0]][automata_states[1]].push_back(transition_condn);
        }


    }

    ifstream grid_file;
    grid_file.open(grid_info_file);
    int x,y;
    grid_file>>nrow;
    grid_file>>ncol;
    sz = max(nrow,ncol);

    grid_file>>num_obstacles;
    grid_sz = sz;
    sz = max(sz,qtot);
    /*storing obstacle coordinates*/
    for(i=0;i<num_obstacles;i++)
    {
        grid_file>>x;
        grid_file>>y;
        grid[y*sz+x]=-1;
    }
    
    int num_pos_system_states;
    vector<int> grid_state;
    grid_file>>num_pos_system_states;
    
    /** reading coordinates of states and proposition true at it **/
    for(i=0;i<num_pos_system_states;i++)
    {
        grid_state = vector<int>(2);
        grid_file>>grid_state[0];
        grid_file>>grid_state[1];
        grid_file>>literal;
        s = conv_vec_to_string(grid_state);
        /**pushing the proposition true at state s**/
        prop[s].push_back(literal);
        /**storing cells associated with a proposition**/
        vector<int> ivec(2);
        ivec[0]=grid_state[1];
        ivec[1]=grid_state[0];
        prop_sys_states[literal].push_back(ivec);
        /**storing the list of states with a proposition true at it**/
        if(prop[s].size()==1)
            pos_system_state.push_back(grid_state);
    }
    grid_file.close();
}



void calsystates()
{
    int i,j,z,cnt=0;
    vector<int> v;
    string s;
    //TS states satifying an incoming transition to an automaton state
    qsystate = vector< vector<int> >(qtot);
    
    //TS states satifying an automaton transition
    qsytrans = vector< vector< vector<int> > >(qtot,vector< vector<int> >(qtot) );

    for(int nstate=0;nstate<qtot;nstate++)
    {
        for(j=0;j<qtot;j++)
        {
            for(z=0;z<trans[nstate][j].size();z++)
            {
                vector<int> transition_req = trans[nstate][j][z];
                if(transition_req.size()==0)
                    continue;

                for(i=0;i<pos_system_state.size();i++)
                {

                    point pt;
                    pt.y =  pos_system_state[i][1];
                    pt.x =  pos_system_state[i][0];

                    v = vector<int>(0);
                    v.push_back(pt.x);
                    v.push_back(pt.y);

                    vector<int> prop_cur;
                    s = conv_vec_to_string(v);
                    if(prop.find(s)==prop.end())
                    { 
                        prop_cur.push_back(0);
                    }
                    else
                    {
                        prop_cur = prop[s];
                    }


                    int pos=0;
                    cnt=0;
                    for(int k=0;k<transition_req.size();k++)
                    {
                	    if(transition_req[k]>=0)
                            for(int l=0;l<prop_cur.size();l++)
                            {
                                if(transition_req[k]==prop_cur[l])
                                    {cnt++; pos=1; break;}
                            }
                            else
                            {
                                int tmp=1;
                                for(int l=0;l<prop_cur.size();l++)
                                {
                                    if(-1*transition_req[k]==prop_cur[l])
                                        {tmp=0; break;}
                                }
                                cnt+=tmp;
                                if(tmp==0)
                                    break;
                            }
                    }
                   if(pos && cnt==transition_req.size())
                    {
                        qsystate[j].push_back(i);
                        qsytrans[nstate][j].push_back(i);
                    }

                }
            }
        }
    }
}



//**adding/removing obstacles to/from the grid**/
void copy_map(unordered_map<long long int,int> &obs,int v)
{
    unordered_map<long long int,int>::iterator it;
    for(it = obs.begin();it!=obs.end();it++)
    {
        if(it->second)
        grid[it->first]=v;
    }
}




float calc_cost(cycle_node* source,cycle_node* dest)
{
   

   if(source->x==dest->x && source->y==dest->y)
   return(0);
   
    unordered_map<long long int,node* > init_vertex;
    unordered_map<long long int,int > vis;

    vector<int> node_inf={source->y,source->x,source->state};

    node* start_nd = new node;
    start_nd->x=source->x;
    start_nd->y=source->y;
    start_nd->state=source->state;
    start_nd->g=0;  
    
    node* goal_nd= new node;
    goal_nd->x=dest->x;
    goal_nd->y=dest->y;
    goal_nd->state=dest->state;
    goal_nd->h=0;
    goal_nd->g=FLT_MAX;
    
    start_nd->h=get_initial_heuristic(start_nd,goal_nd);
    start_nd->f=start_nd->g+start_nd->h;


    node_oq qopen;
    vector<node*> closed;
    qopen.insert(make_pair(start_nd,key(node_inf)));
    vector<int> src_vertex={start_nd->y,start_nd->x,start_nd->state};  
    vector<int> dest_vertex={goal_nd->y,goal_nd->x,goal_nd->state};        
    
    init_vertex[key(src_vertex)]=start_nd;
    init_vertex[key(dest_vertex)]=goal_nd;
    
    while(!qopen.empty())
    {
        node* cur_nd = qopen.begin()->first;
        qopen.erase(qopen.begin());
        int current_automaton_state = cur_nd->state;
        vector<int> cur_node_info{cur_nd->y,cur_nd->x,current_automaton_state};
        
        if(cur_nd->x==dest->x && cur_nd->y==dest->y)
        {
            for(auto it:closed)
            delete it;   
            return(cur_nd->f);
        }

        long long int cur_node_k = key(cur_node_info);
        
        if(vis.find(cur_node_k)!=vis.end())
            continue;
        else
            vis[cur_node_k]=1;

        vector<int> neigh_state;

        point nbh;
        for(int i=0;i<dir;i++)
        {

            if(!valid(cur_nd->y+nb[i][0],cur_nd->x+nb[i][1]) && ((cur_nd->y+nb[i][0])!=dest->y || (cur_nd->x+nb[i][1]!=dest->x)))
                continue;

            nbh.y =  cur_nd->y+nb[i][0];
            nbh.x =  cur_nd->x+nb[i][1];
            neigh_state = vector<int>(0);
            neigh_state.push_back(nbh.y);
            neigh_state.push_back(nbh.x);
            neigh_state.push_back(current_automaton_state);

            node* oldtmp;
              
            if(init_vertex.find(key(neigh_state))==init_vertex.end())
            {   
                node* neigh_node = new_node(nbh.x,nbh.y,current_automaton_state);           
                neigh_node->g = cur_nd->g+1;
                neigh_node->h = get_initial_heuristic(neigh_node,goal_nd);
          
                neigh_node->f = neigh_node->g+neigh_node->h;
                neigh_node->par = cur_nd;
                init_vertex.insert(make_pair(key(neigh_state),neigh_node) );
                qopen.insert(make_pair(neigh_node,key(neigh_state)));

            }
            else
            {

                oldtmp = init_vertex[key(neigh_state)];
                if(oldtmp->g > (cur_nd->g+1))
                {
                    if(qopen.find(oldtmp)!=qopen.end())
                    qopen.erase(qopen.find(oldtmp));

                    oldtmp->g = cur_nd->g+1;
                    oldtmp->f = oldtmp->g+oldtmp->h;
                    oldtmp->par = cur_nd;

                    qopen.insert(make_pair(oldtmp,key(neigh_state)));
                    
                }

            }

        }

    closed.push_back(cur_nd);
    
   }
return(FLT_MAX);

}



void find_non_adj(queue<cycle_node*> &qopen,cycle_node* cur_nd)
{
    int zeroloop=0,negloop=0;
    int current_automaton_state = cur_nd->state;
    vector<int>cur_node_info;
    point nbh;
    point cur_pt={cur_nd->x,cur_nd->y};

    if(trans[current_automaton_state][current_automaton_state].size()==1 && trans[current_automaton_state][current_automaton_state][0].size()==1 && trans[current_automaton_state][current_automaton_state][0][0]==0)
        zeroloop=1;

    if(negtrans[current_automaton_state][current_automaton_state].size()>0) 
        negloop=1;


	    for(int j=0;j<qtot;j++)
	    {
	        if(trans[current_automaton_state][j].size()==1 && trans[current_automaton_state][j][0].size()==1 && trans[current_automaton_state][j][0][0]==0)
	        {   
	            cycle_node* neigh_nd=new_cycle_node(cur_nd->x,cur_nd->y,j);
		    cur_nd->non_adj.push_back(neigh_nd);
	        }
	        else
	        {

	            for(int z=0;z<qsytrans[current_automaton_state][j].size();z++)
	            {
                    int neigh_state = qsytrans[current_automaton_state][j][z];
                    vector<int> neigh_node_info={pos_system_state[neigh_state][0],pos_system_state[neigh_state][1],j};
                    nbh.x = neigh_node_info[0];
                    nbh.y = neigh_node_info[1];
                    if(!valid(nbh.y,nbh.x) || (!negloop && !zeroloop && cal_heuristic_cost(cur_pt,nbh)>=2) )
                         continue;

                      vector<int> node_inf={nbh.y,nbh.x,j};
                      if(cy_created.find(key(node_inf))==cy_created.end())
                       { 

                         cycle_node* neigh_nd = new_cycle_node(nbh.x,nbh.y,j);
                         cy_created.insert(make_pair(key(node_inf),neigh_nd));
                         qopen.push(neigh_nd);
		         cur_nd->non_adj.push_back(neigh_nd);
         
	               }
                      else
                        {  
                           cur_nd->non_adj.push_back(cy_created[key(node_inf)]); 
                           if(cy_visited.find(key(node_inf))==cy_visited.end())
                              qopen.push(cy_created[key(node_inf)]);
                            
                         } 
	            }
	        }
	    }
	    
	    for(int j=0;j<qtot;j++)
	    {
	        if(j==current_automaton_state && !neg_trans_to_neighbour[j])
	            continue;

	        for(int z=0;z<negtrans[current_automaton_state][j].size();z++)
	        {
	            vector<int> transition_req = negtrans[current_automaton_state][j][z];

	            for(int i=0;i<dir;i++)
	            {
	                if(!valid(cur_pt.y+nb[i][0],cur_pt.x+nb[i][1]))
	                    continue;


	                nbh.y =  cur_pt.y+nb[i][0];
	                nbh.x =  cur_pt.x+nb[i][1];

	                vector<int> neigh_grid_cell{nbh.y,nbh.x};
	                string nbh_s = conv_vec_to_string(neigh_grid_cell);
	                vector<int> nxtstate;
	                vector<int> prop_nbh;
	               
	                if(prop.find(nbh_s)==prop.end())
	                {
	                    prop_nbh.push_back(0);
	                }
	                else
	                {
	                    prop_nbh = prop[nbh_s];
	                    prop_nbh.push_back(0);
	                }

	                int trans_literals_satisfied=0;
	                for(int k=0;k<transition_req.size();k++)
	                {
	                    if(transition_req[k]>=0)
	                        for(int l=0;l<prop_nbh.size();l++)
	                        {
	                            if(transition_req[k]==prop_nbh[l])
	                            {
	                            trans_literals_satisfied++; break;
	                            }
	                        }
	                    else
	                    {
	                        int satisfied_neg_literal=1;
	                        for(int l=0;l<prop_nbh.size();l++)
	                        {
	                            if(-1*transition_req[k]==prop_nbh[l])
	                            {
	                                satisfied_neg_literal=0; break;
	                            }
	                        }
	                        trans_literals_satisfied+=satisfied_neg_literal;
	                    }
	                }
	                if(trans_literals_satisfied==transition_req.size())
	                {
                           vector<int>node_inf={nbh.y,nbh.x,j};
                             if(cy_created.find(key(node_inf))==cy_created.end())
                             { 
                               cycle_node* neigh_nd = new_cycle_node(nbh.x,nbh.y,j);
                               cy_created.insert(make_pair(key(node_inf),neigh_nd));
		               cur_nd->non_adj.push_back(neigh_nd);
                      	       qopen.push(neigh_nd); 
	                     }
                             else
                               {
                                 cur_nd->non_adj.push_back(cy_created[key(node_inf)]);
                           	 if(cy_visited.find(key(node_inf))==cy_visited.end())
                              	     qopen.push(cy_created[key(node_inf)]); 
                               }
	                }
	            }
	        }
	    }
}



void mark_obs(int par_aut_state)
{

     int undesired_prop;
     obs1.clear();
      if(negtrans[par_aut_state][par_aut_state].size() > 0)
      {

	vector<int> transition_req = negtrans[par_aut_state][par_aut_state][0];
	for(int it=0;it<transition_req.size();it++)
	{
	    undesired_prop = transition_req[it]*-1;
	    for(int v=0;v<prop_sys_states[undesired_prop].size();v++)
	    {
	        obs1[key(prop_sys_states[undesired_prop][v])]=-1;

	    }
	}
	for(int ind=1;ind<negtrans[par_aut_state][par_aut_state].size();ind++)
	{
	    unordered_map<long long int,int> obs2;
	    vector<int> transition_req = negtrans[par_aut_state][par_aut_state][ind];
	    for(int it=0;it<transition_req.size();it++)
	    {
	        undesired_prop = transition_req[it]*-1;
	        for(int v=0;v<prop_sys_states[undesired_prop].size();v++)
	        {
	            if(obs1[key(prop_sys_states[undesired_prop][v])])
	            obs2[key(prop_sys_states[undesired_prop][v])]=-1;
	        }
	    }
	    obs1 = obs2;
	}
	for(int z=0;z<qsytrans[par_aut_state][par_aut_state].size();z++)
	{
	    int grid_state = qsytrans[par_aut_state][par_aut_state][z];
	    obs1[key(pos_system_state[grid_state])]=0;

	}

	copy_map(obs1,-1); 
    }


}


void get_path(cycle_node* nd,cycle_node* root)
{
  cycle_node* temp =nd;
  suf_cycle.clear();
  suf_cycle.push_back(make_pair(root->f,root));

  while(temp!=root)
  {  
    suf_cycle.push_back(make_pair(temp->f,temp));
    temp=temp->par;

   }
  suf_cycle.push_back(make_pair(0,root));
  reverse(suf_cycle.begin(),suf_cycle.end());
  float cost = root->f;
  for(auto i:suf_cycle)
  {
    i.second->par=NULL;
    i.second->f=FLT_MAX;
  }
  all_suf_cycle.insert(make_pair(root,make_pair(cost,suf_cycle)));
}



void get_min_cycle(cycle_node* root)
{

bool flag=false;
unordered_map<long long int,bool>visited;
vector<int>node_inf1;
vector<int>node_inf2;


vector<cycle_node*> allnodes;
cycle_queue.clear();

node_inf1.push_back(root->y);
node_inf1.push_back(root->x);
node_inf1.push_back(root->state);
root->f=0;


allnodes.push_back(root);
cycle_queue.insert(make_pair(root,key(node_inf1)));


  while(cycle_queue.size()!=0)
  {
    cycle_node* cur_nd = cycle_queue.begin()->first;
    cycle_queue.erase(cycle_queue.begin());
    node_inf1.clear();
    int par_aut_state=cur_nd->state;
    
    node_inf1.push_back(cur_nd->y);
    node_inf1.push_back(cur_nd->x);
    node_inf1.push_back(cur_nd->state);


    if(cur_nd==root && flag)
    {

      get_path(root->par,root);
      copy_map(obs1,0); 
      for(auto it:allnodes)
      {
         it->par=NULL;
         it->f=FLT_MAX;
      }
      break;
    }

    if(visited.find(key(node_inf1))!=visited.end())
    continue;
    else 
    visited.insert(make_pair(key(node_inf1),true));

    mark_obs(cur_nd->state);
    for(int i=0;i<cur_nd->non_adj.size();i++)
    {  
          node_inf2.clear(); 
          float cost=calc_cost(cur_nd,cur_nd->non_adj[i]);
          if(cost==INF)
          continue;
          node_inf2.push_back(cur_nd->non_adj[i]->y);
          node_inf2.push_back(cur_nd->non_adj[i]->x);
          node_inf2.push_back(cur_nd->non_adj[i]->state);         

           cycle_node* temp=cur_nd->non_adj[i];
           if(temp==root)
           {
            if(root->f >0 && (root->f > cur_nd->f+cost))
             { 
                temp->f = cur_nd->f+cost;
                temp->par=cur_nd; 
                cycle_queue.erase(temp);
                cycle_queue.insert(make_pair(temp,key(node_inf2)));
                allnodes.push_back(temp); 
             }
           else if(root->f==0) 
            {
              temp->f = cur_nd->f+cost;
              temp->par=cur_nd; 
              cycle_queue.insert(make_pair(temp,key(node_inf2)));
              allnodes.push_back(temp); 
            }
        }

        else if(temp->f > cur_nd->f+cost)
        {
          temp->f=cur_nd->f+cost;
          temp->par=cur_nd; 
          cycle_queue.insert(make_pair(temp,key(node_inf2)));
          allnodes.push_back(temp); 
        }
      

        
     } 

    copy_map(obs1,0); 
    if(!flag)flag=true;

  }
  for(auto it:allnodes)
  {
    it->par=NULL;
    it->f=FLT_MAX;
  }
  copy_map(obs1,0); 
}


void display_graph(cycle_node* root)
{

   unordered_map<long long int,bool> visited;
   queue<cycle_node*> myq;
   myq.push(root);
   while(!myq.empty())
   {
     cycle_node* cur=myq.front();
     myq.pop();
     cur->f=FLT_MAX;
     vector<int>node_inf={cur->y,cur->x,cur->state};

     if(visited.find(key(node_inf))!=visited.end())
        continue; 

      else
      	visited.insert(make_pair(key(node_inf),true));
         
         cout<<"\nNode is "<<node_inf[1]<<","<<node_inf[0]<<","<<node_inf[2]<<" and neigh are "<<endl;
         for(auto i:cur->non_adj)
         {
           cout<<i->x<<","<<i->y<<","<<i->state<<"----";
           vector<int>neigh_inf={i->y,i->x,i->state};
           if(visited.find(key(neigh_inf))!=visited.end())
             continue;
          else
                myq.push(i); 
         }

   }

}



void create_graph(cycle_node* node)
{
   vector<int> node_inf={node->y,node->x,node->state};
   long long int node_info=key(node_inf);
   queue<cycle_node*> qopen; 

   node_info=key(node_inf);
   if(cy_created.find(node_info)==cy_created.end())
   cy_created.insert(make_pair(node_info,node));
   qopen.push(node);

   while(!qopen.empty())
   {
     cycle_node* cur=qopen.front();
     qopen.pop();
     vector<int>node_inf={cur->y,cur->x,cur->state};
     node_info=key(node_inf);
     if(cy_visited.find(node_info)!=cy_visited.end())
        continue;

      else
      	cy_visited.insert(make_pair(node_info,true));

        find_non_adj(qopen,cur);

   }

  //display_graph(node);
}


void get_destn_nodes()
{
    int i;
    pair<int,int> dest_node;
    vector<cycle_node*> myqueue;
    map<long long int,bool>vis_node;

    /*******for every final state of automaton***************/
    for(int fs=0;fs<dest.size();fs++)
    {
        vector<int> vis_sy_state(pos_system_state.size());
        // for every node through which dest[fs] could be reached 
        for(i=0;i<qsystate[dest[fs]].size();i++)
        {
            //cout<<"--systate--"<<qsystate[dest[fs]].size()<<endl; 
            int systate = qsystate[dest[fs]][i];
            if(vis_sy_state[systate]==1)
                continue;
            else
                vis_sy_state[systate]=1;
	
            dstar_inf st;
            st.x=pos_system_state[systate][0];
            st.y=pos_system_state[systate][1];
            st.state=dest[fs];
            dest_coords.push_back(st);
           
            cycle_node *temp=NULL;
            graph_root[st]=temp;

        }
    }


}





void mark_destn(cycle_node* root)
{

 queue<cycle_node*> tr_queue;
 tr_queue.push(root);
 map<cycle_node*,bool> visited;
  while(!tr_queue.empty())
 {

    cycle_node* cur=tr_queue.front();
    tr_queue.pop();

    if(visited.find(cur)!=visited.end())
        continue;
    else
        visited[cur]=true;

    dstar_inf temp;
    temp.x=cur->x;
    temp.y=cur->y;
    temp.state=cur->state;
    if(find(dest_coords.begin(),dest_coords.end(),temp)!=dest_coords.end())
    {
        graph_root[temp]=cur;
    }

    for(auto i:cur->non_adj)
        tr_queue.push(i);
 }

}



float calc_dy_cost(cycle_node* source, cycle_node* dest,ull from,ull till=plan_till)
{
   if(source->x==dest->x && source->y==dest->y)
   return(0);
   
   if(grid_dy.find(make_pair(source->x,source->y))!=grid_dy.end())
   {

     if(grid_dy[make_pair(source->x,source->y)]>from)
     {
          return(FLT_MAX);     
     }

   }
    unordered_map<long long int,node* > init_vertex;
    unordered_map<long long int,int > vis;
    
    long long wait_time=0;
    vector<int> node_inf={source->y,source->x,source->state};
    ull cur_ts=from; 


    node* start_nd = new node;
    start_nd->x=source->x;
    start_nd->y=source->y;
    start_nd->state=source->state;
    start_nd->g=0;  
    
    node* goal_nd= new node;
    goal_nd->x=dest->x;
    goal_nd->y=dest->y;
    goal_nd->state=dest->state;
    goal_nd->h=0;
    goal_nd->g=FLT_MAX;

    
    start_nd->h=get_initial_heuristic(start_nd,goal_nd);
    start_nd->f=start_nd->g+start_nd->h;


    node_oq qopen;
    vector<node*> closed;
    qopen.insert(make_pair(start_nd,key(node_inf)));
    vector<int> src_vertex={start_nd->y,start_nd->x,start_nd->state};  
    vector<int> dest_vertex={goal_nd->y,goal_nd->x,goal_nd->state};        
    
    init_vertex[key(src_vertex)]=start_nd;
    init_vertex[key(dest_vertex)]=goal_nd;
    
    while(!qopen.empty() && cur_ts<=till)
    {
        node* cur_nd = qopen.begin()->first;
        qopen.erase(qopen.begin());

        cur_ts=from+cur_nd->g;
        int current_automaton_state = cur_nd->state;
        vector<int> cur_node_info{cur_nd->y,cur_nd->x,current_automaton_state};
        
        if(cur_nd->x==dest->x && cur_nd->y==dest->y)
        {
            node* temp=cur_nd->par;
            while(temp!=start_nd)
            {
             temp=temp->par;
            }
            for(auto it:closed)
            delete it;
            return(cur_nd->f);
        }

        long long int cur_node_k = key(cur_node_info);
        
        if(vis.find(cur_node_k)!=vis.end())
            continue;
        else
            vis[cur_node_k]=1;

        vector<int> neigh_state;

        point nbh;
        for(int i=0;i<dir;i++)
        {
            wait_time=0;
            if(!valid(cur_nd->y+nb[i][0],cur_nd->x+nb[i][1]) && ((cur_nd->y+nb[i][0])!=dest->y || (cur_nd->x+nb[i][1]!=dest->x)))
                continue;

            nbh.y =  cur_nd->y+nb[i][0];
            nbh.x =  cur_nd->x+nb[i][1];
            neigh_state = vector<int>(0);
            neigh_state.push_back(nbh.y);
            neigh_state.push_back(nbh.x);
            neigh_state.push_back(current_automaton_state);

            node* oldtmp;
            if(grid_dy.find(make_pair(nbh.x,nbh.y))!=grid_dy.end())
            {
             
             wait_time=(grid_dy[make_pair(nbh.x,nbh.y)]-cur_ts);
             if(wait_time<=0)
             wait_time=0;
            }
              
            if(init_vertex.find(key(neigh_state))==init_vertex.end())
            {   
                node* neigh_node = new_node(nbh.x,nbh.y,current_automaton_state);           
                neigh_node->g = cur_nd->g+1+wait_time;

                neigh_node->h = get_initial_heuristic(neigh_node,goal_nd);           
                neigh_node->f = neigh_node->g+neigh_node->h;
                neigh_node->par = cur_nd;
                init_vertex.insert(make_pair(key(neigh_state),neigh_node) );
                qopen.insert(make_pair(neigh_node,key(neigh_state)));

            }
            else
            {

                oldtmp = init_vertex[key(neigh_state)];
                if(oldtmp->g > (cur_nd->g+1+wait_time) )
                {
                    if(qopen.find(oldtmp)!=qopen.end())
                    qopen.erase(qopen.find(oldtmp));

                    oldtmp->g = cur_nd->g+1+wait_time;
                    oldtmp->f = oldtmp->g+oldtmp->h;
                    oldtmp->par = cur_nd;
                    init_vertex[key(neigh_state)]=oldtmp;
                    qopen.insert(make_pair(oldtmp,key(neigh_state)));
                    
                }

            }

        }

    closed.push_back(cur_nd);
    
   }
return(FLT_MAX);


}




void get_static_edge_cost(cycle_node* root)
{
   unordered_map<long long int,bool> visited;
   queue<cycle_node*> myq;
   myq.push(root);
   while(!myq.empty())
   {
     cycle_node* cur=myq.front();
     myq.pop();
     cur->f=FLT_MAX;
     vector<int>cur_inf={cur->y,cur->x,cur->state};

     if(visited.find(key(cur_inf))!=visited.end())
        continue; 

      else
      	visited.insert(make_pair(key(cur_inf),true));
         
         int par_aut_state=cur->state;
    	 mark_obs(par_aut_state);  
    	    
         for(auto i:cur->non_adj)
         {
               vector<int>neigh_inf={i->y,i->x,i->state};
               float cost=calc_cost(cur,i);
               static_edge_len.insert(make_pair(make_pair(key(cur_inf),key(neigh_inf)),cost));
               myq.push(i);
         }
     
       copy_map(obs1,0); 
   }
}



float get_dy_cy(cycle_node* root,ull from,ull till)
{
bool flag=false;
unordered_map<long long int,bool>visited;
vector<int>node_inf1;
vector<int>node_inf2;
unsigned long long cur_ts=from;

vector<cycle_node*> allnodes;
cycle_queue.clear();

node_inf1.push_back(root->y);
node_inf1.push_back(root->x);
node_inf1.push_back(root->state);
root->f=0;


allnodes.push_back(root);
cycle_queue.insert(make_pair(root,key(node_inf1)));


  while(cycle_queue.size()!=0 && (cur_ts<=till))
  {
    cycle_node* cur_nd = cycle_queue.begin()->first;
    cur_ts=cur_nd->f + from;
    cycle_queue.erase(cycle_queue.begin());
    
    node_inf1.clear();
    node_inf1.push_back(cur_nd->y);
    node_inf1.push_back(cur_nd->x);
    node_inf1.push_back(cur_nd->state);


    if(cur_nd==root && flag)
    {
       suf_dy_path.clear();
       suf_dy_path.push_back(make_pair(cur_nd->f,root));
       cycle_node* temp =root->par;
       copy_map(obs1,0);        
       while(temp!=root)
       {
         suf_dy_path.push_back(make_pair(temp->f,temp));
         temp=temp->par;
       }
       suf_dy_path.push_back(make_pair(0,root));
       reverse(suf_dy_path.begin(),suf_dy_path.end());
       float cost=root->f;
       for(auto it:allnodes)
       {
         it->par=NULL;
         it->f=FLT_MAX;
       }
       return(cost);
    }

    if(visited.find(key(node_inf1))!=visited.end())
    continue;
    
    mark_obs(cur_nd->state);
    for(int i=0;i<cur_nd->non_adj.size();i++)
    {

           node_inf2.clear();       
	  float cost=calc_dy_cost(cur_nd,cur_nd->non_adj[i],cur_ts);
          if(cost==FLT_MAX)continue;

              node_inf2.push_back(cur_nd->non_adj[i]->y);
	      node_inf2.push_back(cur_nd->non_adj[i]->x);
	      node_inf2.push_back(cur_nd->non_adj[i]->state);         

           cycle_node* temp=cur_nd->non_adj[i];
           if(temp==root)
           {
            if(root->f >0 && (root->f > cur_nd->f+cost))
             { 
                temp->f = cur_nd->f+cost;
                temp->par=cur_nd; 
                cycle_queue.erase(temp);
                cycle_queue.insert(make_pair(temp,key(node_inf2)));
                allnodes.push_back(temp); 
             }
           else if(root->f==0) 
            {
              temp->f = cur_nd->f+cost;
              temp->par=cur_nd; 
              cycle_queue.insert(make_pair(temp,key(node_inf2)));
              allnodes.push_back(temp); 
            }

        }

        else if(temp->f > cur_nd->f+cost)
        {
          temp->f=cur_nd->f+cost;
          temp->par=cur_nd; 
          cycle_queue.insert(make_pair(temp,key(node_inf2)));
          allnodes.push_back(temp); 
        }
        
     } 
    copy_map(obs1,0); 
    if(!flag)flag=true;
    visited.insert(make_pair(key(node_inf1),true));

  }

 return(FLT_MAX);

}





float get_dy_pref(cycle_node* start, cycle_node* dest,ull from )
{

 pwh_pref.clear();
if(start->x==dest->x && start->y==dest->y)
return(0);

unordered_map<long long int,bool>visited;
vector<int>node_inf1;
vector<int>node_inf2;
unsigned long long cur_ts=from;

vector<cycle_node*> allnodes;
cycle_queue.clear();

node_inf1.push_back(start->y);
node_inf1.push_back(start->x);
node_inf1.push_back(start->state);

allnodes.push_back(start);
start->f=0;
cycle_queue.insert(make_pair(start,key(node_inf1)));


  while(cycle_queue.size()!=0 && cur_ts<=plan_till)
  {
    cycle_node* cur_nd = cycle_queue.begin()->first;
    cycle_queue.erase(cycle_queue.begin());
    cur_ts=from+cur_nd->f;

    if(cur_ts>plan_till)
    break;
    
    
    node_inf1.clear();
    node_inf1.push_back(cur_nd->y);
    node_inf1.push_back(cur_nd->x);
    node_inf1.push_back(cur_nd->state);


    if(cur_nd==dest)
    {
       pwh_pref.clear();
       pwh_pref.insert(make_pair(cur_nd->f,cur_nd));
       cycle_node* temp =cur_nd->par;
       copy_map(obs1,0); 
       while(temp!=start)
       {
         pwh_pref.insert(make_pair(temp->f,temp));
         temp=temp->par;
       }
       pwh_pref.insert(make_pair(0,start));

       float cost=cur_nd->f;
       for(auto it:allnodes)
       {
         it->par=NULL;
         it->f=FLT_MAX;
       }
       return(cost);
    }

    if(visited.find(key(node_inf1))!=visited.end())
    continue;

    mark_obs(cur_nd->state);
    for(int i=0;i<cur_nd->non_adj.size();i++)
    {

          node_inf2.clear();       
	  float cost=calc_dy_cost(cur_nd,cur_nd->non_adj[i],cur_ts);
          if(cost==FLT_MAX)continue;

	  node_inf2.push_back(cur_nd->non_adj[i]->y);
	  node_inf2.push_back(cur_nd->non_adj[i]->x);
	  node_inf2.push_back(cur_nd->non_adj[i]->state);         

	  cycle_node* temp=cur_nd->non_adj[i];

	  if(temp->f > cur_nd->f+cost)
	  {
	    temp->f=cur_nd->f+cost;
	    temp->par=cur_nd; 
	    if(cycle_queue.find(temp)!=cycle_queue.end())
	    cycle_queue.erase(temp);
	    
	    cycle_queue.insert(make_pair(temp,key(node_inf2)));
	    allnodes.push_back(temp); 
         }

        
     }
 
    copy_map(obs1,0); 
    visited.insert(make_pair(key(node_inf1),true));

  }
 copy_map(obs1,0); 
 for(auto it:allnodes)
 {
     it->par=NULL;
     it->f=FLT_MAX;
  }
 return(FLT_MAX);

}




float get_dy_suff(cycle_node* root,ull from,ull prev_cost=0,ull till=plan_till)
{

bool flag=false;
unordered_map<long long int,bool>visited;
vector<int>node_inf1;
vector<int>node_inf2;
unsigned long long cur_ts=from;

vector<cycle_node*> allnodes;
cycle_queue.clear();

node_inf1.push_back(root->y);
node_inf1.push_back(root->x);
node_inf1.push_back(root->state);
root->f=0;


allnodes.push_back(root);
cycle_queue.insert(make_pair(root,key(node_inf1)));


  while(cycle_queue.size()!=0 )
  {
    cycle_node* cur_nd = cycle_queue.begin()->first;
    cur_ts=cur_nd->f + from;
    if(cur_ts>plan_till)
    break;
    
    cycle_queue.erase(cycle_queue.begin());
    
    node_inf1.clear();
    if(cur_nd->par!=NULL)

    node_inf1.push_back(cur_nd->y);
    node_inf1.push_back(cur_nd->x);
    node_inf1.push_back(cur_nd->state);


    if(cur_nd==root && flag)
    {
       pwh_suf.clear();
       pwh_suf.insert(make_pair(root->f+prev_cost,root));
       cycle_node* temp =root->par;
       copy_map(obs1,0);
       while(temp!=root)
       {
         pwh_suf.insert(make_pair(temp->f+prev_cost,temp));
         temp=temp->par;
       }
       pwh_suf.insert(make_pair(prev_cost,root));
       float cost=root->f;
       for(auto it:allnodes)
       {
         it->par=NULL;
         it->f=FLT_MAX;
       }
       return(cost);
    }

    if(visited.find(key(node_inf1))!=visited.end())
    continue;
    
    mark_obs(cur_nd->state);
    for(int i=0;i<cur_nd->non_adj.size();i++)
    {

          node_inf2.clear();       
	  float cost=calc_dy_cost(cur_nd,cur_nd->non_adj[i],cur_ts);
          if(cost==FLT_MAX)
          continue;

          node_inf2.push_back(cur_nd->non_adj[i]->y);
	  node_inf2.push_back(cur_nd->non_adj[i]->x);
	  node_inf2.push_back(cur_nd->non_adj[i]->state);         

           cycle_node* temp=cur_nd->non_adj[i];
           if(temp==root)
           {
            if(root->f >0 && (root->f > cur_nd->f+cost))
             { 
                temp->f = cur_nd->f+cost;
                temp->par=cur_nd; 
                
                if(cycle_queue.find(temp)!=cycle_queue.end())
                cycle_queue.erase(temp);
                
                cycle_queue.insert(make_pair(temp,key(node_inf2)));
                allnodes.push_back(temp); 
             }
           else if(root->f==0) 
            {
              temp->f = cur_nd->f+cost;
              temp->par=cur_nd; 
              cycle_queue.insert(make_pair(temp,key(node_inf2)));
              allnodes.push_back(temp); 
            }
        }

        else if(temp->f > cur_nd->f+cost)
        {
          temp->f=cur_nd->f+cost;
          temp->par=cur_nd; 
          
          if(cycle_queue.find(temp)!=cycle_queue.end())
          cycle_queue.erase(temp);          
          
          cycle_queue.insert(make_pair(temp,key(node_inf2)));
          allnodes.push_back(temp); 
        }
        
     } 
    copy_map(obs1,0);
    if(!flag)flag=true;
    visited.insert(make_pair(key(node_inf1),true));

  }
  
  for(auto it:allnodes)
 {
  it->par=NULL;
  it->f=FLT_MAX;
 }

 return(FLT_MAX);

}



bool is_dest(cycle_node* nd)
{

    for(auto i:graph_root)
       if(nd==i.second)
       return(true);
       
 return(false);

}



float get_pref_cost(cycle_node* src,cycle_node* dest)
{

unordered_map<long long int,bool>visited;
vector<int>node_inf1;
vector<int>node_inf2;
vector<cycle_node*> allnodes;
cycle_queue.clear();

node_inf1.push_back(src->y);
node_inf1.push_back(src->x);
node_inf1.push_back(src->state);
src->f=0;

allnodes.push_back(src);
cycle_queue.insert(make_pair(src,key(node_inf1)));


  while(cycle_queue.size()!=0)
  {
    cycle_node* cur_nd = cycle_queue.begin()->first;
    cycle_queue.erase(cycle_queue.begin());
    
    node_inf1.clear();
    node_inf1.push_back(cur_nd->y);
    node_inf1.push_back(cur_nd->x);
    node_inf1.push_back(cur_nd->state);

    if(cur_nd==dest)
    {
       pref_path.clear();
       cycle_node* temp =dest;
       float cost=cur_nd->f;
       copy_map(obs1,0);
       while(temp!=src)
       {
         pref_path.push_back(make_pair(temp->f,temp));
         temp=temp->par;
       }
       pref_path.push_back(make_pair(0,src));
       reverse(pref_path.begin(),pref_path.end());
       for(auto it:allnodes)
       {
         it->par=NULL;
         it->f=FLT_MAX;
       }
       return(cost);
    }

    if(visited.find(key(node_inf1))!=visited.end())
    continue;
    
    mark_obs(cur_nd->state);
    for(int i=0;i<cur_nd->non_adj.size();i++)
    {

          node_inf2.clear();       
          node_inf2.push_back(cur_nd->non_adj[i]->y);
          node_inf2.push_back(cur_nd->non_adj[i]->x);
          node_inf2.push_back(cur_nd->non_adj[i]->state); 

	  float cost=static_edge_len[make_pair(key(node_inf1),key(node_inf2))];
          if(cost==FLT_MAX)continue;      

           cycle_node* temp=cur_nd->non_adj[i];
           if(temp->f > cur_nd->f+cost)
           {
              if(cycle_queue.find(temp)!=cycle_queue.end())
               cycle_queue.erase(temp);
             
              temp->f=cur_nd->f+cost;
              temp->par=cur_nd; 
              cycle_queue.insert(make_pair(temp,key(node_inf2)));
              allnodes.push_back(temp); 
          }
        
     } 
    visited.insert(make_pair(key(node_inf1),true));
    copy_map(obs1,0);
  }
  for(auto it:allnodes)
  {
   it->par=NULL;
   it->f=FLT_MAX;
  }
 return(FLT_MAX);

}



dstar_inf get_position(cycle_node* source, cycle_node* dest,ull from, int index)
{

    unordered_map<long long int,node* > init_vertex;
    unordered_map<long long int,int > vis;
    vector<node*>temp_path;
    dstar_inf pos;
    int i;
    pos.x=source->x;
    pos.y=source->y;
    pos.state=source->state;
    
    
    long long wait_time=0;
    vector<int> node_inf={source->y,source->x,source->state};
    ull cur_ts=from; 

    node* start_nd = new node;
    start_nd->x=source->x;
    start_nd->y=source->y;
    start_nd->state=source->state;
    start_nd->g=0; 
    
    mark_obs(source->state); 
    
    node* goal_nd= new node;
    goal_nd->x=dest->x;
    goal_nd->y=dest->y;
    goal_nd->state=dest->state;
    goal_nd->h=0;
    goal_nd->g=FLT_MAX;

    
    start_nd->h=get_initial_heuristic(start_nd,goal_nd);
    start_nd->f=start_nd->g+start_nd->h;


    node_oq qopen;
    vector<node*> closed;
    qopen.insert(make_pair(start_nd,key(node_inf)));

    vector<int> src_vertex={start_nd->y,start_nd->x,start_nd->state};  
    vector<int> dest_vertex={goal_nd->y,goal_nd->x,goal_nd->state};        
    
    init_vertex[key(src_vertex)]=start_nd;
    init_vertex[key(dest_vertex)]=goal_nd;
    
    while(!qopen.empty())
    {
        node* cur_nd = qopen.begin()->first;
        qopen.erase(qopen.begin());

        cur_ts=from+cur_nd->g;
        int current_automaton_state = cur_nd->state;
        vector<int> cur_node_info{cur_nd->y,cur_nd->x,current_automaton_state};
        
        if(cur_nd->x==dest->x && cur_nd->y==dest->y)
        {
    
            copy_map(obs1,0); 
            node* temp=cur_nd->par;

            temp_path.push_back(cur_nd);
            while(temp!=start_nd)
            {
             temp_path.push_back(temp);              
             temp=temp->par;
            }
            temp_path.push_back(start_nd);
            reverse(temp_path.begin(),temp_path.end());
            
            for(i=0;i<temp_path.size()-1;i++)
            {
             if(temp_path[i]->g == index)
             break;
             else if(temp_path[i]->g<index && temp_path[i+1]->g>index)
             break;
            }
            
            if(i==temp_path.size()-1)
            {
             if(index!=temp_path[i]->g)
             i--;
            
            }
            
            pos.x=temp_path[i]->x;
            pos.y=temp_path[i]->y;
            pos.state=temp_path[i]->state;
            
            for(auto it:closed)
            delete it;
            return(pos);
        }

        long long int cur_node_k = key(cur_node_info);
        
        if(vis.find(cur_node_k)!=vis.end())
            continue;
        else
            vis[cur_node_k]=1;

        vector<int> neigh_state;

        point nbh;
        for(int i=0;i<dir;i++)
        {
            wait_time=0;
            if(!valid(cur_nd->y+nb[i][0],cur_nd->x+nb[i][1]) && ((cur_nd->y+nb[i][0])!=dest->y || (cur_nd->x+nb[i][1]!=dest->x)))
                continue;
                
            nbh.y =  cur_nd->y+nb[i][0];
            nbh.x =  cur_nd->x+nb[i][1];
            neigh_state = vector<int>(0);
            neigh_state.push_back(nbh.y);
            neigh_state.push_back(nbh.x);
            neigh_state.push_back(current_automaton_state);

            node* oldtmp;
            if(grid_dy.find(make_pair(nbh.x,nbh.y))!=grid_dy.end())
            {
             
             wait_time=(grid_dy[make_pair(nbh.x,nbh.y)]-cur_ts);
             if(wait_time<=0)
             wait_time=0;
            }
              
            if(init_vertex.find(key(neigh_state))==init_vertex.end())
            {   
                node* neigh_node = new_node(nbh.x,nbh.y,current_automaton_state);           
                neigh_node->g = cur_nd->g+1+wait_time;

                neigh_node->h = get_initial_heuristic(neigh_node,goal_nd);           
                neigh_node->f = neigh_node->g+neigh_node->h;
                neigh_node->par = cur_nd;
                init_vertex.insert(make_pair(key(neigh_state),neigh_node) );
                qopen.insert(make_pair(neigh_node,key(neigh_state)));

            }
            else
            {

                oldtmp = init_vertex[key(neigh_state)];
                if(oldtmp->g > (cur_nd->g+1+wait_time) )
                {
                    if(qopen.find(oldtmp)!=qopen.end())
                    qopen.erase(qopen.find(oldtmp));

                    oldtmp->g = cur_nd->g+1+wait_time;
                    oldtmp->f = oldtmp->g+oldtmp->h;
                    oldtmp->par = cur_nd;
                    init_vertex[key(neigh_state)]=oldtmp;
                    qopen.insert(make_pair(oldtmp,key(neigh_state)));
                    
                }

            }

        }

    closed.push_back(cur_nd);
    
   }
       copy_map(obs1,0); 
       return(pos);

}





int main(int args,char **argv)
{
    long long obs_ts,check_time,dec_time;
    obs_ts=LONG_MAX;
    int x,y,num_obs;
    ull dur;
    long long plan_comp=0;
    float tcost=FLT_MAX;    
    vector<pair<int,cycle_node*>> static_pref;
    bool plan_found=false;
    float cy_cost=FLT_MAX;
    float pref_cost=FLT_MAX;
    
    fstream res;
    
    res.open("result.txt",ios::out|ios::app);    
    if(!res)
    {
 	cout<<"Open failed for result file\n";
        return(0);
    }
    
    if(args<3)
    {
    cout<<"cmd args are not properly passed through command line\n";
    return(0);
    
    }
    
    bool recheck=false;
    initializegrid(argv[1]);

    calsystates();
    create_graph(origin); // Created a graph with origin as root

    get_destn_nodes(); 
    mark_destn(origin); //created destn nodes here and stored in graph_root
   
    plan_till=stoi(argv[3]);
    for(auto i:graph_root)
       get_min_cycle(i.second);
       
    get_static_edge_cost(origin);// Got all the static edge cost
    static_pref.clear();
    for(auto i:cy_created)
    i.second->f=FLT_MAX;
    for(auto i:all_suf_cycle)
    {
       float cost=get_pref_cost(origin,i.first);
       if(cost==FLT_MAX)
       continue;
       
       if((cost+i.second.first)<tcost)
       {
               pref_cost=cost;       
               current_best=i.first;
               tcost=i.second.first+cost;
               static_pref.clear();
       	       for(auto j:pref_path)
        	static_pref.push_back(j);
         
               
       }
       
    }
   
   if(static_pref.size()==0)
   {
    res<<"No Destn is reachable .. Returning \n";
    return(0);  
   }   
   else
   {
     pwh_plan.clear();
     for(auto i:static_pref)
     pwh_plan[i.first]=i.second;

   }

    reader2.open(argv[2]);  
   if(!reader2.eof())
    reader2>>obs_ts;
    if(obs_ts==-1)
    obs_ts=plan_till+1;
    
     check_time=LONG_MAX; 
     dec_time=(dynamic_obstacles.size()==0?LONG_MAX:dynamic_obstacles.begin()->first);
  
     check_time=min({obs_ts,plan_till,dec_time});
     
     while(1)
     { 
       pref_cost+=all_suf_cycle[current_best].first;
       pwh_plan.insert(make_pair(pref_cost,current_best));
       if(pref_cost>check_time)
       break;

      }
       
	while(global_ts<plan_till)
	{
	  cycle_node* prev=pwh_plan.begin()->second;
	  cycle_node* cur;
	  int ts;
	  recheck=false;
  	  pwh_plan.erase(pwh_plan.begin());
	  while(global_ts<check_time && pwh_plan.size()!=0)
	  {

	    cur=pwh_plan.begin()->second;
	    ts=pwh_plan.begin()->first;
	    if(ts>=plan_till)
            { 
              if(ts==plan_till)
              {
	        res<<"Moving from "<<prev->x<<","<<prev->y<<","<<prev->state<<" to "<<cur->x<<","<<cur->y<<","<<cur->state<<" at ts = "<<ts<<endl;
	        if(is_dest(cur) && prev!=NULL && prev_f==cur )
                {
                 cycle_count++;
                }
              }
               res<<" Greedy2 Number of cycles completed are "<<cycle_count<<" till "<<global_ts<<endl;
               return(0);
            }
	    if(ts>=check_time)   
	     break;
	    else
	    {
	      res<<"Moving from "<<prev->x<<","<<prev->y<<","<<prev->state<<" to "<<cur->x<<","<<cur->y<<","<<cur->state<<" at ts = "<<ts<<endl;
	      if(is_dest(cur))
	      {
	       if(prev_f==NULL)
	       prev_f=cur;
	       else if(prev_f==cur)
               {
                 cycle_count++;
               }
	       else
	       prev_f=cur;
	      }
	      prev=cur;
	      global_ts=ts;
	      pwh_plan.erase(pwh_plan.begin());
	    } 
	  }
	  if(pwh_plan.empty())
	  {
	   recheck=true;
	   new_origin=cur;
	   cur->f=FLT_MAX;
	  }						
	  else if(ts==check_time)
	  {
	    new_origin=cur;
	    new_origin->f=FLT_MAX;
	    global_ts=check_time;
	    res<<"Moving from "<<prev->x<<","<<prev->y<<","<<prev->state<<" to "<<new_origin->x<<","<<new_origin->y<<","<<new_origin->state<<" at ts = "<<global_ts<<endl;
	    if(is_dest(new_origin))
            { 
              if(prev_f!=NULL && prev_f==new_origin)
	        cycle_count++;
                prev_f=new_origin;
             }
	  }
	  else 
	  {
	    if(prev->x==cur->x && prev->y==cur->y)
	    { 
	      float cy_cost=get_dy_cy(prev,global_ts,plan_till);

	      for(int i=1;i<suf_dy_path.size();i++)
	      {
	       if(global_ts+suf_dy_path[i].first==check_time)
		{
		  global_ts=check_time;
		  new_origin=suf_dy_path[i].second;
		  new_origin->f=FLT_MAX;
		  break;
		}
		else if(global_ts+suf_dy_path[i].first>check_time)
		{
		  global_ts=global_ts+suf_dy_path[i-1].first;
		  int index=check_time-global_ts;

		  dstar_inf pos=get_position(suf_dy_path[i-1].second,suf_dy_path[i].second,global_ts,index);
		  vector<int> node_inf={pos.y,pos.x,pos.state};
		  if(cy_created.find(key(node_inf))!=cy_created.end())
		  {
		    new_origin=cy_created[key(node_inf)];
		    new_origin->f=FLT_MAX;
		    
		  }
		  else
		  {
		    new_origin=new_cycle_node(pos.x,pos.y,pos.state);
		    cy_created.insert(make_pair(key(node_inf),new_origin));
		  }
		  global_ts=check_time;
		  break;
		}
	      }
	     
	    }
	    else
	    {
	    
		  int index=check_time-global_ts;

		  dstar_inf pos=get_position(prev,cur,global_ts,index);
		  vector<int> node_inf={pos.y,pos.x,pos.state};
		  if(cy_created.find(key(node_inf))!=cy_created.end())
		  { 
		    new_origin=cy_created[key(node_inf)];
		    new_origin->f=FLT_MAX;
		  }  
		  
		  else
		  {

		    new_origin=new_cycle_node(pos.x,pos.y,pos.state);
		    cy_created.insert(make_pair(key(node_inf),new_origin));
		  }
		  global_ts=check_time;
	    }
	    
             res<<"Moving from "<<prev->x<<","<<prev->y<<","<<prev->state<<" to "<<new_origin->x<<","<<new_origin->y<<","<<new_origin->state<<" at ts = "<<check_time<<endl;
             if(is_dest(new_origin))
             { 
              if(prev_f!=NULL && prev_f==new_origin)
               {
                 cycle_count++;
               }
                 prev_f=new_origin;
             }
          
	  }  
  	  if(check_time==global_ts && check_time==dec_time)
          {
            recheck=true;
	    for(auto it:dynamic_obstacles.begin()->second)
	    {
	      if(grid_dy.find(make_pair(it.first,it.second))!=grid_dy.end() && grid_dy[{it.first,it.second}]<=global_ts)
	      grid_dy.erase(grid_dy.find(make_pair(it.first,it.second)));
	     }	    
             dynamic_obstacles.erase(dynamic_obstacles.begin()); 
             dec_time=(dynamic_obstacles.size()==0?LONG_MAX:dynamic_obstacles.begin()->first);

	  }        
	  if(check_time==global_ts && check_time==obs_ts)
	  {
	      recheck=true;
	      reader2>>num_obs;
	      ull max_time;	      

	      for(int i=0;i<num_obs;i++)
	      {
		reader2>>x;
		reader2>>y;
		reader2>>dur;
		
		vector<int> nd_inf={x,y};
		string s=conv_vec_to_string(nd_inf);
		
	        if(x==new_origin->x && y==new_origin->y)
		{
		 res<<"Greedy2 cant mark current robot position as obstacle "<<x<<","<<y<<" Map is invalid "<<endl;
		 return(0);
		}
		if(grid_dy.find(make_pair(x,y))!=grid_dy.end())
                {
                   if(grid_dy[make_pair(x,y)]<(long long)(obs_ts+dur))
                   {
                      grid_dy[make_pair(x,y)]=obs_ts+dur;
                      if(dynamic_obstacles.find(obs_ts+dur)!=dynamic_obstacles.end())
		      dynamic_obstacles[obs_ts+dur].push_back(make_pair(x,y));
		      else
		      {
		        vector<pair<int,int>> temp;
		        temp.push_back(make_pair(x,y));
		        dynamic_obstacles.insert(make_pair(obs_ts+dur,temp));
		       }                   
                   } 
                
                }
                else
                {
                   grid_dy[make_pair(x,y)]=obs_ts+dur;
                   dynamic_obstacles[(obs_ts+dur)].push_back(make_pair(x,y));
                }
		

	      }

	       reader2>>obs_ts;
	       if(obs_ts==-1)
	       obs_ts=plan_till+1;
	       dec_time=(dynamic_obstacles.size()==0?LONG_MAX:dynamic_obstacles.begin()->first);

	    }           
	    if(recheck)
	    { 
	           if(global_ts==plan_till)
	           {
	           res<<"Greedy2 Number of cycle traversed = "<<cycle_count<<endl;
	           return(0);
	           }
	           
	           if(pwh_plan.size()!=0)
	           create_graph(new_origin);
	           
		   float best_cy_cost=FLT_MAX;
		   pwh_plan.clear();
		   pwh_pref.clear();
		   pwh_suf.clear();
		   cycle_node* cur_best;
		   float pref_cost;
		   float plan_cost=FLT_MAX;

		   for(auto i:graph_root)
		   {

		       pref_cost=get_dy_pref(new_origin,i.second,global_ts);
		       
		       if(pref_cost==FLT_MAX)
		       continue;

		       float suff_cost=get_dy_suff(i.second,global_ts+pref_cost,pref_cost);
		       if(suff_cost==FLT_MAX)
		       continue;
		       if((pref_cost+suff_cost)<plan_cost)
		       {   
		             pwh_plan.clear();
		             if(pref_cost==0)
		             pwh_plan.insert(make_pair(global_ts,new_origin));
		             
		             else
		             for(auto i:pwh_pref)
		             pwh_plan.insert(make_pair(global_ts+i.first,i.second));             
		             plan_cost=pref_cost+suff_cost;
		             cur_best=i.second;
		           
		       }
		       
		     
		   }
		   if(plan_cost==FLT_MAX)
		   {

		      while(obs_ts<=plan_till)
		      {
			      reader2>>num_obs;
			      ull max_time;	      
			      for(int i=0;i<num_obs;i++)
			      {
				reader2>>x;
				reader2>>y;
				reader2>>dur;

				if(x==new_origin->x && y==new_origin->y)
				{
				 res<<"Greedy2 cant mark current robot position as obstacle "<<x<<","<<y<<" Map is invalid "<<endl;
				 return(0);
				}

			       }
			       reader2>>obs_ts;
			       if(obs_ts==-1)
			       obs_ts=plan_till+1;
		      } 
		      res<<"Greedy2 Number of cycle traversed = "<<cycle_count<<endl;
		      return(0); 
		   }
                   else
                    pwh_plan.insert(make_pair(global_ts+plan_cost,cur_best));
                    
                    dec_time=(dynamic_obstacles.size()==0?LONG_MAX:dynamic_obstacles.begin()->first);
                    check_time=min({plan_till,obs_ts,dec_time}); 
	     } 

	  
	}
	res<<"Greedy2 Number of cycle traversed = "<<cycle_count<<endl;

return(0);
}


