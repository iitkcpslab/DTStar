#include"z3_opt.h"
#include<unistd.h>
#define mx_steps 50000
#define INF 100000 
#define ull unsigned long long
#include<ctime>

using namespace std;
using namespace planner_info;

fstream reader2,fptr2;
cycle_node* origin=new_cycle_node(0,0,0);
cycle_node* new_origin;
int md_count=0;
string model_file;
int tcomp=1;
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

    if ((fp = popen(ltl_query, "r")) == NULL)
    {
        printf("Error opening pipe!\n");
        return automata_info;
    }
    
    while (fgets(buf, BUFSIZE, fp) != NULL) 
    {
        table_text = table_text+buf;
    }

    if(pclose(fp))  
    {
        printf("Command not found or exited with error status\n");
        return automata_info;
    }
    unordered_map<string,int> mp;
    unordered_map<string,int> state;

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
    if (ft == NULL) 
    {
      perror("Failed: ");
      return;
    }
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
    
    /**reading automata transitions**/
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
        /**storing transition condition  in automata**/
        trans[automata_states[0]][automata_states[1]].push_back(transition_condn);

      

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
        
        prop[s].push_back(literal);

        vector<int> ivec(2);
        ivec[0]=grid_state[1];
        ivec[1]=grid_state[0];
        prop_sys_states[literal].push_back(ivec);

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

    qsystate = vector< vector<int> >(qtot);
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
   
   vector<int> src_vertex={source->y,source->x,source->state};  
   vector<int> dest_vertex={dest->y,dest->x,dest->state}; 
   vector<vector<int>> path;

   if(source->x==dest->x && source->y==dest->y)
    {
      static_edge_path.insert(make_pair(make_pair(key(src_vertex),key(dest_vertex)),path));    	
      return(0);
    }
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
    
    init_vertex[key(src_vertex)]=start_nd;
    init_vertex[key(dest_vertex)]=goal_nd;
    
    while(!qopen.empty())
    {
        node* cur_nd = qopen.begin()->first;
        qopen.erase(qopen.begin());
        int current_automaton_state = cur_nd->state;
        //map<pair<unsigned long long,unsigned long long>,vector<pair<int,int>>>static_edge_path;
        vector<int> cur_node_info{cur_nd->y,cur_nd->x,current_automaton_state};
        
        if(cur_nd->x==dest->x && cur_nd->y==dest->y)
        {
            node* temp=cur_nd;
            vector<vector<int>> path;
            //path.push_back({cur_nd->x,cur_nd->y});
             
            while(temp!=start_nd) 
            {
              path.push_back({temp->x,temp->y,temp->state});
              temp=temp->par;
            }
            path.push_back({temp->x,temp->y,temp->state});
            path[0][2]=dest->state;		
            reverse(path.begin(),path.end());

            static_edge_path.insert(make_pair(make_pair(key(src_vertex),key(dest_vertex)),path)); 
            
            for(auto it:closed)
            delete it;   
            return(cur_nd->f);
        }

        long long int cur_node_k = key(cur_node_info);
        
        if(vis.find(cur_node_k)!=vis.end())
            continue;
            
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
  static_edge_path.insert(make_pair(make_pair(key(src_vertex),key(dest_vertex)),path));    	
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

	                vector<int> neigh_grid_cell{nbh.x,nbh.y};
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
      //copy_map(obs1,0); 
      for(auto it:allnodes)
      {
         it->par=NULL;
         it->f=FLT_MAX;
      }
      break;
    }

    if(visited.find(key(node_inf1))!=visited.end())
    continue;

    visited.insert(make_pair(key(node_inf1),true));
    
    // mark_obs(par_aut_state);
    for(int i=0;i<cur_nd->non_adj.size();i++)
    {  
          node_inf2.clear(); 
          node_inf2.push_back(cur_nd->non_adj[i]->y);
          node_inf2.push_back(cur_nd->non_adj[i]->x);
          node_inf2.push_back(cur_nd->non_adj[i]->state);
          float cost=static_edge_len[make_pair(key(node_inf1),key(node_inf2))];
          if(cost==INF)
            continue;         

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

    //copy_map(obs1,0); 
    if(!flag)flag=true;

  }
  for(auto it:allnodes)
  {
    it->par=NULL;
    it->f=FLT_MAX;
  }
   //copy_map(obs1,0); 
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
    
         for(auto i:cur->non_adj)
         {

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

}


void get_destn_nodes()
{
    int i;
    pair<int,int> dest_node;
    vector<cycle_node*> myqueue;
    map<long long int,bool>vis_node;

    for(int fs=0;fs<dest.size();fs++)
    {
        vector<int> vis_sy_state(pos_system_state.size());
        for(i=0;i<qsystate[dest[fs]].size();i++)
        {

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




void final_path_gen(int start)
{
  
  int tm=0;
  //cout<<"in final path gen .....path size is "<<temp_plan.size()<<endl;
  for(int i=0;i<=temp_plan.rbegin()->first;i++)
  {

    if(temp_plan.find(i)==temp_plan.end())
    {
      vector<int> path{temp_plan[tm][0],temp_plan[tm][1],temp_plan[tm][2]};
      final_plan.insert({start+i,path});

    }
    else
    {
      tm=i;
      vector<int> path{temp_plan[tm][0],temp_plan[tm][1],temp_plan[tm][2]};
      final_plan.insert({start+i,path});

    }      

  } 

}



float calc_dy_cost(cycle_node* source, cycle_node* dest,ull from,ull till=plan_till,bool gen_path=false)
{
    
    //cout<<" in calc_dy_cost() between "<<source->x<<","<<source->y<<" to "<<dest->x<<","<<dest->y<<endl;
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
    
    //cout<<"cur_ts and till is "<<cur_ts<<","<<till<<endl;
    while(!qopen.empty() && cur_ts<=till)
    {
        node* cur_nd = qopen.begin()->first;
        qopen.erase(qopen.begin());

        cur_ts=from+cur_nd->g;
        int current_automaton_state = cur_nd->state;
        vector<int> cur_node_info{cur_nd->y,cur_nd->x,current_automaton_state};
        //cout<<"popping "<<cur_node_info[1]<<","<<cur_node_info[0]<<","<<cur_node_info[2]<<endl;
        
        if(cur_nd->x==dest->x && cur_nd->y==dest->y)
        {
            //cout<<" found solution \n";
            if(gen_path)
            {
              temp_plan.clear();
              node* temp=cur_nd;
              temp_plan.insert({temp->g,{temp->x,temp->y,dest->state}});              
              while(temp!=start_nd)
              {
                //cout<<"inserting "<<temp->x<<","<<temp->y<<","<<temp->state<<endl;
                temp=temp->par;
                temp_plan.insert({temp->g,{temp->x,temp->y,temp->state}});

              }
                //temp_plan.insert({temp->g,{temp->x,temp->y,temp->state}});
                
                final_path_gen(from);
            }

            for(auto it:closed)
            delete it;
            return(cur_nd->f);
        }

        long long int cur_node_k = key(cur_node_info);
        
        if(vis.find(cur_node_k)!=vis.end())
            continue;
            vis[cur_node_k]=1;

        vector<int> neigh_state;

        point nbh;
        for(int i=0;i<dir;i++)
        {
            wait_time=0;
            nbh.y =  cur_nd->y+nb[i][0];
            nbh.x =  cur_nd->x+nb[i][1];

            if(!valid(cur_nd->y+nb[i][0],cur_nd->x+nb[i][1]) && ((cur_nd->y+nb[i][0])!=dest->y || (cur_nd->x+nb[i][1]!=dest->x)))
                continue;

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
   //cout<<"returning float max \n";
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


float ed_cost(cycle_node* src, cycle_node* dest,ull ts)
{

  vector<int>src_inf={src->y,src->x,src->state};
  vector<int>dest_inf{dest->y,dest->x,dest->state};
  pair<ull,ull> node_inf(key(src_inf),key(dest_inf));
  int i=0;
  for(i=0;i<dyn_edge_len[node_inf].size()-1;i++)
  {
    if(dyn_edge_len[node_inf][i].first<=ts && ts<dyn_edge_len[node_inf][i+1].first) 
     {
       return(dyn_edge_len[node_inf][i].second);
      
      }
   
  }
  return(dyn_edge_len[node_inf][i].second);

}





float get_dy_cy2(cycle_node* root,ull from,ull till)
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

    if(cur_ts>till)
    break;
    if(cur_nd==root && flag)
    {
       suf_dy_path.clear();
       suf_dy_path.push_back(make_pair(cur_nd->f,root));
       cycle_node* temp =root->par;
       while(temp!=root)
       {
         suf_dy_path.push_back(make_pair(temp->f,temp));
         temp=temp->par;
       }
       suf_dy_path.push_back(make_pair(0,root));
       reverse(suf_dy_path.begin(),suf_dy_path.end());

      //for(auto i:suf_dy_path)
        //cout<<i.second->x<<","<<i.second->y<<","<<i.second->state<<" at "<<i.first<<endl;

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
    

    for(int i=0;i<cur_nd->non_adj.size();i++)
    {

          node_inf2.clear();
	       float cost=ed_cost(cur_nd,cur_nd->non_adj[i],cur_ts);
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

    if(!flag)flag=true;
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


float get_dy_pref2(cycle_node* start, cycle_node* dest,ull from )
{

pwh_pref.clear();
if(start->x==dest->x && start->y==dest->y)
{
  pwh_pref.insert(make_pair(0,dest));
  return(0);
}

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
      
    for(int i=0;i<cur_nd->non_adj.size();i++)
    {

          node_inf2.clear();       
	        float cost=ed_cost(cur_nd,cur_nd->non_adj[i],cur_ts);
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
 
    visited.insert(make_pair(key(node_inf1),true));

  }
  
 for(auto it:allnodes)
 {
  it->par=NULL;
  it->f=FLT_MAX;
 } 
 return(FLT_MAX);
}







float get_dy_pref(cycle_node* start, cycle_node* dest,ull from )
{

pwh_pref.clear();
if(start->x==dest->x && start->y==dest->y)
{
  pwh_pref.insert(make_pair(0,dest));
  return(0);
}

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
    node_inf1.push_back(cur_nd->y);
    node_inf1.push_back(cur_nd->x);
    node_inf1.push_back(cur_nd->state);


    if(cur_nd==root && flag)
    {
       pwh_suf.clear();
       pwh_suf.insert(make_pair(root->f+prev_cost,root));
       cycle_node* temp =root->par;
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
       copy_map(obs1,0);   
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


void update_state(int tm)
{
  
  for(auto i=chng.begin();i!=chng.end();i++)
  {
    if(i->first==goal && i->second.second>0)
    {
      i->second.second--;
      i->first->rhs--;
      get_key(i->first);
      if(dyq.find(i->first)!=dyq.end())
       dyq.erase(i->first);
      dyq.insert(i->first);      
    
    }
    else if( i->second.second>0)
    { 
      i->second.second--;
      i->first->rhs--;
      get_key(i->first);
      if(dyq.find(i->first)!=dyq.end())
       dyq.erase(i->first);
       dyq.insert(i->first);

    }

  }
   	

}


void lla(ull cur_from,ull cur_till,cycle_node* cur, cycle_node* dest, pair<ull,ull> &var)
{


 int cnt=0;
 if(grid_dy.find(make_pair(cur->x,cur->y))!=grid_dy.end())
 { 
   int wait_time=grid_dy[make_pair(cur->x,cur->y)]-cur_from;
   vector<pair<ull,float>> temp;
   pair<ull,float> entry (cur_from,FLT_MAX);
   
   if(wait_time>0 && cur_from+wait_time<cur_till)
   {
 	 temp.push_back(entry);
	 dyn_edge_len.insert(make_pair(var,temp));
	 cur_from+=wait_time;

   }
   else if(wait_time>0)
   {
 	 temp.push_back(entry);
	 dyn_edge_len.insert(make_pair(var,temp));
	 return;
   }
 }  
 
 start = new_node(cur->x,cur->y,cur->state);
 start->rhs=0; 
 start->g=FLT_MAX;

 
 goal = new_node(dest->x,dest->y,dest->state); 
 goal->g=FLT_MAX;
 goal->h=0;
 goal->rhs=FLT_MAX;
 
 start->h=get_initial_heuristic(start,goal);
 get_key(start); 
 get_key(goal);
 
 bool flag=false;
 
 dyq.clear();
 chng.clear();
 unordered_map<int,node*> created;
 vector<int> vec(2);
 
 vec[0]=start->x;
 vec[1]=start->y;
 
 created[key(vec)]=start;
 
 vec[0]=dest->x;
 vec[1]=dest->y;

 
 created[key(vec)]=goal; 
 float best_cost=static_edge_len[var];
 int horizon=cur_till-cur_from;
 dyq.insert(start);
 
 while(cur_from<=cur_till)
 {

    while( dyq.size()!=0 && ((*dyq.begin())->k<goal->k || goal->rhs!=goal->g))
    {
      node* cur=*dyq.begin();
      dyq.erase(dyq.begin());

      if(cur->g > cur->rhs)
 	   cur->g=cur->rhs;
      if(cur_from + cur->g > cur_till)
      {
         pair<ull,float> entry (cur_from,FLT_MAX);
         dyn_edge_len[var].push_back(entry);
         return;
      }	   
      if(cur==goal)
      break;	 
 	    for(int i=1;i<dir;i++)
  	    {
    	      int x=cur->x+nb[i][0];
    	      int y=cur->y+nb[i][1];
              vector<int> nvec{x,y};

            if(!valid(y,x) && (y!=goal->y || (x!=goal->x)))
                continue;
   	      int wait=0;

    
	      if(grid_dy.find(make_pair(x,y))!=grid_dy.end())
      	        wait=max(0,(int)(grid_dy[make_pair(x,y)] - cur->rhs-cur_from));
      	      
   	    node* nd;
    	    if(created.find(key(nvec))!=created.end())
    	    {

		nd=created[key(nvec)];
		int g=cur->g+1+wait;

		if(nd->rhs>g)
		{
    	  	  nd->rhs=g;
		  nd->h=get_initial_heuristic(nd,goal);      	
		  if(dyq.find(nd)!=dyq.end())
		  dyq.erase(dyq.find(nd));
		  get_key(nd);
		  dyq.insert(nd);
		  if(wait>0) 
	          {
		    chng[nd]={cur->rhs+cur_from,wait};
	          }         		  
		  
    	        } 
             }	
    	     else	 
     	     {
  		nd=new_node(x,y,cur->state);
  		nd->h=get_initial_heuristic(nd,goal);
  		nd->rhs=cur->g+1+wait;
  		get_key(nd);

   		dyq.insert(nd);

   		created[key(nvec)]=nd;
		if(wait>0) 
		{
		  chng[nd]={cur->rhs+cur_from,wait};
		} 
    
             }   

     }

 }

 if(goal->rhs==best_cost)
 {
   pair<ull,float> entry (cur_from,best_cost);
   dyn_edge_len[var].push_back(entry);
   entry.first=cur_till-best_cost+1;
   entry.second=FLT_MAX;
   dyn_edge_len[var].push_back(entry);
   return;
 
 }
 else if(goal->rhs == FLT_MAX)
 {
   pair<ull,float> entry (cur_from,FLT_MAX);
   dyn_edge_len[var].push_back(entry);    
   return; 
 }
 else
 {
   pair<ull,float> entry (cur_from,goal->rhs);
   if(dyn_edge_len.find(var)==dyn_edge_len.end())
      dyn_edge_len[var].push_back(entry);    
   else
   {
      float prev_cost=(dyn_edge_len[var].end()-1)->second;
      if(prev_cost!=goal->rhs)
        dyn_edge_len[var].push_back(entry);            
      
    }   
  }
  cur_from++;
  update_state(cur_from);
 }

}


void get_dy_edge_cost(cycle_node* root,ull from,ull horizon)
{
   dyn_edge_len.clear();
   unordered_set<long long int> visited;
   map<long long int,cycle_node*> myq;
   
   vector<int> root_inf={root->y,root->x,root->state};
   myq.insert(make_pair(key(root_inf),root));
   while(!myq.empty())
   {
     cycle_node* cur=myq.begin()->second;
     myq.erase(myq.begin());
     vector<int>cur_inf={cur->y,cur->x,cur->state};

     if(visited.find(key(cur_inf))!=visited.end())
        continue; 

      else
      	visited.insert(key(cur_inf));
         
         int par_aut_state=cur->state;
    	 mark_obs(par_aut_state);  
         for(auto i:cur->non_adj)
         {
               vector<int>neigh_inf={i->y,i->x,i->state};
               pair<ull,ull> var(key(cur_inf),key(neigh_inf));
               float prev_cost=FLT_MAX;
           
               long long cur_from=from;
               long long cur_till=from+horizon;
               float best_cost=static_edge_len[make_pair(key(cur_inf),key(neigh_inf))]; 
               
               lla(cur_from,cur_till,cur,i,var);

               if(visited.find(key(neigh_inf))==visited.end())
               myq.insert(make_pair(key(neigh_inf),i));
         }    	 

         copy_map(obs1,0); 
     }
    
  
} 





bool is_dest(cycle_node* nd)
{

    for(auto i:graph_root)
       if(nd==i.second)
       return(true);
       
 return(false);

}






void get_model_file(cycle_node* new_origin,ull from,ull horizon)
{

str_to_cy.clear();
obj.clear();
cons.clear();
last_cons.clear();
gen_cons.clear();
all_var.clear();
int_pos.clear();
int_last_pos.clear();
map<int,unordered_set<cycle_node*> >myq;
last_cy.clear();


unordered_set<cycle_node*> qlist;
qlist.insert(new_origin);
myq.insert(make_pair(from,qlist));

ull ind_a=0;
ull ind_b=0;
ull till=from+horizon;
ull con_id=1;
int sos_id=1;
fstream fptr;
string prev_cons;

fptr.open("modelfile.txt",ios::out);
if(!fptr)
{
 cout<<"Open failed";
 return;
}

  
   while(!myq.empty())
   {
     unordered_set<cycle_node*> cur_nd = myq.begin()->second;
     ull time=myq.begin()->first;

     myq.erase(myq.begin());


     for(auto i:cur_nd)
     { 
       string s1,s2;
       vector<string> s3;
       unordered_set<string> tempset;
       s2="X_"+state_to_str(i)+"_"+to_string(time);
       str_to_cy[s2]=i;
       
       if(all_var.find(time)!=all_var.end())
       all_var[time].insert(s2);      
       else
       {
         tempset.clear();
         tempset.insert(s2);
         all_var.insert(make_pair(time,tempset));
       }
       
       if(is_dest(i))
       {
         string temp;
         string cy_cons;

         float cost=get_dy_cy2(i,time,till);

         if(cost!=FLT_MAX)
          {
             int icost=cost;
             s1="C_"+state_to_str(i)+"_"+to_string((int)time+icost)+"_"+to_string(icost);
             if(obj.find((int)time+icost)==obj.end())
             {
                unordered_set<string> fin;                
                fin.insert(s1);
                obj.insert(make_pair((int)time+icost,fin));
             }  
             else    
             obj[(int)time+icost].insert(s1);
             
             s1="X_"+state_to_str(i)+"_"+to_string(time+icost);   
             str_to_cy[s1]=i; 
               
             if(all_var.find(time+icost)!=all_var.end())
               all_var[time+icost].insert(s1);      
              else
              {
                tempset.clear();
         	tempset.insert(s1);
                all_var.insert(make_pair((time+icost),tempset));
             }

           cy_cons="(assert(= C_"+state_to_str(i)+"_"+to_string((int)time+icost)+"_"+to_string(icost)+" (and X_"+state_to_str(i)+"_"+to_string(time)+"  X_"+state_to_str(i)+"_"+to_string((int)time+icost)+" ) ) )";
           gen_cons.insert(cy_cons);

           
           temp="X_"+state_to_str(i)+"_"+to_string((int)time+icost);
	   if(myq.find((int)time+icost)==myq.end())
           {
             unordered_set<cycle_node*> list;
             list.insert(i);
             myq.insert(make_pair(((int)time+icost),list));
           }
           else
           myq[((int)time+icost)].insert(i);
           s3.push_back(temp);

          }
       
       }
       for(auto j:graph_root)
       {

         if(i==j.second)
         continue;

          float cost=get_dy_pref2(i,j.second,time);

          if(cost==FLT_MAX || (time+cost)>till)
          continue;
          
          int icost=cost;
          string temp="X_"+state_to_str(j.second)+"_"+to_string((int)time+icost);
          str_to_cy[temp]=j.second;
          
          if(myq.find((int)time+icost)==myq.end())
          {
            unordered_set<cycle_node*> list;
            list.insert(j.second);
            myq.insert(make_pair(((int)time+icost),list));
          }
          else
          myq[((int)time+icost)].insert(j.second);
          
          if(all_var.find(time+icost)!=all_var.end())
             all_var[time+icost].insert(temp);      
          else
           {
                tempset.clear();
         	tempset.insert(temp);
                all_var.insert(make_pair((time+icost),tempset));
           }


          s3.push_back(temp);
       
        }
        if(s3.size()!=0)
          cons.insert(make_pair(s2,s3));
        else
        last_cons.insert(s2);
       

    }
   
   
   }

  for(auto i:obj)
   for(auto j:i.second)
   fptr<<"(declare-const "<<j<<"  Bool)"<<endl;
  
  for(auto i:all_var)
  for(auto j:i.second)
   fptr<<"(declare-const "<<j<<"  Bool)"<<endl;
   
   for(auto i:last_cy)
   {
     for(auto j:i.second)
     fptr<<"(declare-const "<<j.second<<"  Int)"<<endl;
   
   }
   
 if(obj.size()!=0)
 {

   for(auto i:obj)
   {  
      map<int,string> temp;
      last_cy.insert(make_pair(i.first,temp));
      for(auto j:i.second)
      {
        unsigned pos=j.find_last_of("_");   
        unsigned len=stoi(j.substr(pos+1,j.size()-1));
        last_cy[i.first].insert(make_pair(len,j));
      }

   }

   
   
   fptr<<"(declare-const L_len Int)"<<endl;
   
   int count=2;
   fptr<<"(assert(= L_len ";

   for(auto i=last_cy.begin();i!=last_cy.end();i++)
   {
    for(auto j:i->second) 
    {
      count++;
      fptr<<"(ite(= "<<j.second<<" true)"<<j.first<<" ";
    }
   }

   fptr<<till<<" ";
   while(count!=0)
   {
     fptr<<" )";
     count--;
   }
   fptr<<endl;
  
 } 
   
   
   
   
   for(auto i:obj)
   {

     fptr<<"(declare-const P_"<<i.first<<"  Int)"<<endl;
     if(i.second.size()==1)
     fptr<<"(assert (= P_"<<i.first<<" (ite(= "<<*i.second.begin()<<" true) "<<i.first<<" "<<2*till<<" ) ) )\n";
     else
     {
       fptr<<"(assert (= P_"<<i.first<<"(ite(="<<" (or ";
       for(auto j:i.second)
       fptr<<j<<" ";
       fptr<<")true) "<<i.first<<" "<<2*till<<" ) ) )\n";
     }
   }

  
  
 fptr<<"(assert X_"+state_to_str(new_origin)+"_"+to_string(from)+" )\n";
 fptr<<"(declare-const T_max Int)"<<endl;
 if(obj.size()!=0)
 {
   int count=2;
   fptr<<"(assert(= T_max ";
   for(auto i=obj.rbegin();i!=obj.rend();i++)
   {
    count++;
    fptr<<"(ite(= P_"<<i->first<<" "<<i->first<<") "<<i->first<<" ";
   }
   fptr<<till<<" ";
   while(count!=0)
   {fptr<<" )";count--;}
   fptr<<endl;
 }
  
 
 for(auto i:cons)
  {
    string str="(assert(<= (+ ";
    if(i.second.size()==0)
    continue;
    fptr<<"(assert(=> "<<i.first<<" (or ";
    
    for(auto j=0;j<i.second.size();j++)
    {
     fptr<<i.second[j]<<" ";
     str.append(" (ite(= ");
     str.append(i.second[j]);
     str.append(" true) 1 0) ");
    }
    fptr<<") ) )\n";
    str.append(") 1 ) )");
    fptr<<str<<endl;
  } 
 
  for(auto i:cons)
  {

    string lhs,rhs;
    for(auto j=0;j<i.second.size();j++)
    {
      string and_str;
      string from,to;
      from=i.first;
      to=i.second[j];
      
      lhs="A_"+to_string(ind_a);
      rhs="B_";
     
        and_str="(assert(= "+lhs+" (and "+i.first+ "  "+i.second[j]+" ) ) )";
      
        unsigned start = from.find_last_of("_");
        unsigned start_time=stoi(from.substr(start+1,from.size()-1));
        unsigned end = to.find_last_of("_");
        unsigned end_time=stoi(to.substr(end+1,to.size()-1));
        
        if(end_time-start_time<=1)
        continue;
        rhs.append(to_string(start_time));
        rhs.append("_");
        rhs.append(to_string(end_time));
       
        if(int_pos.find(rhs)==int_pos.end())
        {
              unordered_set<string> temp;
              auto it=all_var.find(start_time);
              it++;
     	      for(;it!=all_var.end();it++)
              {  
                 if(it->first==end_time)break;
                 for(auto jt:it->second)
                   temp.insert(jt);      
              }
              int_pos.insert(make_pair(rhs,temp));
              if(temp.size()>0)
              {
                  fptr<<"(declare-const "<<lhs<<"  Bool)"<<endl;
                  fptr<<"(declare-const "<<rhs<<"  Bool)"<<endl;            
         	  string or_str="(assert(= "+rhs+" (or ";      
         	   for(auto i:int_pos[rhs])
          	   {     
          	    or_str.append(i);              
              	    or_str.append(" "); 
         	   }  
         	   or_str.append(") ) ) ");
          	   fptr<<or_str<<endl;
   		   fptr<<and_str<<endl;          	   
                   fptr<<"(assert(<= (+ (ite(= "<<lhs<<" true) 1 0) (ite(= "<<rhs<<" true) 1 0 ) ) 1 ) ) \n";                    
	           ind_a++;              
              }


         }            
         else if(int_pos[rhs].size()>0)
         {       
            fptr<<"(declare-const "<<lhs<<"  Bool)"<<endl; 
            fptr<<and_str<<endl;
            fptr<<"(assert(<= (+ (ite(= "<<lhs<<" true) 1 0) (ite(= "<<rhs<<" true) 1 0 ) ) 1 ) ) \n";                    
            ind_a++;
         }
      

      }
       
  }

  for(auto i:last_cons)
  {
  
     unsigned start = i.find_last_of("_");
     unsigned start_time=stoi(i.substr(start+1,i.size()-1));
     
     
     string and_str,or_str;
     string lhs,rhs;

     
     lhs="A_"+to_string(ind_a);
     rhs="E_"+to_string(start_time)+"_"+to_string(till);
    
     or_str="(assert(= "+rhs+" (or ";
     and_str="(assert(= "+lhs+" (and "+i+ "  "+i+" ) ) )";
     
     
     if(int_last_pos.find(rhs)==int_last_pos.end())
     {
         unordered_set<string> temp;
         auto it=all_var.find(start_time);
         it++;
     	 for(;it!=all_var.end();it++)
      	 {  
             for(auto jt:it->second)
                temp.insert(jt);      
          }
           int_last_pos.insert(make_pair(rhs,temp));
           if(temp.size()>0)
           {
               fptr<<"(declare-const "<<lhs<<"  Bool)"<<endl;
               fptr<<"(declare-const "<<rhs<<"  Bool)"<<endl;
               string or_str="(assert(= "+rhs+" (or ";      
               for(auto i:int_last_pos[rhs])
               {     or_str.append(i);              
                     or_str.append(" "); 

                }  
                or_str.append(") ) ) ");
                fptr<<and_str<<endl;
                fptr<<or_str<<endl;
                ind_a++;
                fptr<<"(assert(<= (+ (ite(= "<<lhs<<" true) 1 0) (ite(= "<<rhs<<" true) 1 0 ) ) 1 ) ) \n";                    

           }
      }
     else if(int_last_pos[rhs].size()>0)
     {
            fptr<<"(declare-const "<<lhs<<"  Bool)"<<endl; 
            fptr<<and_str<<endl;
            fptr<<"(assert(<= (+ (ite(= "<<lhs<<" true) 1 0) (ite(= "<<rhs<<" true) 1 0 ) ) 1 ) ) \n";                    
            ind_a++;

      }
            
  }   
  

  for(auto i:all_var)
   {
      
         if(i.second.size()<=1)
         continue;
         
         unordered_set<string> imap;
         unordered_map<string,unordered_set<string>> red;
         imap.insert(*(i.second.begin()));
         
         unsigned pos_temp = (i.second.begin()->find_last_of("_"));         
         string temp="Y_"+i.second.begin()->substr(2,pos_temp);
         unsigned pos=(temp.find_last_of("_"));
         temp=temp.substr(0,pos)+i.second.begin()->substr(pos_temp+1,i.second.begin()->size()-1);
         
         red.insert(make_pair(temp,imap)); 
         
           auto j=i.second.begin();
           j++;
 	   for(;j!=i.second.end();j++)
  	   {  

              unsigned pos_temp = (j->find_last_of("_"));
              string temp="Y_"+j->substr(2,pos_temp);
              unsigned pos=(temp.find_last_of("_")); 
              temp=temp.substr(0,pos)+j->substr(pos_temp+1,j->size()-1);

  	      if(red.find(temp)==red.end())
  	      {
	        unordered_set<string> myset;
                myset.insert(*j);
                red.insert(make_pair(temp,myset)); 
  	      }
             else
  	      red[temp].insert(*j);
  	   }
  	    auto k=red.begin();
  	    auto l=k->second.begin();
            string constr;
  	    if(k->second.size()==1)
            {   
               constr="(assert(<= (+ (ite(= "+*l+" true) 1 0) ";
               k++;    
            }	
  	    else
  	    { 

  	      fptr<<"(declare-const "<<k->first<<"  Bool)"<<endl;
              constr="(assert(<= (+ (ite(= "+k->first+" true) 1 0) ";

              string str="(assert(= "+k->first+" or( ";
  	      str.append(*l);
              l++;
              for(;l!=k->second.end();l++)
               {
                 str.append(*l);               
                 str.append(" ");

               }
               str.append(" ) ) )");
               gen_cons.insert(str);
  	       k++;
  	    
  	    }
  	   
  	   
            for(;k!=red.end();k++)
            { 
              if(k->second.size()<2)
              {
               constr.append(" (ite(= ");
               constr.append(*k->second.begin());
               constr.append(" true) 1 0 ) ");
               continue;
              }
              else
              {
               constr.append(" (ite(= ");
               constr.append(k->first);
               constr.append(" true) 1 0 ) ");
  	       fptr<<"(declare-const "<<k->first<<"  Bool)"<<endl;
  	       
               string str="(assert(= "+k->first+" or( ";
  	       str.append(*l);
               l++;
               for(;l!=k->second.end();l++)
               {
                 str.append(*l);               
                 str.append(" ");

               }
               str.append(" ) ) )");
               gen_cons.insert(str);
  	       k++;

              }
            }

         constr.append(" ) 1 ) )");
         fptr<<constr<<endl;
  	 con_id++;
     }
  
 for(auto i:gen_cons)
  fptr<<i<<endl;  
  
  

   
 // 1- *********OBJECTIVE***********  

 auto itr=obj.begin();
 if(obj.size()!=0)
 {
   fptr<<"(maximize(+";  
   auto inner=itr->second.begin();
   fptr<<" (ite(= "<<*inner<<" true) 1 0 )";
   inner++;
   for(;inner!=itr->second.end();inner++)
   fptr<<" (ite(= "<<*inner<<" true) 1 0 )";
   itr++;
   for(;itr!=obj.end();itr++)
     for(inner=itr->second.begin();inner!=itr->second.end();inner++)
          fptr<<" (ite(= "<<*inner<<" true) 1 0 )";
          fptr<<" ) )"<<endl;
 }
 
 
fptr<<"(minimize L_len)\n"; 
fptr<<"(minimize T_max)\n"; 
fptr<<"(check-sat)\n";
fptr<<"(get-model)\n";
 fptr.close();
}


void generate_path(cycle_node* src,cycle_node* dest, int time)
{
  
  if(src==dest)
  {
    //suffix cycle logic
    //cout<<"suffix cycle logic";
    get_dy_cy2(src,time,plan_till);
    auto itr1=suf_dy_path.begin();

      if(itr1==suf_dy_path.end())
      {
        //cout<<"returning";
        return;
      }

    auto itr2=std::next(itr1);
    while(itr2!=suf_dy_path.end())
    {
      calc_dy_cost(itr1->second,itr2->second,time+itr1->first,plan_till,true);
      itr1++;
      itr2++;

    }    

  }
  else
  {
   //prefix path logic
    //cout<<"Prefix path logic";

    get_dy_pref2(src,dest,time);
    auto itr1=pwh_pref.begin();

      if(itr1==pwh_pref.end())
      {
        //cout<<"returning";
        return;
      }


    auto itr2=std::next(itr1);
    while(itr2!=pwh_pref.end())
    {
      calc_dy_cost(itr1->second,itr2->second,time+itr1->first,plan_till,true);
      itr1++;
      itr2++;

    }

  }



}


pair<bool,int> plan_within_horizon(cycle_node* new_origin,ull from,ull len)
{

 get_model_file(new_origin,from,len);
 pwh_plan.clear();
 
 fstream fptr;
 fptr.open("model_soln.txt",ios::in | ios::ate);
 
 system("z3 modelfile.txt");

 fptr2.open("z3_op2.txt",ios::out | ios::app);
 if(!fptr)
 {
   cout<<"Cant open z3 output file\n";
   return(make_pair(false,INT_MAX));   
 }
 string word,final_cy;
 map<int,string> plan;
 int final_time=0;
 while(!fptr.eof())
 {
   fptr>>word;
   while(word!="|->")
   {
     if(word=="(error" || fptr.eof() || word=="Warning")
     {
       cout<<"returning because error1\n";
       return(make_pair(false,INT_MAX));
     }
     fptr>>word;        
   }
    fptr>>word;
    if(word=="0")
    {
      cout<<"returning because error2\n";
      return(make_pair(false,INT_MAX));
    }
   while(word!="T_max")
   {
     fptr>>word;
     if(word=="(error" || fptr.eof())
     {
      cout<<"returning because error3\n";      
      return(make_pair(false,INT_MAX));        
     }
   }
   fptr>>word;
   fptr2<<word<<endl;
   if(!fptr.eof()){fptr>>word;   
    fptr2<<word<<endl;}
   if(!fptr.eof() )
   final_time=stoi(word);
    while(!fptr.eof())
    {
      string var;
      int clk;
      fptr>>word;
      if(word[0]=='X')
      {
        var=word;
        while(word!="true)" && word!="false)")
        fptr>>word;
        if(word=="true)")
        {
          unsigned pos=var.find_last_of("_");
          clk=stoi(var.substr(pos+1,var.size()-1));
          if(clk>final_time)
          continue;
          //cout<<"inserting\n";
          pwh_plan.insert(make_pair(clk,str_to_cy[var]));
        }
      }

    }
 }
  fptr2<<"Final cy completion time is "<<final_time<<endl;
  for(auto i:pwh_plan)
  fptr2<<i.first<<"---"<<i.second->x<<","<<i.second->y<<","<<i.second->state<<endl;
  final_plan.clear();
  auto itr1=pwh_plan.begin();
  auto itr2=std::next(itr1);
  while(itr2!=pwh_plan.end())
  {
    generate_path(itr1->second,itr2->second,itr1->first);
    //fptr2<<i.first<<"---"<<i.second->x<<","<<i.second->y<<","<<i.second->state<<endl;
    itr1++;
    itr2++;
  }
  //for(auto i:final_plan)
   // cout<<i.first<<" --- "<<i.second[0]<<","<<i.second[1]<<","<<i.second[2]<<endl;
 fptr.close();
 fptr2.close();
 fptr2<<"returning true\n";
 return(make_pair(true,final_time));
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
       //copy_map(obs1,0); 
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
    
    //mark_obs(cur_nd->state);
    for(int i=0;i<cur_nd->non_adj.size();i++)
    {
          node_inf2.clear();       
          node_inf2.push_back(cur_nd->non_adj[i]->y);
          node_inf2.push_back(cur_nd->non_adj[i]->x);
          node_inf2.push_back(cur_nd->non_adj[i]->state); 

	        float cost=static_edge_len[make_pair(key(node_inf1),key(node_inf2))];
          if(cost==FLT_MAX)
            continue;      

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
    //copy_map(obs1,0); 
  }
  for(auto it:allnodes)
  {
   it->par=NULL;
   it->f=FLT_MAX;
  }
 return(FLT_MAX);

}



void put_new_cost(cycle_node* cur)
{

  vector<int>cur_inf={cur->y,cur->x,cur->state};
  for(auto i:cur->non_adj)
  {
              
       
       vector<int>neigh_inf={i->y,i->x,i->state};
       pair<ull,ull> temp_pair;
       temp_pair.first=key(cur_inf);
       temp_pair.second=key(neigh_inf);
       if(static_edge_len.find(temp_pair)!=static_edge_len.end())
       continue;
       
       float cost=calc_cost(cur,i);
      static_edge_len.insert(make_pair(make_pair(key(cur_inf),key(neigh_inf)),cost));


  
  }

}


double mx(double a,double b)
{
 if(a>b)
 return a;
 return b;

}


void remove_obs()
{
  while(global_ts>=dec_time)
  {

    for(auto it:dynamic_obstacles.begin()->second)
    {
      if(grid_dy.find(make_pair(it.first,it.second))!=grid_dy.end())
      {
        pair<int,int> xy(it.first,it.second);
        if(grid_dy[xy]<=global_ts) 
        grid_dy.erase(grid_dy.find(make_pair(it.first,it.second)));
      }
    }      
   dynamic_obstacles.erase(dynamic_obstacles.begin()); 
   dec_time=(dynamic_obstacles.size()==0?LONG_MAX:dynamic_obstacles.begin()->first);
  }
}



bool rdest(vector<int> &coords)
{
 
      for(auto i:graph_root)
       if(coords[0]==i.first.x && coords[1]==i.first.y && coords[2]==i.first.state)
       return(true);
       
      return(false);


}



int main(int args,char **argv)
{
    
    int x,y,num_obs;
    ull dur;
    long long plan_comp=LONG_MAX;
    float suff_best=FLT_MAX;    
    vector<pair<int,cycle_node*>> static_pref;
    bool plan_found=false;
    float cy_cost=FLT_MAX;
    float pref_cost=FLT_MAX;
    double t2=0.0;
    double t3=0.0;
    double mxt=0.0,t3m=0.0,t2m=0.0;
    fstream res;
    int ptime=1;
    int cycle_count=0;
    global_ts=0;
    
    res.open("result.txt",ios::out );
    if(!res)
    {
 	      cout<<"Open failed for result file\n";
        return(0);
    }
    
    if(args<4)
    {
    cout<<"Command line args needs to b properly passed (static file, dynamic file, horizon length and plan till)\n";
    return(0);
    
    }
    
    bool recheck=false;
    initializegrid(argv[1]);

    calsystates();
    create_graph(origin); 
    display_graph(origin);
    get_destn_nodes(); 
    mark_destn(origin);
   
    horizon_len=stoi(argv[3]);
    plan_till=stoi(argv[4]);
    get_static_edge_cost(origin);    

    for(auto i:graph_root)
     {
        get_min_cycle(i.second);
     }
    

    static_pref.clear();
    for(auto i:cy_created)
    i.second->f=FLT_MAX;
    for(auto i:all_suf_cycle)
    {
       float cost=get_pref_cost(origin,i.first);
       if(cost==FLT_MAX)
       continue;
       
       else if((i.second.first<suff_best)||(i.second.first==suff_best && cost<pref_cost))
       {
         pref_cost=cost;
         current_best=i.first;
         suff_best=i.second.first;
         static_pref.clear();
         
         for(auto j:pref_path)
         static_pref.push_back(j);
        
       }
       
    }
   
   if(static_pref.size()==0)
   {
    cout<<"No Destn is reachable .. Returning \n";
    return(0);  
   }   
   else
   {
     final_plan.clear();
     final_plan[0]={0,0,0};
     
     for(auto i=0;i<static_pref.size()-1;i++)
     {
        vector<int> inf1{static_pref[i].second->y,static_pref[i].second->x,static_pref[i].second->state};
        vector<int> inf2{static_pref[i+1].second->y,static_pref[i+1].second->x,static_pref[i+1].second->state};
        
        
        for(int j=1;j<static_edge_path[{key(inf1),key(inf2)}].size();j++)
           final_plan[ptime++]=static_edge_path[{key(inf1),key(inf2)}][j];

     }

   }
   
    reader2.open(argv[2]);  
    if(!reader2.eof())
    reader2>>obs_ts;
    if(obs_ts==-1)
    obs_ts=plan_till+1;
  
    check_time=min(obs_ts,plan_till);
     
     while(ptime<=check_time)
     { 
        for(int i=0;i<all_suf_cycle[current_best].second.size()-1;i++)
        {
          cycle_node* src=all_suf_cycle[current_best].second[i].second;
          cycle_node* dest=all_suf_cycle[current_best].second[i+1].second;

          vector<int> inf1{src->y,src->x,src->state};
          vector<int> inf2{dest->y,dest->x,dest->state};          
         
           for(int j=1;j<static_edge_path[{key(inf1),key(inf2)}].size();j++)
           final_plan[ptime++]=static_edge_path[{key(inf1),key(inf2)}][j];


        }

     }
      vector<int> prev_pos,next_pos;

      while(global_ts<plan_till)
      {
      	while(global_ts<plan_till && final_plan.find(global_ts)!=final_plan.end() && final_plan.find(global_ts+1)!=final_plan.end())
      	{
          //res<<"here!!!";
          prev_pos=final_plan[global_ts];
          next_pos=final_plan[global_ts+1]; 
          if(global_ts==check_time)
            break;

          if(global_ts==dec_time)
            remove_obs();

          if(rdest(next_pos))
          {
            //i.e. the position at which the robot is going to move is a destination node
            if(prev_final.size()==0)
              prev_final=next_pos;
            else
            {
             if(prev_final==next_pos)
             {
              //i.e the next position at which the robot is moving is the last destination node the robot visited 
              if(prev_pos!=next_pos)//Not waiting
                {
                  res<<"completed one cycle here \n";
                  cycle_count++;
                }
             }
             else
              prev_final=next_pos;//change the cycle origin node
            }

          }
         
          res<<prev_pos[0]<<" "<<prev_pos[1]<<endl;
          ++global_ts;
          //res<<"Moving from "<<prev_pos[0]<<","<<prev_pos[1]<<","<<prev_pos[2]<<" to "<<next_pos[0]<<","<<next_pos[1]<<","<<next_pos[2]<<"at ts="<<++global_ts<<endl;
          prev_pos=next_pos;
          //global_ts++;
          //Insert ros code here
        }



        if(global_ts>=plan_till)
        {
           res<<"1)Z3's number of cycle is "<<cycle_count<<endl;
           return(0);
        }
        //Replan logic
        vector<int> v{prev_pos[1],prev_pos[0],prev_pos[2]};
        if(cy_created.find(key(v))!=cy_created.end())
        {
          new_origin=cy_created[key(v)];
          new_origin->f=FLT_MAX;
          //res<<"new_origin is "<<new_origin->x<<","<<new_origin->y<<","<<new_origin->state<<endl;
          
        }
        else
        {
          new_origin=new_cycle_node(prev_pos[0],prev_pos[1],prev_pos[2]);
          cy_created.insert(make_pair(key(v),new_origin));
          //res<<"new_origin is "<<new_origin->x<<","<<new_origin->y<<","<<new_origin->state<<endl;

        }
        if(final_plan.find(global_ts)==final_plan.end() || final_plan.find(global_ts+1)==final_plan.end())
          recheck=true;
        if(global_ts==obs_ts)
        {
          recheck=true;
          reader2>>num_obs;

          for(int i=0;i<num_obs;i++)
          {
            reader2>>x;
            reader2>>y;
            reader2>>dur;
            ull max_time;

            vector<int> nd_inf={x,y};
            string s=conv_vec_to_string(nd_inf);

            if(x==new_origin->x && y==new_origin->y)
            {
              res<<"Algo cant mark current robot position as obstacle "<<x<<","<<y<<" Map is invalid "<<endl;
              return(0);
            }

            if(grid_dy.find(make_pair(x,y))!=grid_dy.end())
            {
               if(grid_dy[make_pair(x,y)]<(long long)(obs_ts+dur))
               {
                  grid_dy[make_pair(x,y)]=obs_ts+dur;
                  dynamic_obstacles[(obs_ts+dur)].push_back(make_pair(x,y));                   
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
           
          check_time=min(obs_ts,plan_till);           
          dec_time=(dynamic_obstacles.size()==0?LONG_MAX:dynamic_obstacles.begin()->first); 
        }
         if(global_ts==dec_time)
          remove_obs();
        //Repeatative SMT Solver Planning Logic

        if(recheck)
        {

        create_graph(new_origin);
        put_new_cost(new_origin);

        while(!plan_found)
        {
          global_ts+=tcomp;
          int tc=tcomp;
          // while(tc--)
          //   res<<new_origin->x<<" "<<new_origin->y<<endl;
          //res<<"Waiting at "<<new_origin->x<<" "<<new_origin->y<<new_origin->state<<endl;
          final_plan.clear();
          if(global_ts>=plan_till)
          {
            res<<" 2)Z3's  Number of cycle traversed = "<<cycle_count<<" global_ts "<<global_ts<<" plan till is "<<plan_till<<endl;
            return(0);
          }   

          while(obs_ts<=global_ts)
          {
            reader2>>num_obs;
            for(int i=0;i<num_obs;i++)
            {
              reader2>>x;
              reader2>>y;
              reader2>>dur;
              ull max_time;

              vector<int> nd_inf={x,y};
              string s=conv_vec_to_string(nd_inf);

              if(x==new_origin->x && y==new_origin->y)
              {
                res<<"Algo cant mark current robot position as obstacle "<<x<<","<<y<<" Map is invalid "<<endl;
                return(0);
              }

              if(grid_dy.find(make_pair(x,y))!=grid_dy.end())
              {
                if(grid_dy[make_pair(x,y)]<(long long)(obs_ts+dur))
                {
                  grid_dy[make_pair(x,y)]=obs_ts+dur;
                  dynamic_obstacles[(obs_ts+dur)].push_back(make_pair(x,y));                   
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
            check_time=min({obs_ts,plan_till});     
          }
          
          while(dec_time<=global_ts)
           remove_obs();



          long long len=min(horizon_len-1,plan_till-global_ts);        
          get_dy_edge_cost(new_origin,global_ts,len);

          pair<bool,int>temp=plan_within_horizon(new_origin,global_ts,len);
          plan_found=temp.first;
          plan_comp=temp.second;
          check_time=min({obs_ts,plan_till,plan_comp});



      }
      plan_found=false;
      recheck=false;

     }

   }
  // This must be the end of infinite planning
  res<<"3)Z3's  Number of cycle traversed = "<<cycle_count<<" till "<<global_ts<<endl;
  res.close();
  return(0);
}
