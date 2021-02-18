# DTStar

REQUIREMENTS:
ltl2tgba tools , Z3 solver

Installation for spot ltl2ba library:
https://spot.lrde.epita.fr/install.html

Install Spot-2.6/2.8 tool for LTL2TGBA converter
Check if the tool is working by running example command in the current folder:
./ltl2tgba --spin '<>p && <>q'
Copy the ltl2tgba file from the bin folder of spot installation to the folder with z3_opt.cpp 



Files in the working directory and syntax of input
......................................................
1)query.dat --> LTL query command for LTL to BA converter 
Syntax :
./ltl2tgba --spin '[]( <>p1 & <>p2) & []( p2 -> X(!p2 U p1)) & []( p1 -> X(!p1 U p2)) '


2)Workspace files:
static.txt --> file containing the static environment information:
Syntax:
nr - number of rows
nc - number of columns
nobs - no of obstacles followed by coordinates of obstacles
np - no of coordinates where a certain proposition is true followed by coordinates and the proposition true at it.

Example file:
10
10
5
1 2
3 4
5 6
9 5
3 9
4
5 5 1
6 6 2
2 2 3 
3 3 4


3)Workspace files:
dynamic.txt--> file containing the dynamic environment information:
Syntax:
ts1- timestamp at which obstacles appears
nobs - no of obstacles followed by coordinates of obstacles and duration for which these obstacles will stay
ts2-  next timestamp at which obstacles appears
nobs - no of obstacles followed by coordinates of obstacles 
.
.
.
-1

Example file:
201  --> At time stamp 201
2    --> 2 locations are becoming unreachalbe
3 18 32  -->(3,18) unreachable till 201+32=233
10 12 38 -->(10,12) unreachable till 201+38=239
313
2
10 12 31
3 18 24
408
4
3 18 27
11 18 39
1 5 28
6 15 36
-1  -->no more obstacles


To compile :
............
1)g++ z3_opt -std=c++11 -o z3_opt
2)g++ greedy1.cpp -std=c++11 -o greedy1
3)g++ greedy2.cpp -std=c++11 -o greedy2


To run :
.........
1)./z3_opt static.txt dynamic.txt Horizon_length Total_planning_time |> model_soln.txt
2)./greedy1 static.txt dynamic.txt Total_planning_time
3)./greedy2 static.txt dynamic.txt Total_planning_time

eg:
1)./z3_opt static dynamic 50 500|> model_soln.txt
2)./greedy1 static dynamic 500
3).//greedy2 static dynamic 500

This means that based on the query given in query.dat, environment in static file and dynamic changes mentioned in dynamic.txt, the z3 generates a solution given in "result.txt". "Model_soln.txt" will contain the model solution for each plan_in_horizon call. Horizon length for above example is 200 timestamps and total planning time is 500 timestamps. For greedy algorithms we dont need horizon length.


OUTPUT FORMAT : 
The results are stored in "result.txt" file. One could check the decision sequence and the final number of cycles travesed till Total_planning_time in this file for all the three algorithms. 

© 2020 GitHub, Inc.
Terms
Privacy
Security
Status
Help
Contact GitHub
Pricing
API
Training
Blog
