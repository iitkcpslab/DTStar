#!/bin/sh

rm -f result.txt
#g++ env_gen.cpp -std=c++11 -o env_gen


echo "file is $file"
for i in `seq 1 100`; 
do
         echo "iteration $i"
        ./env_gen $1 

	./z3_opt $3 dy.txt $2 $1 > model_soln.txt
        ./greedy1 $3 dy.txt $1
        ./greedy2 $3 dy.txt $1	
        #./z3_ll_123 $3 dy.txt 70 $1 > ll_op.txt
	#./z3_ll_123 $3 dy.txt 100 $1 > ll_op.txt
	#./z3_ll_123 $3 dy.txt 120 $1 > ll_op.txt
	#./z3_ll_123 $3 dy.txt 150 $1 > ll_op.txt			        				        

done		        				        

done
