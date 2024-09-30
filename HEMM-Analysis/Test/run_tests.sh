#!/bin/bash

# Folder where Tests are defined
test_folder='./'

# Folder to store the matrices
matrices_folder=test_folder+'Matrices'
mkdir $matrices_folder

# Folder to store the results
results_folder=test_folder+'Results'
mkdir $results_folder

# Test parameters
dimensions=(1 2 4 8 16 32 64 128) # Matrix dimensions 
key_size=(128 256) # Key size
prec=(40) # Precision
random_matrices=false
reps=10 # Number of executions of each test (repetitions)

# Test matrix generation
# - If matrix files exist, nothing happens
# Matrix files are named [m1|m2]_[max_number_rows]x[max_number_columns].csv
# - If random_matrices, then random matrices will be generated and used at each execution
if ! $random_matrices
then
    cd $matrices_folder
    for dim in ${dimensions[@]}
    do
    	file_matrix1=m1_${dim}x${dim}.csv
    	if [ -f "$file_matrix1" ]
    	then
    		echo $file_matrix1 already existed
    	else
       		echo python3 $test_folder -s -fm1 $dim -cm1 $dim -sm1 $file_matrix1 -fm2 1 -cm2 1 >> generation.log
        	python3 $test_folder -s -fm1 $dim -cm1 $dim -sm1 $file_matrix1 -fm2 1 -cm2 1 >> generation.log
    	fi
    	file_matrix2=m2_${dim}x${dim}.csv
    	if [ -f "$file_matrix2" ]
    	then
    		echo $file_matrix2 already existed 
    	else
       		echo python3 $test_folder -s -fm1 1 -cm1 1 -fm2 $dim -cm2 $dim -sm2 $file_matrix2 >> generation.log
        	python3 $test_folder -s -fm1 1 -cm1 1 -fm2 $dim -cm2 $dim -sm2 $file_matrix2 >> generation.log
    	fi
    done
    cd ..
fi

# Launch tests: plain computation
cd $results_folder
for dim in ${dimensions[@]}
do
    for enc in ${key_size[@]}
    do
        for pre in ${prec[@]}
        do
            for i in $(seq 1 1 $reps)
            do
                    res_file='res_P'${dim}'_'${i}'.csv'
                    m1_file=$matrices_folder'/m1_'${dim}'x'${dim}'.csv'
                    m2_file=$matrices_folder'/m2_'${dim}'x'${dim}'.csv'
                    echo python3 $test_folder -s -P -csv $res_file -enc $enc -pre $pre -fm1 $dim -cm1 $dim -fm2 $dim -cm2 $dim -lm1 $m1_file -lm2 $m2_file >> execution.log
                    python3 $test_folder -s -P -csv $res_file -enc $enc -pre $pre -fm1 $dim -cm1 $dim -fm2 $dim -cm2 $dim -lm1 $m1_file -lm2 $m2_file
            done
        done
    done
done
cd $test_folder

# Launch tests: EVA and SEAL
cd $results_folder
for dim in ${dimensions[@]}
do
    for enc in ${key_size[@]}
    do
        for pre in ${prec[@]}
        do
            for i in $(seq 1 1 $reps)
            do
                    res_file='res_'${enc}'_'${pre}'_'${dim}'_'${i}'.csv'
                    m1_file=$matrices_folder'/m1_'${dim}'x'${dim}'.csv'
                    m2_file=$matrices_folder'/m2_'${dim}'x'${dim}'.csv'
                    echo python3 $test_folder -s -E -S -csv $res_file -enc $enc -pre $pre -fm1 $dim -cm1 $dim -fm2 $dim -cm2 $dim -lm1 $m1_file -lm2 $m2_file >> execution.log
                    python3 $test_folder -s -E -S -csv $res_file -enc $enc -pre $pre -fm1 $dim -cm1 $dim -fm2 $dim -cm2 $dim -lm1 $m1_file -lm2 $m2_file
            done
        done
    done
done
cd $test_folder
