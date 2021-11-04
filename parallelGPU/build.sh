# Build the FMU 
num_of_tanks=$1
num_of_replications=$2
num_of_cpus=$3	
# Compile our simulator. You must have in the include path
#    The FMI for model exchange header files
# You must link to system library libdl
g++ -Wall -fopenmp -I../../../../include main_par.cpp -std=c++11 -fpermissive -ldl
# Run the model
#for j in {4..4}
#	do
	for i in {1..1}
		do	
		./a.out $num_of_tanks $num_of_replications $num_of_cpus
		done
	#done
