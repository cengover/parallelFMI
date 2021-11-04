# Build the FMU 
num_of_cpus=$1
# Compile our simulator. You must have in the include path
# The FMI for model exchange header files
# You must link to system library libdl
g++ -Wall -O -fopenmp -I${HOME}/Documents/FMI/FMI_for_ModelExchange_and_CoSimulation_v2.0 main.cpp -fpermissive -std=c++11 -ldl -o run
# Run the model
# The number of cores
for i in {1..4}
	do
	# The number of replications
	for j in {1..10}
		do
		./run $i
		done
	echo " DONE with $i"
	done
