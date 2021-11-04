# Build the FMU 
num_of_tanks=$1
num_of_replications=$2
omc sim_tank.mos
# Cleanup the junk produced by the omc compiler
rm -f Tank_*
rm -f Tank.c
rm -f Tank.h
# Unpack the fmu data 
unzip -o -qq Tank.fmu
rm -f Tank.fmu
rm -f *.json
rm -f *.libs
rm -f *.log
rm -rf sources 
mv *.xml Tank
# Compile our simulator. You must have in the include path
# The FMI for model exchange header files
# You must link to system library libdl
g++ -Wall -fopenmp -I../../../../include -I${HOME}/Documents/FMI/FMI_for_ModelExchange_and_CoSimulation_v2.0 main_par.cpp -ldl
# Run the model
./a.out $num_of_tanks $num_of_replications


