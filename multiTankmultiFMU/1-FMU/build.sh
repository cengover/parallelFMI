# Build the FMU 
num_of_tanks=$1
num_of_replications=$2
# Change the .mo file here!
sed "s/TanksInSeries TT(n=5);/TanksInSeries TT(n=$num_of_tanks);/g" TankMParse.mo > TankM.mo
# Build the FMU 
omc sim_tankM.mos
rm -f SerialTanks_*
rm -f SerialTanks.c
rm -f SerialTanks.h
# Unpack the fmu data 
unzip -o -qq SerialTanks.fmu
#rm -f SerialTanks.fmu
rm -f *.json
rm -rf sources 
mv *.xml SerialTanks	
# Compile our simulator. You must have in the include path
#    The FMI for model exchange header files
# You must link to system library libdl
g++ -Wall -fopenmp -I../../../../include -I${HOME}/Documents/FMI/FMI_for_ModelExchange_and_CoSimulation_v2.0 main_par.cpp -std=c++11 -fpermissive -ldl
# Run the model
./a.out $num_of_tanks $num_of_replications

