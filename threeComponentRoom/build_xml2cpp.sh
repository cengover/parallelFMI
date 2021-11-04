# Clean up stale files
rm -f a.out
rm -rf binaries
# Cleanup the junk produced by the omc compiler - The file name (ParallelComponents_Components_X) should be replaced with the fmu file name. 
rm -f ParallelComponents_Components_Controller_*
rm -f ParallelComponents_Components_Controller.c
rm -f ParallelComponents_Components_CoolingSystem_*
rm -f ParallelComponents_Components_CoolingSystem.c
rm -f ParallelComponents_Components_Outside_*
rm -f ParallelComponents_Components_Outside.c
rm -f ParallelComponents_Components_Zone_*
rm -f ParallelComponents_Components_Zone.c
# Unpack FMUs and generate header files
unzip -o -qq ParallelComponents_Components_Controller.fmu
python xml2cpp.py -r modelDescription.xml -type double -f $PWD/binaries/linux64/ParallelComponents_Components_Controller.so -o Controller
unzip -o -qq ParallelComponents_Components_CoolingSystem.fmu
python xml2cpp.py -r modelDescription.xml -type double -f $PWD/binaries/linux64/ParallelComponents_Components_CoolingSystem.so -o Cooling
unzip -o -qq ParallelComponents_Components_Zone.fmu
python xml2cpp.py -r modelDescription.xml -type double -f $PWD/binaries/linux64/ParallelComponents_Components_Zone.so -o Zone
