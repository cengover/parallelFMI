# Clean up stale files
rm -f a.out
rm -rf binaries
# Cleanup the junk produced by the omc compiler - The file name (ParallelComponents_Components_X) should be replaced with the fmu file name. 
rm -f ParallelComponents_Components_OneZoneOneCooler1_*
rm -f ParallelComponents_Components_OneZoneOneCooler1.c
rm -f ParallelComponents_Components_OneZoneOneCooler2_*
rm -f ParallelComponents_Components_OneZoneOneCooler2.c
rm -f ParallelComponents_Components_OneZoneOneCooler3_*
rm -f ParallelComponents_Components_OneZoneOneCooler3.c
rm -f ParallelComponents_Components_OneZoneOneCooler4_*
rm -f ParallelComponents_Components_OneZoneOneCooler4.c
# Unpack FMUs and generate header files
unzip -o -qq ParallelComponents_Components_OneZoneOneCooler1.fmu
python xml2cpp.py -r modelDescription.xml -type double -f $PWD/binaries/linux64/ParallelComponents_Components_OneZoneOneCooler1.so -o ParallelComponentsComponentsOneZoneOneCooler1
unzip -o -qq ParallelComponents_Components_OneZoneOneCooler2.fmu
python xml2cpp.py -r modelDescription.xml -type double -f $PWD/binaries/linux64/ParallelComponents_Components_OneZoneOneCooler2.so -o ParallelComponentsComponentsOneZoneOneCooler2
unzip -o -qq ParallelComponents_Components_OneZoneOneCooler3.fmu
python xml2cpp.py -r modelDescription.xml -type double -f $PWD/binaries/linux64/ParallelComponents_Components_OneZoneOneCooler3.so -o ParallelComponentsComponentsOneZoneOneCooler3
unzip -o -qq ParallelComponents_Components_OneZoneOneCooler4.fmu
python xml2cpp.py -r modelDescription.xml -type double -f $PWD/binaries/linux64/ParallelComponents_Components_OneZoneOneCooler4.so -o ParallelComponentsComponentsOneZoneOneCooler4

