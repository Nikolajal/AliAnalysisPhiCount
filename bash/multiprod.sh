for i in {1..10}
do
   ./exe/GeneratorMC_BKG a0.root >a0.log &
done

#  $ROOT_SYS hadd outGeneratorMC.root *.root
