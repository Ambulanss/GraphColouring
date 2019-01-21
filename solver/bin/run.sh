# to run tests, run this script in bin directory
# you can skip one, two, three or four last arguments

# pattern:

#./solver ../instances/INSTANCE_NAME POPULATION_SIZE BREEDERS_Q BREEDERS_P NUMBER_OF_GENERATIONS 

printf "\n[INSTANCE 0]\n"
./solver ../instances/simple.txt 100 0.4 0.1 100
printf "\n[INSTANCE 1]:\n"
./solver ../instances/gc500.txt 100 0.4 0.1 50 # 83 colours!
