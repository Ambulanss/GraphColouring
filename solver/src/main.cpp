#include <iostream>
#include <vector>
#include <bits/stdc++.h> 
#include "graphloading.hpp"
#include <algorithm>
#include <random>

int lowest_colour(std::vector<bool> &row, std::vector<int> &result){
    /* -----
    Returns the lowest colour for colouring given vertex.

    @row - row from adjacency matrix, containing neighbours of vertex that is being processed
    @result - vector containing current number of used colours (at index 0)
            and colours of vertices (index means a vertex number and value is colour)
    ----- */

    std::vector<bool> available_colours (result[0]+2, true); //vector of all colours used so far + 1
    for (unsigned int i = 1; i < row.size(); i++){
        if (row[i]){ //if a neighbour is found
            if (result[i]){ //if a neighbour is coloured - his value in result vector is != 0
                available_colours[result[i]] = false; //set this colour as impossible to use
            }
        }
    }
    for (unsigned int i = 1; i < available_colours.size(); i++){
        if (available_colours[i]){ //if a colour is possible to use
            if (i == available_colours.size()-1){ //if chosen colour is new
                result[0]++; //then the number of colours must be incremented
            }
            return i; 
        }
    }
    return -1;
}


std::vector<int> greedy(std::vector<std::vector<bool> > &matrix){
    /* -----
    Returns a vector with number of used colours as the first element
    and colours of consecutive vertices.

    @matrix - adjacency matrix of input graph
    ----- */

    int colour;
    std::vector<int> result (matrix.size(), 0);
    for (unsigned int i = 1; i < matrix.size(); i++){
        colour = lowest_colour(matrix[i], result);
        result[i] = colour;
    }
    return result;
}


void print_result(std::vector<int> vec){
    std::cout << "\nNumber of colours: " << vec[0] << std::endl;
    std::cout << "Consecutive colour values: " << std::endl;
    for (unsigned int i = 1; i < vec.size(); i++)
        if (i != vec.size()-1)
            std::cout << vec[i] << ", ";
        else
            std::cout << vec[i] << std::endl; 
}


bool sortbysec(const std::pair<int,int> &a, const std::pair<int,int> &b) {   
    /* ----- 
    Returns the third argument for sort() function
    - for sorting a vector of pairs by second element of pairs.
    ----- */
    return (a.second < b.second); 
} 


bool sortbysecdesc(const std::pair<int,int> &a, const std::pair<int,int> &b) { 
    /* -----
    Returns the third argument for sort() function - for sorting in descending order
    a vector of pairs by second element of pairs.
    ----- */
        return a.second>b.second; 
} 



std::vector< std::pair <int,int> > sort_degrees(std::vector<std::vector<bool> > &matrix, 
bool increasing = true){
    /* -----
    Returns a sorted vector of pairs, where every pair consists of a vertex and its degree.
    A vector is sorted by a vertex's degree.

    @matrix - adjacency matrix of input graph
    @increasing - if true, function sorts by descending order, if false - in descending
    ----- */
    unsigned n = matrix.size();
    std::vector< std::pair <int,int> > result (n);
    unsigned degree = 0;
    for (unsigned int i = 1; i < n; i++){
        result[i].first = i;
        for (unsigned int j = 1; j < n; j++){
            if (i != j && matrix[i][j]){ //if a neighbour is found
                degree++; 
            }
        }
        result[i].second = degree;
        degree = 0;
    }
    if (increasing){
        sort(result.begin(), result.end(), sortbysecdesc);
    }
    else {
        sort(result.begin(), result.end(), sortbysec); 
    }
    return result;
}


std::vector<int> largest_first(std::vector<std::vector<bool> > &matrix){
    /* -----
    Returns a vector with number of used colours as the first element
    and colours of consecutive vertices. Uses largest first algorithm - 
    processes vertices with highest degree first.

    @matrix - adjacency matrix of input graph
    ----- */
    unsigned n = matrix.size();
    int colour;
    std::vector<int> result (n, 0);
    std::vector< std::pair<int, int> > indexes (n);
    indexes = sort_degrees(matrix);

    for (unsigned int i = 1; i < n; i++){
        colour = lowest_colour(matrix[indexes[i-1].first], result); //indexes is indexed from 0
        result[indexes[i-1].first] = colour;
    }
    return result;
}

auto parse_input(int argc, char const *argv[])
{
    std::vector<std::vector<bool> > matrix;
    std::string filename;
    if(argc > 1){
        filename = static_cast<std::string>(argv[1]);
        matrix = load_matrix(filename);
        return matrix; 
    }
    else
    {
        std::cout<<"Usage:\nsolver <instance_name>"<<std::endl;
        exit(-1);
    }
}


class genetic_parameters{
    public:
        int population_size;
        double breeders_p; //percent of all individuals; parents chosen as the best of all - TODO: it has to be double
        double breeders_q; //percent of all individuals; parents chosen randomly - TODO: it has to be double
        double mutation_chance; //chance for mutation, expressed as percent - TODO: it has to be double
        /* parents_choosing - number expressing when should we use 
        first parents choosing method and when the second one. 
        It means after which generation should we switch to the second method. */
        int parents_choosing; 
        int number_of_generations; //number of generations before stopping the algorithm
        int create_child_method; //if 1, use first mutation method, if 2, use the second one
        double mutation_size; //percent of values in child that should be mutated; expressed as a fraction

        genetic_parameters(){}

        genetic_parameters (int pop_size, double p = 0.25, double q = 0.25, double mut_chance = 0.0, 
                            int par_choosing = 1, int number_of_gen = 100, int breeding_method = 1, 
                            double mut_size = 0.0){
                population_size = pop_size;
                breeders_p = p;
                breeders_q = q;
                mutation_chance = mut_chance;
                parents_choosing = par_choosing;
                number_of_generations = number_of_gen;
                create_child_method = breeding_method;
                mutation_size = mut_size;
            }
};


class individual {
        /* value vector keeps the order of vertices to process by typical greedy algorithm.
        We iterate not by i = 1, 2, 3, ..., n-1, but by values of value vector.
        For example, value[1] = 3 means that the first vertex we should colour is vertex nr 3; 
        value[2] = 7 means that the second vertex we should colour is vertex nr 7; and so on... */

    public: 
        std::vector <int> value;            //vertices
        int fitness; 
        int n;                              //size of main adjacency matrix
        std::vector <int> result_colours;   //result_colours[0] is fitness - after running count_fitness

        individual(){}                      //default constructor

        individual (int my_n){ //engine cannot be passed as a constructor parameter;
                                                                //it has to be only shuffle method argument
            n = my_n;
            fitness = 0;
            //initialize_value(); - I left it to remember about running this method when generating population 
            //from external class - it shouldn't be always ran by a constructor
            initialize_colours();
            //generate(); - like above
        }

        void initialize_colours (){ //fill colours vector with n zeroes
            for (int i = 0; i < n; i++){
                result_colours.push_back(0);
            }
        }

        void initialize_value (){ //fill value vector with n zeroes
            for (int i = 0; i < n; i++){
                value.push_back(0);
            }
        }
        //std::default_random_engine rng
        void generate() { //generate random permutation - random vertices' order
            for (int i = 1; i < n; i++){ //first assign value equal to index to later shuffle it
                value[i] = i;
            }
            std::random_shuffle(std::begin(value)+1, std::end(value)); 
            /* shuffle everything apart from the first element
            because the first element is dedicated for/to number of colours (which version is gramatically correct?)
            and is initially equal to 0 */

        }

        void count_fitness (std::vector<std::vector<bool> > &matrix){
            /* This function is the proper algorithm of graph colouring, which uses vertex order
            coded in value vector.

            fitness value is kept as result_colours[0] - it's the number of colours used for graph colouring. */
            
            int colour;
            for (int i = 1; i < n; i++){ 
                colour = lowest_colour(matrix[value[i]], result_colours); 
                result_colours[value[i]] = colour;
            }
            fitness = result_colours[0];
        }

        void print_individual() { //every class should have something like this for debugging reasons I think
            std::cout << "\nIndividual: ";
            for (int i = 1; i < n; i++){
                std::cout << value[i] << " ";
            }
            std::cout << "\nFitness: " << result_colours[0] << "\n";
        }
};

void display_population (std::vector< individual > population, int pop_size, int n){
    for (int i = 0; i < pop_size; i++){
        std::cout << i << ".      Fitness: " << population[i].fitness << std::endl;
        std::cout << "Order:  ";
        for (int j = 0; j < n; j++){
            std::cout << population[i].value[j] << " ";
        }
        std::cout << "\nColour: ";
        for (int j = 0; j < n; j++){
            std::cout << population[i].result_colours[j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}



class parents_pair{
    /* one pair of parents */
    public:
        individual parent_A;
        individual parent_B;
        individual child;
        int n;

        parents_pair(individual A, individual B, int input_n): parent_A(A), parent_B(B), n(input_n){ 
        /* individual needs two parameters for constructor. Random engine cannot be one of them (it has to be changed),
        because it makes a mess. I would leave only n. Because of that, we need to pass n as a constructor
        parameter for parents_pair class I guess, and then set parent_A and parent_B parameters in function
        like "set_parents". */ 
        //child = individual(n);
        }

        void create_child1(std::default_random_engine rng){
            /* We take first values from parent_A (until crosspoint index) and the rest
            from parent_B (the ones that are possible to use without repeating). 
            We then will the gaps with random values chosen from available ones.

            I would add functionality to randomly use parent_A's values at the beginning 
            of the child (before crosspoint) or at the end (after crosspoint).
            */
            std::cout << "Create child\n";
            child = individual(n);
            child.initialize_value(); //now the child vector is filled with zeroes
            std::vector <bool> available_values (child.n, true); //should create vector with n true values
            std::cout << "available values: ";
            for (int i = 0; i < n; i++){
                std::cout << available_values[i] << " ";
            }
            std::cout << std::endl;
            available_values[0] = false; //we don't have vertex nr 0
            
            std::vector <int> left_indexes; //we will use them to know where the gaps are
            bool fill_the_gaps = false; //if this is true, we know that there are some gaps (it will fore sure always be true)
            
            int crosspoint;
            std::random_device rd;
            std::mt19937 mt(rd());
            auto dist = std::uniform_int_distribution<int>(1, child.n-1);
            std::cout << "child.n: " << child.n << std::endl;
            crosspoint = dist(mt);
            std::cout << "Crosspoint: " << crosspoint << std::endl;
            std::cout << "Parent_A values: ";

            for (int i = 0; i < parent_A.value.size(); i++)
                std::cout << parent_A.value[i] << " ";
            std::cout << "\nParent_B values: ";
            for (int i = 0; i < parent_B.value.size(); i++)
                std::cout << parent_B.value[i] << " ";
            std::cout << std::endl;

            for (int i = 1; i < crosspoint; i++){ //we take first values from parent A
                child.value[i] = parent_A.value[i];
                available_values[parent_A.value[i]] = false; //we mark used vertex as impossible to use
                std::cout << "Taking from parent A value " << parent_A.value[i] << " and put it on i = " << i << std::endl;
                /* reminder: value field contains vertices' numbers, so parent_A.value[2] = 3 means
                that parent_A is vertex order in which the second vertex processed by greedy algorithm
                is the third vertex in the input graph */
            }
            for (int i = crosspoint; i < child.n; i++){ //the rest is filled with possible values from parent_B
                if (available_values[parent_B.value[i]] == true){ //if vertex hasn't been used yet
                    child.value[i] = parent_B.value[i];    //then insert it into child
                    available_values[parent_B.value[i]] = false;
                    std::cout << "Taking from parent B value " << parent_B.value[i] << " and put it on i = " << i << std::endl;
                }
                else {                                  //remember that you have to fill this gap later
                    left_indexes.push_back(i);          //index i has no vertex
                    fill_the_gaps = true;
                } 
            }
            
            if (fill_the_gaps){ //if there are some gaps to fill

                std::cout << "left indexes: ";
                for (int i = 0; i < left_indexes.size(); i++)
                    std::cout << left_indexes[i] << " ";
                std::cout << std::endl;
            
                std::cout << "Child before filling gaps:\n";
                for (int i = 0; i < child.value.size(); i++){
                    std::cout << child.value[i] << " ";
                }
                std::cout << std::endl;
                std::random_shuffle(std::begin(left_indexes), std::end(left_indexes));
                for (int i = 1; i < child.n; i++){
                    std::cout << "Im in fatal for.\n";
                    if(available_values[i]){ //if unused vertex is found
                        std::cout << "I found unused value = " << i << std::endl;
                        child.value[left_indexes[0]] = i; 
                        std::cout << "child.value[left_indexes[0]] = available_values[i] \n";
                        std::cout << "left_indexes[0] = " << left_indexes[0] << std::endl;
                        std::cout << "child.value[left_indexes[0]] = " << child.value[left_indexes[0]] << std::endl;
                        /* assign first available vertex to random gap in child's values - random, because 
                        free left indexes were shuffled. */
                        available_values[i] = false; //set this vertex as used
                        std::cout << "Now erase first element from left indexes, left it like: \n";
                        left_indexes.erase(std::begin(left_indexes));
                        for (int i = 0; i < left_indexes.size(); i++){
                            std::cout << left_indexes[i] << " ";
                        }
                        std::cout << std::endl;
                    }
                }
            }
            for (int i = 0; i < n; i++){
                std::cout << child.value[i] << " ";
            }
            std::cout << std::endl;
        }

        void create_child2(){
            /* for every position randomly take value 
            from first or second parent;
            if there are repeated values, take the nearest;
            
            maybe there sould be a parameter like 
            int = {1, 2, 3}, saying that you can 
            insert not one, but two or three values
            at once?
            */
        }

        void mutate_child1(double mutation_size, std::default_random_engine rng){
            /* Function swaps two random values in child. 
            We don't use it yet. */
            int mutation_number = static_cast<int>(floor(((child.n-1) * mutation_size))); 
            auto dist = std::uniform_int_distribution<int>(1, child.n);
            for (int i = 0; i < mutation_number; i++){
                random1 = dist(rng);
                random2 = dist(rng);
                std::cout <<"[DEBUG] Random numbers in mutate_child1: " << random1 << " " << random2<<"\n";
                swap(random1, random2);
            } 
        }

        void mutate_child2(){
            //maybe later
        }

    private:
        int random1;
        int random2;

        void swap(int a, int b){
            /* Function swaps two values in child (swaps order of two vertices) */
            int temp;
            temp = child.value[a];
            child.value[a] = child.value[b];
            child.value[b] = temp;
        }

};

class parents {
    public:
        std::vector <individual> sorted_population; //population sorted by fitness, ascending
        std::vector <individual> all_parents;
        int breeders_N;                             //number of parents to choose as the best individuals
        int breeders_M;                             //number of parents to choose as random individuals
        std::default_random_engine rng;
        int parents_number;

        parents(){}

        parents(double p, double q, std::vector <individual> input_population, std::default_random_engine my_rng){
            /* Constructor initializing object keeping all parents in population.
            Parameters:
            @p - percent of parents to choose as the best individuals --- it has to be double ;_; 
                                                                        //it's treated as fraction later
            @q - percent of parents to choose as random individuals --- just like above
            @input_population - sorted population vector */
            rng = my_rng;
            sorted_population = input_population;
            if (p + q > 100){
                std::cout << "Error: wrong number of breeders. \n";
                exit(-1);
            }
            breeders_N = static_cast<int>(floor(p*sorted_population.size()));
            std::cout << "Breeders_N: " << breeders_N << std::endl;
            breeders_M = static_cast<int>(floor(q*sorted_population.size()));
            std::cout << "Breeders_M: " << breeders_M << std::endl;
            parents_number = breeders_N + breeders_M;
            std::cout << "Parents number: " << parents_number << std::endl;
        }

        void choose1(){
            /* Take the first N as the best from the population and M as lucky randoms from the rest. */
            int temp;
            for (int i = 0; i < breeders_N; i++){ //put the first N individuals to all_parents vector
                all_parents.push_back(sorted_population[i]); 
            }
            auto dist = std::uniform_int_distribution<int>(breeders_N, static_cast<int>(sorted_population.size()));
            for (int i = 0; i < breeders_M; i++){
                temp = dist(rng);
                std::cout << "[Choose1] random index: " << temp << std::endl;
                all_parents.push_back(sorted_population[temp]);
            }
            std::cout << "Unshuffled parents: \n";
            display_population(all_parents, parents_number, sorted_population[1].n);
            std::random_shuffle(std::begin(all_parents), std::end(all_parents)); //suffle all parents
            std::cout << "Shuffled parents:   \n";
            display_population(all_parents, parents_number, sorted_population[1].n);
        }

        void choose2(){
            //pick two random individuals and take better of them;
            //repeat it N+M times
        }

        void choose3(){
            //choose N+M best individuals
        }

};

bool sort_by_fitness(individual a, individual b){
    /* -----
    Returns the third argument for sort() function - for sorting in ascending order
    a vector of individual objects by value of 'fitness' field.
    ----- */
    return a.fitness < b.fitness;
}

class genetic_algorithm {

    //private:
    public:
        int population_size; //I think it shouldn't be here (because it's inside genetic_parameters params), 
                            //but it may have some references in other classes - that has to be checked 
        int n; //number of vertices + 1 - size of graph adjacency matrix - also in params, but I would leave it here for comfort
        std::vector< individual > population; //vector of vectors with vertices' order
        std::vector< individual > new_population;
        std::vector<std::vector<bool> > matrix; //graph
        genetic_parameters params; //parameters, I hope it would be nice to set all parameters in one place
                                    //and then just run solver
        parents chosen_parents; //parents chosen for the next breeding

        std::default_random_engine rng;


        void generate_population(){
            for (int i = 0; i < population_size; i++){
                population.push_back(individual(n)); 
                population[i].initialize_value(); //initialize vector with vertices order by 0
                population[i].generate(); //generate an individual
                population[i].count_fitness(matrix); //run greedy algorithm on the individual
            }
        }




        void sort_population(){
            /* Function sorting population by ascending fitness function. 
            Modifies population vector to make it ordered ascending by population[i].fitness value
            */
            std::sort(std::begin(population), std::end(population), sort_by_fitness);

        }

        void choose_parents(){
            chosen_parents = parents(params.breeders_p, params.breeders_q, population, rng); 
            if (params.parents_choosing == 1){ //I plan to have more than one parents choosing method
                chosen_parents.choose1();
            }
            else {
                chosen_parents.choose2();
            }
            for (int i = 0; i < chosen_parents.all_parents.size(); i++){
                new_population.push_back(chosen_parents.all_parents[i]);
                std::cout << "wpycham element do new_population, i = " << i << std::endl;
            }
        }


        void create_children(){
            //create as many children as rejected individuals
            int size = chosen_parents.all_parents.size();
            int goal = params.population_size - new_population.size();
            std::vector<parents_pair> pairs;
            std::cout << "size = " << size << std::endl;
            std::cout << "size/2 = " << size/2 << std::endl;
            for (int i = 0; i < size/2; i++){
                pairs.push_back(parents_pair(chosen_parents.all_parents[i], chosen_parents.all_parents[size-1-i], n));
                std::cout << "creating pairs, i = " << i << std::endl;
            }
            int i = 0;
            std::cout << "Create children dupa\n" << std::endl;
            while (goal--){
                std::cout << "i = " << i << std::endl;
                pairs[i].create_child1(rng);
                new_population.push_back(pairs[i].child);
                i++;
                if (i == pairs.size()){
                    std::cout << "Reseting i\n" << std::endl;
                    i = 0;
                }
            }
            if (new_population.size() != population_size) {
                std::cout << "Size of new population is incorrect" << std::endl;
            }
            /* 
            Number of children for pair is:
            (population_size - number_of_parents) / (number_of_parents/2)
            where number_of_parents = N + M.
            N + M has to be even. We take pairs from chosen_parents randomly and use 
            chosen_parents.create_child1() or create_child2(), depending on
            create_child_method parameter in parameters object. The result should be added to
            another population vector (I haven't created it yet.) 
            */ 
            

        }

        void mutate(){
            //swap random pairs 
            //don't do it now, we have to think about it
            //anyway, if we decide to mutate, we have to distinguish between methods
        }


    //public:


        genetic_algorithm (std::vector<std::vector<bool> > &input_matrix, genetic_parameters my_params){
            /* Maybe constructor should get an object with parameters?
             Like:
             - population size
             - breeder's number (percent)
                * N (number of best specimen)
                * M (number of random individuals)
             - chance of mutation
             - number saying at what level which 
             parents choosing method should be used
             - mutation method?
             Yes, I added it partially. Thanks. */

            rng = std::default_random_engine();
            n = input_matrix.size();
            matrix = input_matrix;
            params = my_params;
            population_size = my_params.population_size;
            //initialize parents' vector with zeroes
        } 
        ~genetic_algorithm(){};

         
        std::vector <int> solve (){
            while (false){ /*terminating condition is not reached */
                /* Maybe there is another class needed, higher class? 
                We have to:
                1. Generate population
                2. Sort it by fitness.
                3. Choose parents (with proper method)
                4. Create children and add them to new population
                5. Mutate them - not now.
                6. Start algorithm once again with new population and the same parameters.
                 */
            }
        }
        
       

};



void display_parameters (genetic_parameters params){
    std::cout << "Population size: " << params.population_size << std::endl;
    std::cout << "p: " << params.breeders_p << " q: " << params.breeders_q << std::endl;
    std::cout << "Mutation chance: " << params.mutation_chance << "\nParents choosing: " << params.parents_choosing << std::endl;
    std::cout << "Number of generations: " << params.number_of_generations << std::endl;
    std::cout << "Breeding method: " << params.create_child_method << std::endl;
    std::cout << "Mutation size: " << params.mutation_size << std::endl;
}

void display_matrix (std::vector<std::vector<bool> > matrix){
    for (int i = 0; i < matrix.size(); i++){
        for (int j = 0; j < matrix.size(); j++){
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std:: endl;
}

int main(int argc, char const *argv[])
{
    std::vector<std::vector<bool> >matrix(parse_input(argc, argv));

    std::vector<int> result;
    result = greedy(matrix);
    std::cout << "\n----- Greedy algorithm -----";
    print_result(result);
    std::vector<int> result2;
    result2 = largest_first(matrix);
    std::cout << "\n----- Largest first algorithm -----";
    print_result(result2);
    std::vector<int> result3;
    genetic_parameters params = genetic_parameters(10);
    std::cout << "\n----- Genetic algorithm -----\n";

    /* Display input matrix */
    display_matrix(matrix);
    
    /* Display genetic parameters */
    std::cout << "----- Parameters -----\n";
    display_parameters(params);

    /* Display population values, used colours and fitness */
    genetic_algorithm genetic = genetic_algorithm(matrix, params);
    genetic.generate_population();
    std::cout << "\n----- Population -----\n";
    display_population(genetic.population, params.population_size, matrix.size());

    /* Display sorted population */
    genetic.sort_population();
    std::cout << "----- Sorted population -----\n";
    display_population(genetic.population, params.population_size, matrix.size());

    /* Choose and display parents */
    genetic.choose_parents();

    /* Create children */
    genetic.create_children();
    for (int i = 0; i < params.population_size; i++){
        genetic.new_population[i].count_fitness(matrix);
    }
    display_population(genetic.new_population, params.population_size, matrix.size());


    return 0;
}
