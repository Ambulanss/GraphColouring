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
        int breeders_N; //percent of all individuals; parents chosen as the best of all - TODO: it has to be double
        int breeders_M; //percent of all individuals; parents chosen randomly - TODO: it has to be double
        int mutation_change; //chance for mutation, expressed as percent - TODO: it has to be double
        /* parents_choosing - number expressing when should we use 
        first parents choosing method and when the second one. 
        It means after which generation should we switch to the second method. */
        int parents_choosing; 
        int number_of_generations; //number of generations before stopping the algorithm
        int create_child_method; //if 1, use first mutation method, if 2, use the second one
        double mutation_size; //percent of values in child that should be mutated; expressed as a fraction
};

class parents_pair{
    /* one pair of parents */
    public:
        individual parent_A;
        individual parent_B;
        individual child;
        int n;

        parents_pair(individual A, individual B, int input_n){ //it cannot be like this;
        /* individual needs two parameters for constructor. Random engine cannot be one of them (it has to be changed),
        because it makes a mess. I would leave only n. Because of that, we need to pass n as a constructor
        parameter for parents_pair class I guess, and then set parent_A and parent_B parameters in function
        like "set_parents". */ 
            n = input_n;
            parent_A = A;
            parent_B = B;
        }

        void create_child1(){
            /* We take first values from parent_A (until crosspoint index) and the rest
            from parent_B (the ones that are possible to use without repeating). 
            We then will the gaps with random values chosen from available ones.

            I would add functionality to randomly use parent_A's values at the beginning 
            of the child (before crosspoint) or at the end (after crosspoint).
            */
            child.initialize_value(); //now the child vector is filled with zeroes
            std::vector <bool> available_values (child.n, true); //should create vector with n true values
            available_values[0] = false; //we don't have vertex nr 0
            int crosspoint;
            std::vector <int> left_indexes; //we will  use them to know where the gaps are
            bool fill_the_gaps = false; //if this is true, we know that there are some gaps (it will fore sure always be true)
            /* TODO: crosspoint = random value from [1, child.n) */
            for (int i = 1; i < crosspoint; i++){ //we take first values from parent A
                child.value[i] = parent_A.value[i];
                available_values[parent_A.value[i]] = false; //we mark used vertex as impossible to use
                /* reminder: value field contains vertices' numbers, so parent_A.value[2] = 3 means
                that parent_A is vertex order in which the second vertex processed by greedy algorithm
                is the third vertex in the input graph */
            }
            for (int i = crosspoint; i < child.n; i++){ //the rest is filled with possible values from parent_B
                if (available_values[parent_B.value[i]] == true){ //if vertex hasn't been used yet
                    child.value[i] == parent_B.value[i]; //then insert it into child
                }
                else { //remember that you have to fill this gap later
                    left_indexes.push_back(i); //index i has no vertex
                    fill_the_gaps = true;
                } 
            }
            if (fill_the_gaps){ //if there are some gaps to fill
                /* TODO: shuffle left_indexes vertex... oh my, I think our random engine is needed as a parameter 
                to create_child function*/
                for (int i = 1; i < child.n; i++){
                    if(available_values[i]){ //if unused vertex is found
                        child.value[left_indexes[0]] = available_values[i]; 
                        /* assign first available vertex to random gap in child's values - random, because 
                        free left indexes were shuffled. */
                        available_values[i] = false; //set this vertex as used
                        /* TODO: remove first element (on index 0) from left_indexes vector;
                        it can't be used anymore. */
                    }
                }
            }
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

        void mutate_child1(double mutation_size){
            /* Function swaps two random values in child. 
            We don't use it yet. */
            int mutation_number = floor(((child.n-1) * mutation_size)); 
            for (int i = 0; i < mutation_number; i++){
            /* TODO: Generate two random indexes from [1, child.n);
            Assign them to random1 and random2. */
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
        std::vector <individual> sorted_population; //population sorted by fitness, descending
        std::vector <individual> all_parents;
        int breeders_N; //number of parents to choose as the best individuals
        int breeders_M; //number of parents to choose as random individuals
        std::default_random_engine rng;

        parents(int N, int M, std::vector <individual> input_population, std::default_random_engine my_rng){
            /* Constructor initializing object keeping all parents in population.
            Parameters:
            @N - percent of parents to choose as the best individuals --- it has to be double ;_; 
                                                                        //it's treated as fraction later
            @M - percent of parents to choose as random individuals --- just like above
            @input_population - sorted population vector */
            rng = my_rng;
            sorted_population = input_population;
            if (N + M > 100){
                std::cout << "Error: wrong number of breeders. \n";
            }
            breeders_N = N*sorted_population.size();
            breeders_M = M*sorted_population.size();
        }

        void choose1(){
            /* Take the first N as the best from the population and M as lucky randoms from the rest. */
            int temp;
            for (int i = 0; i < breeders_N; i++){ //put the first N individuals to all_parents vector
                all_parents.push_back(sorted_population[i]); 
            }
            for (int i = 0; i < breeders_M; i++){
                /* TODO: temp = random value from [N, sorted_population.size()-1] */
                all_parents.push_back(sorted_population[temp]);
            }
            std::shuffle(std::begin(all_parents), std::end(all_parents), rng); //suffle all parents
        }

        void choose2(){
            //pick two random individuals and take better of them;
            //repeat it N+M times
        }

        void choose3(){
            //choose N+M best individuals
        }

};

class individual {
        /* value vector keeps the order of vertices to process by typical greedy algorithm.
        We iterate not by i = 1, 2, 3, ..., n-1, but by values of value vector.
        For example, value[1] = 3 means that the first vertex we should colour is vertex nr 3; 
        value[2] = 7 means that the second vertex we should colour is vertex nr 7; and so on... */

    public: 
        std::vector <int> value; //vertices
        int fitness; 
        int n; //size of main adjacency matrix
        std::default_random_engine rng;
        std::vector <int> result_colours; //result_colours[0] is fitness xD - after running count_fitness

        individual (std::default_random_engine my_rng, int my_n){ //engine cannot be passed as a constructor parameter;
                                                                //it has to be only shuffle method argument
            rng = my_rng;
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

        void generate() { //generate random permutation - random vertices' order
            for (int i = 1; i < n; i++){ //first assign value equal to index to later shuffle it
                value[i] = i;
            }
            std::shuffle(std::begin(value)+1, std::end(value), rng); 
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

class genetic_algorithm {
    /* Here goes our super class, superhero or god. */
    private:
        int population_size; //I think it shouldn't be here (because it's inside genetic_parameters params), 
                            //but it may have some references in other classes - that has to be checked 
        int n; //number of vertices + 1 - size of graph adjacency matrix - also in params, but I would leave it here for comfort
        std::vector< individual > population; //vector of vectors with vertices' order
        std::vector<std::vector<bool> > matrix; //graph
        genetic_parameters params; //parameters, I hope it would be nice to set all parameters in one place
                                    //and then just run solver
        parents chosen_parents; //parents chosen for the next breeding

        //here started my neurosis
        std::default_random_engine rng;


        void generate_population(){
            for (int i = 0; i < population_size; i++){
                population.push_back(individual(rng, n)); 
                population[i].initialize_value(); //initialize vector with vertices order by 0
                population[i].generate(); //generate an individual
                population[i].count_fitness(matrix); //run greedy algorithm on the individual
            }
        }


        bool sort_by_fitness(individual a, individual b){
            /* -----
            Returns the third argument for sort() function - for sorting in descending order
            a vector of individual objects by value of 'fitness' field.
            ----- */
            return a.fitness > b.fitness;
}

        void sort_population(){
            /* Function sorting population by descending fitness function. 
            Modifies population vector to make it ordered descending by population[i].fitness value
            */
            std::sort(std::begin(population), std::end(population), sort_by_fitness);

        }

        void choose_parents(){
            parents chosen_parents(params.breeders_N, params.breeders_M, population, rng); //fucking constructor with many params
            if (params.parents_choosing == 1){ //I plan to have more than one parents choosing method
                chosen_parents.choose1();
            }
            else {
                chosen_parents.choose2();
            }
        }


        void create_children(){
            //create as much children as rejected individuals
            /* I don't have power to write it today. 
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


    public:


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
            
            std::random_device r;
            rng = rng(r()); //arghhh
            n = input_matrix.size();
            matrix = input_matrix;
            params = my_params;
            population_size = my_params.population_size;
            //initialize parents' vector with zeroes
        } 
        ~genetic_algorithm(){};

         
        std::vector <int> solve (){
            while (/*terminating condition is not reached */){
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
    return 0;
}
