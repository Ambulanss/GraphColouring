#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <tuple>
#include <stdexcept>
#include <algorithm>

class AdjacencyMatrix{    
    private:
        std::vector <std::vector <bool> > matrix;
        int sparsity;
        int number_of_edges;

    public:
        
        AdjacencyMatrix(unsigned long long int size, int spars)
        {
            sparsity = spars; //minimum sparsity of graph given in number of vertices
            size = size + 1; //because vertex indexes start from 1
                    
            // init matrix with zeros meaning no vertices
            matrix = std::vector <std::vector <bool> >(size, std::vector<bool>(size)); 

            //fill matrix with edges according to arguments
            std::default_random_engine generator; //random number generating engine
            std::vector<std::tuple<int, int> > unused_edges; //vector of unused edges so far
            number_of_edges = 0; //holder for current number of vertices, debugging purposes mostly
            int roll = 0; //holder for a random value
            //make a connected graph first
            for(unsigned int i = 1; i < matrix.size(); i++)
            {
                //n operations like this
                if(i != 1)
                {
                    std::uniform_int_distribution<int> distribution(1, i-1);
                    roll = distribution(generator);
                    matrix[i][roll] = true;
                    matrix[roll][i] = true;
                    number_of_edges++;
                }
            }

            for(int x = 1; x < matrix.size(); x++)
            {
                
                for(int y = 1; y < x; y++)
                {
                    //n * (n-1)/2 operations like this
                    if(matrix[x][y] == false)
                    {
                        auto temp = std::make_tuple(x, y);
                        unused_edges.push_back(temp);
                    }
                }
            }

            int maximum_number_of_edges = ((size-1) * (size - 2))/2;// n*(n-1)/2 because graph is undirected
            float edge_goal = static_cast<float>(sparsity)/100.0 * static_cast<float>(maximum_number_of_edges);
            std::cout<<std::endl<<edge_goal<<std::endl;
            
            std::shuffle(unused_edges.begin(), unused_edges.end(), generator);//linear complexity
            while(number_of_edges < edge_goal && !unused_edges.empty())
            {
                auto tmp = unused_edges.back();
                matrix[std::get<0>(tmp)][std::get<1>(tmp)] = true;
                matrix[std::get<1>(tmp)][std::get<0>(tmp)] = true;
                number_of_edges++;
                unused_edges.pop_back(); 

            }
            std::cout<<"Number of edges achieved: "<< number_of_edges<<std::endl;
        }
        
        void print_matrix()
        {
            //prints the hole matrix starting from index number 1
            for(unsigned int i = 1; i < matrix.size(); i++)
            {
                for(unsigned int j = 1; j < matrix[i].size(); j++){
                    std::cout<<(int)matrix[i][j] << " ";
                }
            std::cout<<std::endl;
            }
        }
        
        void write_to_file(std::string filename)
        {
            std::ofstream stream(filename, std::ofstream::out);
            stream << matrix.size() - 1<<std::endl;
            for(unsigned int i = 1; i < matrix.size(); i++)
            {
                for(unsigned int j = 1; j < i; j++)
                {
                    if(matrix[i][j] == true)
                    {
                        stream<<j<<" "<<i<<"\n";
                    }
                }
            }
               
        }

        std::vector<bool>& operator[](int x) //overloading [] operator to be able to access matrix data from outside
        {
            return matrix[x];
        }

};

void parse_and_check_input(int argc, char const *argv[], long long int & size, int &sparsity, std::string & filename)
{
    if(argc != 4){
        throw "Usage:\ngenerator.exe <number_of_edges> <sparsity_as_full_percents> <filename>\n";
    }
    size = atol(argv[1]);
    std::cout<<"Size: "<< size<<std::endl;
    sparsity = atoi(argv[2]);
    std::cout<<"Sparsity: "<< sparsity <<std::endl;
    filename = argv[3];

    if(size < 1 || size > 25000)
    {
        throw "Wrong size: must be from 1 to 25000\n";
    }
    if(sparsity < 1 || sparsity > 100)
    {
        throw "Sparsity must be given as an integer ranging from 1 to 100.\n";
    }
    std::cout<<"Input correct"<<std::endl;
}    //std::vector < std::vector<bool> >test;
    //auto max_sz = static_cast<long long>(test.max_size());

int main(int argc, char const *argv[])
{
    long long int size;
    int sparsity;
    std::string filename;
    
    try
    {
        parse_and_check_input(argc, argv, size, sparsity, filename);
    }
    catch(const char* e)
    {
        std::cout << e;
        exit(-1);
    }
    
    try
    {
        AdjacencyMatrix graph(size, sparsity);
        if(size <= 50)
            graph.print_matrix();
        graph.write_to_file(filename);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        exit(-2);
    }
    
    

    
    return 0;
}
