#include <iostream>
#include <fstream>
#include <vector>
#include <random>

class AdjacencyMatrix{    
    private:
        std::vector <std::vector <bool> > matrix;
        int sparsity;
        int number_of_edges;

    public:
        
        AdjacencyMatrix(int size, int spars)
        {
            sparsity = spars; //minimum sparsity of graph given in number of vertices
            
            try
            {
                // init matrix with zeros meaning no vertices
                matrix = std::vector <std::vector <bool> >(size, std::vector<bool>(size)); 
            }
            catch(const std::exception& e)
            {
                std::cerr << e.what() << '\n';
            }
            
            
            std::default_random_engine generator;
            number_of_edges = 0; //holder for current number of vertices
            int roll = 0; //holder for a random value
            for(int i = 1; i < size; i++)
            {
                if(i != 1)
                {
                    std::uniform_int_distribution<int> distribution(1, i-1);
                    roll = distribution(generator);
                    matrix[i][roll] = true;
                    matrix[roll][i] = true;
                    number_of_edges++;
                }
            }

            int maximum_number_of_edges = (size * (size -1))/2;
            float edge_goal = static_cast<float>(sparsity)/100.0 * static_cast<float>(maximum_number_of_edges);
            std::cout<<std::endl<<edge_goal<<std::endl;
            
            while(number_of_edges < edge_goal)
            {
                std::uniform_int_distribution<int> distribution(1, size - 1);
                int roll_x = distribution(generator);
                int roll_y = distribution(generator);
                if(matrix[roll_x][roll_y] != true)
                {
                    matrix[roll_x][roll_y] = true;
                    matrix[roll_y][roll_x] = true;
                    number_of_edges++;
                }
            }

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
            std::ofstream file(filename, std::ofstream::out);
            std::string output = "";
            file << matrix.size() - 1<<std::endl;
            for(unsigned int i = 1; i < matrix.size(); i++)
            {
                for(unsigned int j = 1; j < i; j++)
                {
                    if(matrix[i][j] == true)
                    {
                        file<<j<<" "<<i<<std::endl;
                    }
                }
            }
            file.close();
            
        }
        
        std::vector<bool>& operator[](int x) //overloading [] operator to be able to access matrix data from outside
        {
            return matrix[x];
        }

};

void parse_and_check_input(int argc, char const *argv[], int & size, int &sparsity, std::string & filename)
{
    if(argc != 4){
        std::cout<<"Usage:\ngenerator.exe <number_of_edges> <sparsity_as_full_percents> <filename>\n";
        exit(-1);
    }
    size = atoi(argv[1]) + 1;
    sparsity = atoi(argv[2]);
    filename = argv[3];
    if(size < 1)
    {
        std::cout<<"Size must be bigger than 0\n";
        exit(-2);
    }
    if(sparsity < 1 || sparsity > 100)
    {
        std::cout<<"Sparsity must be given as an integer ranging from 1 to 100.\n";
        exit (-2);
    }
}

int main(int argc, char const *argv[])
{
    int size, sparsity;
    std::string filename;
    parse_and_check_input(argc, argv, size, sparsity, filename);
    try
    {
        AdjacencyMatrix graph(size, sparsity);
        graph.print_matrix();
        graph.write_to_file(filename);
    }
    catch(const std::exception& e)
    {
        std::cerr << e.what() << '\n';
    }
    
    return 0;
}
