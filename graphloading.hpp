#include <iostream>
#include <fstream>
#include <vector>

void print_matrix(std::vector<std::vector<bool>> matrix)
{

    for(unsigned int i = 0; i < matrix.size(); i++)
    {
        for(unsigned int j = 0; j < matrix[i].size(); j++){
            std::cout<<(int)matrix[i][j] << " ";
        }
    std::cout<<std::endl;
}
}



//for data where the lowest vertex index equals 1
std::vector<std::vector<bool> > load_matrix_from_file(std::string filename){
        std::ifstream file;
        file.open(filename);
        int size;
        file >> size;
        //increase the size by one to avoid indexing mess
        ++size;
        std::vector<std::vector<bool> > matrix(size, std::vector<bool>(size));
        int a, b;
        while(file >> a >> b)
        {
            matrix[a][b] = true;
            matrix[b][a] = true;
        }
        file.close();
        return matrix;
}



