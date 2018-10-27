#include <iostream>
#include <fstream>
#include <vector>

//for data where the lowest vertex index equals 1
void load_file( std::string filename, std::vector<std::vector<bool>> &matrix){
    std::ifstream file;
    file.open(filename);
    int size;
    file >> size;
    //increase the size by one to avoid indexing mess
    size += 1;
    matrix.resize(size);
    //from 1 to size - 1

    for(int i = 1; i < size; i++)
    {
        matrix[i].resize(size, false);
    }
    int a, b;
    while(file >> a >> b)
    {
        matrix[a][b] = true;
        matrix[b][a] = true;
    }

    file.close();

}