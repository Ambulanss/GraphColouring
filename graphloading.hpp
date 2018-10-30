#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

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
auto load_matrix_from_txt(std::string filename){
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

auto load_matrix_from_col(std::string filename)
{
    std::ifstream file(filename);
    char line_marker;   //storage for the first character of a line 
    int a, b;           //storage for vertex ids
    std::string line;   //storage for a line of the source file
    std::vector<std::vector<bool> > matrix;
    std::string graph_type;
    while(std::getline(file, line))
    {
        std::istringstream iss(line);
        iss >> line_marker;
        if(line_marker == 'c'){
            continue;
        }
        else
        if(line_marker == 'p')
        {
            iss >> graph_type >> a >> b;
            int size = a + 1;
            matrix.resize(size, std::vector<bool>(size));
        }
        else
        if(line_marker == 'e'){
            iss >> a >> b;
            matrix[a][b] = true;
            matrix[b][a] = true;
        }
        else
        {
            std::cout<<"\nCorrupted graph source file:" << filename<<"\n Line: "<< line <<std::endl;
        }
    }
    std::cout<< "Check if matrix is ok: \n";
    print_matrix(matrix);
    return matrix;
}

std::string getFileExt(const std::string& s) {

   size_t i = s.rfind('.', s.length());
   if (i != std::string::npos) {
      return(s.substr(i+1, s.length() - i));
   }

   return("");
}

auto load_matrix(std::string filename)
{
 /*   if(filetype == 0)
    {
        return load_matrix_from_txt(filename);
    }
    else
    {
        return load_matrix_from_col(filename);
    }*/
    if(getFileExt(filename).compare("col") == 0){
        return load_matrix_from_col(filename);
    }
    else
    if(getFileExt(filename).compare("txt") == 0)
    {
        return load_matrix_from_txt(filename);
    }
    else
    {
        std::cout<<"Bledny format pliku"<<std::endl;
        exit(-2);
    }
}
