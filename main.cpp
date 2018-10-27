#include <iostream>
#include "graphloading.hpp"

int main(int argc, char const *argv[]){
    if (argc > 1)
    {
        auto test_mat(load_matrix_from_file(static_cast <std::string>(argv[1])));
    }
    else
    {
        std::cout<<"Podaj nazwe pliku zrodlowego."<<std::endl;
    }
    return 0;
}
