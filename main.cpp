#include <iostream>
#include <vector>
#include <bits/stdc++.h> 
#include "graphloading.hpp"

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
    //int filetype;
    if(argc > 1){
        filename = static_cast<std::string>(argv[1]);
        /*if(argc > 2)
        {
            filetype = atoi(argv[2]);
            matrix = load_matrix(filename, filetype);
        }
        else
        {
            matrix = load_matrix(filename, 0);
        }*/
        matrix = load_matrix(filename);
        return matrix; 
    }
    else
    {
        std::cout<<"Podaj przynajmniej nazwe pliku"<<std::endl;
        exit(-1);
    }
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
    return 0;
}
//TODO
/*
1. Generator instancji
2. 
*/