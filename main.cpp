#include <iostream>
#include <vector>

int lowest_colour(std::vector<bool> &row, std::vector<int> &result){
    /* -----
    Returns the lowest colour for colouring given vertex.

    @row - row from adjacency matrix, containing neighbours of vertex that is being processed
    @result - vector containing current number of used colours (at index 0)
            and colours of vertices (index means a vertex number and value is colour)
    ----- */

    std::vector<bool> possible (result[0]+2, true); //vector of all colours used so far + 1
    for (int i = 1; i < row.size(); i++){
        if (row[i]){ //if a neighbour is found
            if (result[i]){ //if a neighbour is coloured - his value in result vector is != 0
                possible[result[i]] = false; //set this colour as impossible to use
            }
        }
    }
    for (int i = 1; i < possible.size(); i++){
        if (possible[i]){ //if a colour is possible to use
            if (i == possible.size()-1){ //if chosen colour is new
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
    for (int i = 1; i < matrix.size(); i++){
        colour = lowest_colour(matrix[i], result);
        result[i] = colour;
    }
    return result;
}

int main(){
    
    return 0;
}
