/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* system header files */
# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>
#include <vector>

/* name space */
using namespace std;

/* header files */
# include "gmsh_io.hpp"

void test();

int main(){
    timestamp();
    cout << "\n";
    cout << "GMSH_IO_TEST\n";
    cout << "  C++ version\n";
    cout << "  Test the GMSH_IO library.\n";

    test();

    cout << "\n";
    cout << "GMSH_IO_TEST\n";
    cout << "  Normal end of execution.\n";
    cout << "\n";
    timestamp();
    return 0;
}

void test(){
    string gmsh_filename = "billet-cap.3D.Q1.12K.msh";
    int dim,node_num;
    int element_num;

    cout << "\nTEST:\n Read data from a file.\n";

    //  Get the data size.
    std::vector<int> elem_counts(TYPE_NUM_MAX,0);
    gmsh_size_read(gmsh_filename,node_num,dim,element_num,elem_counts.data());

    cout << " \n Node data read from file \"" << gmsh_filename << "\"\n\n";
    cout << "  Number of nodes = "      << node_num << "\n";
    cout << "  Spatial dimension = "    << dim << "\n";
    cout << "  Number of elements: "    << "\n";
    for(int i = 1; i < TYPE_NUM_MAX; i++){
        printf("    %s %d\n", TYPE_NAME[i],elem_counts[i]);
    }
    cout << endl;

    //  Allocate memory.
    std::vector<double> node_x(dim*node_num);
    std::vector<std::vector<int>> element_node;

    element_node.resize(TYPE_NUM_MAX);
    for(int i = 1; i < TYPE_NUM_MAX; i++){
        element_node[i].resize(elem_counts[i]*TYPE_NNODES[i]);
      //printf("ALLOCATING %d %d %d\n",i,elem_counts[i],TYPE_NNODES[i]);
    }

    //  Get the data.
    gmsh_data_read(gmsh_filename,dim,
                   node_x.data(),
                   element_node);

//    for(int i = 1; i < TYPE_NUM_MAX; i++){
//        printf("ELEMENT TYPE %s\n",TYPE_NAME[i]);
//        for(int j = 0; j < elem_counts[i]; j++){
//            for(int k = 0; k < TYPE_NNODES[i]; k++){
//                cout << " " << element_node[i][TYPE_NNODES[i]*j + k];
//            }
//            cout << endl;
//        }
//    }
}