/**
 * File:   gmsh_io.hpp
 * Author: akirby
 *
 * Created on November 24, 2023, 3:13 PM
 */

#ifndef GMSH_IO_HPP
#define GMSH_IO_HPP

/* system header files */
# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <fstream>
# include <iomanip>
# include <iostream>
# include <vector>

#ifdef __cplusplus
extern "C" {
#endif

// Element types
#define TYPE_PNT     1
#define TYPE_LIN     2
#define TYPE_TRI     3
#define TYPE_QUA     4
#define TYPE_TET     5
#define TYPE_PYR     6
#define TYPE_PRI     7
#define TYPE_HEX     8
#define TYPE_NUM_MAX 9

#define TYPE_NAME (char[][TYPE_NUM_MAX]){"",        \
                                         "POINT  ", \
                                         "LINE   ", \
                                         "TRI    ", \
                                         "QUAD   ", \
                                         "TET    ", \
                                         "PYRAMID", \
                                         "PRISM  ", \
                                         "HEX    "  }

// Globals
std::vector<double> node_x;
std::vector<int> element_node[TYPE_NUM_MAX];
std::vector<int> wbc;
std::vector<int> obc;

#define TYPE_NNODES (int[9]){0,1,2,3,4,4,5,6,8}
static
void elementTypeName(int type){
    switch(type){
        case 1:  cout << "[1D] ElemType: 2-node line.\n"; break;
        case 2:  cout << "[2D] ElemType: 3-node triangle. \n"; break;
        case 3:  cout << "[2D] ElemType: 4-node quadrangle. \n"; break;
        case 4:  cout << "[3D] ElemType: 4-node tetrahedron. \n"; break;
        case 5:  cout << "[3D] ElemType: 8-node hexahedron. \n"; break;
        case 6:  cout << "[3D] ElemType: 6-node prism. \n"; break;
        case 7:  cout << "[3D] ElemType: 5-node pyramid. \n"; break;
        case 8:  cout << "[1D] ElemType: 3-node second order line (2 nodes associated with the vertices and 1 with the edge). \n"; break;
        case 9:  cout << "[2D] ElemType: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges). \n"; break;
        case 10: cout << "[2D] ElemType: 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face). \n"; break;
        case 11: cout << "[3D] ElemType: 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges). \n"; break;
        case 12: cout << "[3D] ElemType: 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume). \n"; break;
        case 13: cout << "[3D] ElemType: 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces). \n"; break;
        case 14: cout << "[3D] ElemType: 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face). \n"; break;
        case 15: cout << "[1D] ElemType: 1-node point. \n"; break;
        case 16: cout << "[2D] ElemType: 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges). \n"; break;
        case 17: cout << "[3D] ElemType: 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges). \n"; break;
        case 18: cout << "[3D] ElemType: 15-node second order prism (6 nodes associated with the vertices and 9 with the edges). \n"; break;
        case 19: cout << "[3D] ElemType: 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges). \n"; break;
        case 20: cout << "[2D] ElemType: 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges) \n"; break;
        case 21: cout << "[2D] ElemType: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face) \n"; break;
        case 22: cout << "[2D] ElemType: 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges) \n"; break;
        case 23: cout << "[2D] ElemType: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face) \n"; break;
        case 24: cout << "[2D] ElemType: 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges) \n"; break;
        case 25: cout << "[2D] ElemType: 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face) \n"; break;
        case 26: cout << "[1D] ElemType: 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge) \n"; break;
        case 27: cout << "[1D] ElemType: 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge) \n"; break;
        case 28: cout << "[1D] ElemType: 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge) \n"; break;
        case 29: cout << "[3D] ElemType: 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces) \n"; break;
        case 30: cout << "[3D] ElemType: 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume) \n"; break;
        case 31: cout << "[3D] ElemType: 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume) \n"; break;
        case 92: cout << "[3D] ElemType: 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume) \n"; break;
        case 93: cout << "[3D] ElemType: 125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume) \n"; break;
        default: cout << "ELEMENT TYPE NOT RECOGNIZED! " << type << endl; break;
    }
}

static
int getElementType(int type){
    switch(type){
        case 1:  return TYPE_LIN; //cout << "[1D] ElemType: 2-node line.\n"; break;
        case 2:  return TYPE_TRI; //cout << "[2D] ElemType: 3-node triangle. \n"; break;
        case 3:  return TYPE_QUA; //cout << "[2D] ElemType: 4-node quadrangle. \n"; break;
        case 4:  return TYPE_TET; //cout << "[3D] ElemType: 4-node tetrahedron. \n"; break;
        case 5:  return TYPE_HEX; //cout << "[3D] ElemType: 8-node hexahedron. \n"; break;
        case 6:  return TYPE_PRI; //cout << "[3D] ElemType: 6-node prism. \n"; break;
        case 7:  return TYPE_PYR; //cout << "[3D] ElemType: 5-node pyramid. \n"; break;
        case 8:  return TYPE_LIN; //cout << "[1D] ElemType: 3-node second order line (2 nodes associated with the vertices and 1 with the edge). \n"; break;
        case 9:  return TYPE_TRI; //cout << "[2D] ElemType: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges). \n"; break;
        case 10: return TYPE_QUA; //cout << "[2D] ElemType: 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face). \n"; break;
        case 11: return TYPE_TET; //cout << "[3D] ElemType: 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges). \n"; break;
        case 12: return TYPE_HEX; //cout << "[3D] ElemType: 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume). \n"; break;
        case 13: return TYPE_PRI; //cout << "[3D] ElemType: 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces). \n"; break;
        case 14: return TYPE_PYR; //cout << "[3D] ElemType: 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face). \n"; break;
        case 15: return TYPE_PNT; //cout << "[1D] ElemType: 1-node point. \n"; break;
        case 16: return TYPE_QUA; //cout << "[2D] ElemType: 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges). \n"; break;
        case 17: return TYPE_HEX; //cout << "[3D] ElemType: 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges). \n"; break;
        case 18: return TYPE_PRI; //cout << "[3D] ElemType: 15-node second order prism (6 nodes associated with the vertices and 9 with the edges). \n"; break;
        case 19: return TYPE_PYR; //cout << "[3D] ElemType: 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges). \n"; break;
        case 20: return TYPE_TRI; //cout << "[2D] ElemType: 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges) \n"; break;
        case 21: return TYPE_TRI; //cout << "[2D] ElemType: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face) \n"; break;
        case 22: return TYPE_TRI; //cout << "[2D] ElemType: 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges) \n"; break;
        case 23: return TYPE_TRI; //cout << "[2D] ElemType: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face) \n"; break;
        case 24: return TYPE_TRI; //cout << "[2D] ElemType: 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges) \n"; break;
        case 25: return TYPE_TRI; //cout << "[2D] ElemType: 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face) \n"; break;
        case 26: return TYPE_LIN; //cout << "[1D] ElemType: 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge) \n"; break;
        case 27: return TYPE_LIN; //cout << "[1D] ElemType: 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge) \n"; break;
        case 28: return TYPE_LIN; //cout << "[1D] ElemType: 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge) \n"; break;
        case 29: return TYPE_TET; //cout << "[3D] ElemType: 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces) \n"; break;
        case 30: return TYPE_TET; //cout << "[3D] ElemType: 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume) \n"; break;
        case 31: return TYPE_TET; //cout << "[3D] ElemType: 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume) \n"; break;
        case 92: return TYPE_HEX; //cout << "[3D] ElemType: 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume) \n"; break;
        case 93: return TYPE_HEX; //cout << "[3D] ElemType: 125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume) \n"; break;
        default: cout << "ELEMENT TYPE NOT RECOGNIZED! " << type << endl; return 0; break;
    }
}

static
int getElementNumNodes(int type){
    switch(type){
        case 1:  return 2;  //cout << "[1D] ElemType: 2-node line.\n"; break;
        case 2:  return 3;  //cout << "[2D] ElemType: 3-node triangle. \n"; break;
        case 3:  return 4;  //cout << "[2D] ElemType: 4-node quadrangle. \n"; break;
        case 4:  return 4;  //cout << "[3D] ElemType: 4-node tetrahedron. \n"; break;
        case 5:  return 8;  //cout << "[3D] ElemType: 8-node hexahedron. \n"; break;
        case 6:  return 6;  //cout << "[3D] ElemType: 6-node prism. \n"; break;
        case 7:  return 5;  //cout << "[3D] ElemType: 5-node pyramid. \n"; break;
        case 8:  return 3;  //cout << "[1D] ElemType: 3-node second order line (2 nodes associated with the vertices and 1 with the edge). \n"; break;
        case 9:  return 6;  //cout << "[2D] ElemType: 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges). \n"; break;
        case 10: return 9;  //cout << "[2D] ElemType: 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face). \n"; break;
        case 11: return 10; //cout << "[3D] ElemType: 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges). \n"; break;
        case 12: return 27; //cout << "[3D] ElemType: 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume). \n"; break;
        case 13: return 18; //cout << "[3D] ElemType: 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces). \n"; break;
        case 14: return 14; //cout << "[3D] ElemType: 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face). \n"; break;
        case 15: return 1;  //cout << "[1D] ElemType: 1-node point. \n"; break;
        case 16: return 8;  //cout << "[2D] ElemType: 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges). \n"; break;
        case 17: return 20; //cout << "[3D] ElemType: 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges). \n"; break;
        case 18: return 15; //cout << "[3D] ElemType: 15-node second order prism (6 nodes associated with the vertices and 9 with the edges). \n"; break;
        case 19: return 13; //cout << "[3D] ElemType: 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges). \n"; break;
        case 20: return 9;  //cout << "[2D] ElemType: 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges) \n"; break;
        case 21: return 10; //cout << "[2D] ElemType: 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face) \n"; break;
        case 22: return 12; //cout << "[2D] ElemType: 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges) \n"; break;
        case 23: return 15; //cout << "[2D] ElemType: 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face) \n"; break;
        case 24: return 15; //cout << "[2D] ElemType: 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges) \n"; break;
        case 25: return 21; //cout << "[2D] ElemType: 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face) \n"; break;
        case 26: return 4;  //cout << "[1D] ElemType: 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge) \n"; break;
        case 27: return 5;  //cout << "[1D] ElemType: 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge) \n"; break;
        case 28: return 6;  //cout << "[1D] ElemType: 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge) \n"; break;
        case 29: return 20; //cout << "[3D] ElemType: 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces) \n"; break;
        case 30: return 35; //cout << "[3D] ElemType: 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume) \n"; break;
        case 31: return 56; //cout << "[3D] ElemType: 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume) \n"; break;
        case 92: return 64; //cout << "[3D] ElemType: 64-node third order hexahedron (8 nodes associated with the vertices, 24 with the edges, 24 with the faces, 8 in the volume) \n"; break;
        case 93: return 125; //cout << "[3D] ElemType: 125-node fourth order hexahedron (8 nodes associated with the vertices, 36 with the edges, 54 with the faces, 27 in the volume) \n"; break;
        default: cout << "ELEMENT TYPE NOT RECOGNIZED! " << type << endl; return 0; break;
    }
}

char ch_cap(char ch);

bool ch_eqi(char ch1,char ch2);

int ch_to_digit(char ch);

void gmsh_size_read(string gmsh_filename,int &node_num,int &node_dim,
                    int &element_num,int elem_counts[]);

void gmsh_data_read(string gmsh_filename,int node_dim,
                    std::vector<double> &node_x,
                    std::vector<int> element_nodes[]);

bool s_begin(string s1, string s2 );

int s_len_trim(string s );

int s_to_i4(string s, int &last, bool &error );

double s_to_r8(string s, int &lchar, bool &error );

void timestamp();

#ifdef __cplusplus
}
#endif
#endif /* GMSH_IO_HPP */