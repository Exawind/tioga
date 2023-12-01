/**
 * \file    gmsh_io.C
 */

/* name space */
using namespace std;

/* header files */
# include "gmsh_io.hpp"

#define BASE 1
#define _X 0
#define _Y 1
#define _Z 2

/* ========================================================================== */
/*                               GMSH INTERFACE                               */
/* ========================================================================== */
extern "C" {
  void tioga_gmsh_reader_(int *nnodes,double **xyz,
                          int *nwbc, int *nobc,
                          int **wbc_t,int **obc_t,
                          int *n4,int *n5,int *n6,int *n8,
                          int **ndc4,int **ndc5,int **ndc6,int **ndc8){
    std::string gmsh_filename = "billet-cap-big.3D.Q1.60K.tet.msh";
    int dim,node_num;
    int element_num;
    double eps=1.0E-12;
    double x,y,z;

    // remove excess memory
    obc.clear();
    wbc.clear();
    node_x.clear();
    for(int i = 0; i < TYPE_NUM_MAX; i++) element_node[i].clear();

    cout << "\nGMSH:\n Read data from a file.\n";

    //  Get the data size.
    std::vector<int> elem_counts(TYPE_NUM_MAX,0);
    gmsh_size_read(gmsh_filename,node_num,dim,element_num,elem_counts.data());

    cout << " Node data read from file \"" << gmsh_filename << "\"\n\n";
    cout << "  Number of nodes = "      << node_num << "\n";
    cout << "  Spatial dimension = "    << dim << "\n";
    cout << "  Elements counts: "    << "\n";
    for(int i = 1; i < TYPE_NUM_MAX; i++) printf("    %s %d\n",TYPE_NAME[i],elem_counts[i]);
    cout << endl;

    //  allocate memory
    node_x.resize(dim*node_num);
    for(int i = 1; i < TYPE_NUM_MAX; i++){
        element_node[i].resize(elem_counts[i]*TYPE_NNODES[i]);
      //printf("ALLOCATING %d %d %d\n",i,elem_counts[i],TYPE_NNODES[i]);
    }

    // read gmsh mesh data
    gmsh_data_read(gmsh_filename,dim,node_x,element_node);

    // fixed geometry boundaries for msh: DO NOT CHANGE
    double outerboxTop = 1.60; //y-direction
    double innerboxTop = 1.30; //y-direction
    double outerboxh = 0.6;
    double ibh = 0.3;
    double rad = 1.0;

    // count overset bc nodes
    std::vector<char> bcall(node_num,0);
    for(int i=0;i<node_num;i++) bcall[i] =  (node_x[3*i+_Y] >= outerboxTop-eps);// top overset face
    for(int i=0;i<node_num;i++) bcall[i] += (node_x[3*i+_X] <= -outerboxh+eps); // xlo overset face
    for(int i=0;i<node_num;i++) bcall[i] += (node_x[3*i+_X] >=  outerboxh-eps); // xhi overset face
    for(int i=0;i<node_num;i++) bcall[i] += (node_x[3*i+_Z] <= -outerboxh+eps); // zlo overset face
    for(int i=0;i<node_num;i++) bcall[i] += (node_x[3*i+_Z] >=  outerboxh-eps); // zhi overset face

    // fill overset bc nodes
    int nobc_ = 0;
    for(int i=0;i<node_num;i++) nobc_ += (bcall[i] > 0);

    obc.resize(nobc_); nobc_ = 0;
    for(int i=0;i<node_num;i++) {
      if(bcall[i] > 0) obc[nobc_++] = i+BASE;
    }

    // count wall bc nodes
    std::fill(bcall.begin(), bcall.end(), 0);
    for(int i=0;i<node_num;i++) {
        x = node_x[3*i+_X]; y = node_x[3*i+_Y]; z = node_x[3*i+_Z];
        double radiusSqr = x*x + y*y + z*z;

         // sphere surface
        bcall[i] = (radiusSqr <= rad+eps);
//        if(x >= -1.5*ibh-eps && x <= 1.5*ibh+eps){
//            if(z >= -1.5*ibh-eps && z <= 1.5*ibh+eps){
//                bcall[i] = (radiusSqr <= rad+eps);
//            }
//        }

        if(y <= innerboxTop+eps){
            if(x >= -ibh-eps && x <= ibh+eps){
                if(z >= -ibh-eps && z <= ibh+eps){
                    bcall[i] = 1;
                }
            }
        }
    }

    // fill wall bc nodes
    int nwbc_ = 0;
    for(int i=0;i<node_num;i++) nwbc_ += (bcall[i] > 0);

    wbc.resize(nwbc_); nwbc_ = 0;
    for(int i=0;i<node_num;i++) {
      if(bcall[i] > 0) wbc[nwbc_++] = i+BASE;
    }

    // assign mesh statistics
    *nnodes = node_num;
    *xyz = node_x.data();
    *nwbc = nwbc_;
    *nobc = nobc_;
    *wbc_t = wbc.data();
    *obc_t = obc.data();
    *n4 = elem_counts[TYPE_TET];
    *n5 = elem_counts[TYPE_PYR];
    *n6 = elem_counts[TYPE_PRI];
    *n8 = elem_counts[TYPE_HEX];
    *ndc4 = element_node[TYPE_TET].data();
    *ndc5 = element_node[TYPE_PYR].data();
    *ndc6 = element_node[TYPE_PRI].data();
    *ndc8 = element_node[TYPE_HEX].data();
  }
}
/* ========================================================================== */


char ch_cap(char ch){
    if(97 <= ch && ch <= 122) ch -= 32;
    return ch;
}

bool ch_eqi(char ch1, char ch2){
    if(97 <= ch1 && ch1 <= 122) ch1 -= 32;
    if(97 <= ch2 && ch2 <= 122) ch2 -= 32;
    return (ch1 == ch2);
}

int ch_to_digit(char ch){
    int digit;

    if('0' <= ch && ch <= '9'){
        digit = ch - '0';
    } else
    if(ch == ' '){
        digit = 0;
    } else {
        digit = -1;
    }

    return digit;
}

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_DATA_READ reads data from a GMSH file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, character *GMSH_FILENAME, the GMSH filename.
//
//    Input, int NODE_DIM, the spatial dimension.
//
//    Input, double NODE_X[NODE_DIM*NODE_NUM], the node coordinates.
//
//    Input, int ELEMENT_ORDER, the order of the elements.
//
//    Input, int ELEMENT_NUM, the number of elements.
//
//    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM], 
//    the nodes that make up each element.
//
//****************************************************************************80
void gmsh_data_read(string gmsh_filename, int node_dim,
                    std::vector<double> &node_x,
                    std::vector<int> element_nodes[]){
    int i,j,k;
    int level;
    int length;
    int elem_type;
    int element_counts[TYPE_NUM_MAX] = {0};
    double x;
    bool ierror;
    ifstream input;
    string text;

    input.open(gmsh_filename.c_str());

    if (!input) {
        cerr << "\n";
        cerr << "GMSH_DATA_READ - Fatal error!\n";
        cerr << "  Could not open input file \"" << gmsh_filename << "\"\n";
        exit(1);
    }

    //
    // read node information
    //
    level = 0;
    for(;;){
        getline(input,text);
        if(input.eof()) break;

        if(level == 0){
            if(s_begin(text, "$Nodes")) level = 1;
        } else if(level == 1){
            s_to_i4(text, length, ierror);
            level = 2;
            j = 0;
        } else if(level == 2) {
            if(s_begin(text, "$EndNodes")){
                break;
            } else {
                s_to_i4(text, length, ierror);
                text.erase(0, length);
                for(i = 0; i < node_dim; i++) {
                    x = s_to_r8(text, length, ierror);
                    text.erase(0, length);
                    node_x[i+j*node_dim] = x;
                }
                j++;
            }
        }
    }

    //
    //  read element information.
    //
    level = 0;
    for(;;){
        getline(input,text);
        if(input.eof()) break;

        if(level == 0){
          if(s_begin(text, "$Elements")) level = 1;
        } else if(level == 1){
            // read element number
            s_to_i4(text, length, ierror);
            level = 2;
        } else if(level == 2){
            if(s_begin(text, "$EndElements")){
                break;
            } else {
                // read element preamble info
                int em_id = s_to_i4(text,length,ierror); text.erase(0,length); // element id
                elem_type = s_to_i4(text,length,ierror); text.erase(0,length); // element type
                            s_to_i4(text,length,ierror); text.erase(0,length); // number of tags (2)
                            s_to_i4(text,length,ierror); text.erase(0,length); // phys tag
                            s_to_i4(text,length,ierror); text.erase(0,length); // entity tag

                int element_count = element_counts[getElementType(elem_type)];
                int num_nodes = getElementNumNodes(elem_type);

                //printf("Filling Elem %d of type %d: count %d, num_nodes %d\n",em_id,elem_type,element_count,num_nodes);
                int * const elem_nodes = &(element_nodes[getElementType(elem_type)][num_nodes*element_count]);

                // read element nodes
                for(i = 0; i < num_nodes; i++){
                    k = s_to_i4(text, length, ierror);
                    text.erase(0, length);
                    //element_nodes[getElementType(elem_type)][element_order*element_count+i] = k;
                    elem_nodes[i] = k;
                }
                // increase element counter
                element_counts[getElementType(elem_type)]++;
            }
        }
    }
    input.close();
}

//****************************************************************************80
//
//  Purpose:
//
//    GMSH_SIZE_READ reads sizes from a GMSH file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    21 October 2014
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string GMSH_FILENAME, the GMSH filename.
//
//    Output, int &NODE_NUM, the number of nodes.
//
//    Output, int &NODE_DIM, the spatial dimension.
//
//    Output, int &ELEMENT_NUM, the number of elements.
//
//    Output, int &ELEMENT_ORDER, the order of the elements.
//
//****************************************************************************80
void gmsh_size_read(string gmsh_filename, int &node_num, int &node_dim,
                    int &element_num, int elem_counts[]){
    int elem_type;
    int level;
    int length;
    bool ierror;
    string text;
    ifstream input;
    const double r8_big = 1.0E+30;
    double x,y,z;
    double x_max;
    double x_min;
    double y_max;
    double y_min;
    double z_max;
    double z_min;

    node_num = 0;
    node_dim = 0;

    x_max = - r8_big;
    x_min = + r8_big;
    y_max = - r8_big;
    y_min = + r8_big;
    z_max = - r8_big;
    z_min = + r8_big;

    input.open(gmsh_filename.c_str());

    if(!input){
        cerr << "\n";
        cerr << "GMSH_SIZE_READ - Fatal error!\n";
        cerr << "  Could not open input file \"" << gmsh_filename << "\"\n";
        exit(1);
    }

    level = 0;
    for(;;){
        getline(input, text);
        if(input.eof()) break;

        if(level == 0){
            if(s_begin(text, "$MeshFormat")) level = 1;
        } else if(level == 1){
            if(s_begin(text, "$EndMeshFormat")) {
                break;
            } else {
                double fmt_a = s_to_r8(text, length, ierror);
                text.erase(0, length);

                if((int) 10*fmt_a >= 30){
                    cout << "ERROR - GMSH FORMAT 2.X ONLY: read format=" << fmt_a << endl;
                    exit(1);
                }
            }
        }
    }

    level = 0;
    for(;;){
        getline(input, text);
        if(input.eof()) break;

        if(level == 0){
            if(s_begin(text, "$Nodes")) level = 1;\
        } else if(level == 1){
            node_num = s_to_i4(text, length, ierror);
            level = 2;
        } else if(level == 2){
            if(s_begin(text, "$EndNodes")) {
                break;
            } else {
                s_to_i4(text, length, ierror);
                text.erase(0, length);

                x = s_to_r8(text, length, ierror);
                x_min = fmin(x_min, x);
                x_max = fmax(x_max, x);
                text.erase(0, length);

                y = s_to_r8(text, length, ierror);
                y_min = fmin(y_min, y);
                y_max = fmax(y_max, y);
                text.erase(0, length);

                z = s_to_r8(text, length, ierror);
                z_min = fmin(z_min, z);
                z_max = fmax(z_max, z);
                text.erase(0, length);
            }
        }
    }
    //
    //  Make a very simple guess as to the dimensionality of the data.
    //
    node_dim = 3;
    if(z_max == z_min){
        node_dim = 2;
        if(y_max == y_min) node_dim = 1;
    }

    //
    //  Now read element information.
    //
    level = 0;
    for(;;){
        getline(input, text);
        if(input.eof()) break;

        if(level == 0){
          if(s_begin(text, "$Elements")) level = 1;
        } else if(level == 1){
            // read element number
            element_num = s_to_i4(text, length, ierror);
            level = 2;
        } else if(level == 2){
            if(s_begin(text, "$EndElements")) {
                break;
            } else {
                // read element info: elem_id elem_type
                s_to_i4(text, length, ierror); text.erase(0, length);
                elem_type = s_to_i4(text, length, ierror); text.erase(0, length);

                // increase element type counter
                elem_counts[getElementType(elem_type)]++;

                for(;;){
                    s_to_i4(text, length, ierror);
                    text.erase(0, length);
                    if(ierror != 0) break;
                }
            }
        }
    }
    input.close();
}


//****************************************************************************80
//
//  Purpose:
//
//    S_BEGIN reports whether string 1 begins with string 2, ignoring case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    25 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, string S2, two strings.
//
//    Output, bool S_BEGIN, is true if S1 is the same as S2 up to
//    the end of S2, and false otherwise.
//
//****************************************************************************80
bool s_begin(string s1, string s2){
  int i;
  int n1;
  int n2;

  n1 = s1.length();
  n2 = s2.length();

  if(n1 < n2)
  {
    return false;
  }

  for(i = 0; i < n2; i++)
  {
    if(ch_cap(s1[i]) != ch_cap(s2[i]))
    {
      return false;
    }
  }
  return true;
}
//****************************************************************************80

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
//****************************************************************************80
int s_len_trim(string s){
  int n;

  n = s.length();

  while(0 < n) 
  {
    if(s[n-1] != ' ' && s[n-1] != '\n')
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_I4 reads an I4 from a string.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string to be examined.
//
//    Output, int &LAST, the last character of S used to make IVAL.
//
//    Output, bool &ERROR is TRUE if an error occurred.
//
//    Output, int *S_TO_I4, the integer value read from the string.
//    If the string is blank, then IVAL will be returned 0.
//
//****************************************************************************80
int s_to_i4(string s, int &last, bool &error){
  char c;
  int i;
  int isgn;
  int istate;
  int ival;

  error = false;
  istate = 0;
  isgn = 1;
  i = 0;
  ival = 0;

  for(; ;)
  {
    c = s[i];
    i = i + 1;
//
//  Haven't read anything.
//
    if(istate == 0)
    {
      if(c == ' ')
      {
      }
      else if(c == '-')
      {
        istate = 1;
        isgn = -1;
      }
      else if(c == '+')
      {
        istate = 1;
        isgn = + 1;
      }
      else if('0' <= c && c <= '9')
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        error = true;
        return ival;
      }
    }
//
//  Have read the sign, expecting digits.
//
    else if(istate == 1)
    {
      if(c == ' ')
      {
      }
      else if('0' <= c && c <= '9')
      {
        istate = 2;
        ival = c - '0';
      }
      else
      {
        error = true;
        return ival;
      }
    }
//
//  Have read at least one digit, expecting more.
//
    else if(istate == 2)
    {
      if('0' <= c && c <= '9')
      {
        ival = 10 * ival + c - '0';
      }
      else
      {
        ival = isgn * ival;
        last = i - 1;
        return ival;
      }

    }
  }
//
//  If we read all the characters in the string, see if we're OK.
//
  if(istate == 2)
  {
    ival = isgn * ival;
    last = s_len_trim(s);
  }
  else
  {
    error = true;
    last = 0;
  }

  return ival;
}

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int &LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool &ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
//**********************************************************************
double s_to_r8(string s, int &lchar, bool &error){
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim(s);
  error = false;
  r = 0.0;
  lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for(; ;)
  {
    c = s[lchar+1];
    lchar = lchar + 1;
//
//  Blank or TAB character.
//
    if(c == ' ' || c == TAB)
    {
      if(ihave == 2)
      {
      }
      else if(ihave == 6 || ihave == 7)
      {
        iterm = 1;
      }
      else if(1 < ihave)
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if(c == ',' || c == ';')
    {
      if(ihave != 1)
      {
        iterm = 1;
        ihave = 12;
        lchar = lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if(c == '-')
    {
      if(ihave == 1)
      {
        ihave = 2;
        isgn = -1;
      }
      else if(ihave == 6)
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if(c == '+')
    {
      if(ihave == 1)
      {
        ihave = 2;
      }
      else if(ihave == 6)
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if(c == '.')
    {
      if(ihave < 4)
      {
        ihave = 4;
      }
      else if(6 <= ihave && ihave <= 8)
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if(ch_eqi(c, 'E') || ch_eqi(c, 'D'))
    {
      if(ihave < 6)
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if(ihave < 11 && '0' <= c && c <= '9')
    {
      if(ihave <= 2)
      {
        ihave = 3;
      }
      else if(ihave == 4)
      {
        ihave = 5;
      }
      else if(ihave == 6 || ihave == 7)
      {
        ihave = 8;
      }
      else if(ihave == 9)
      {
        ihave = 10;
      }

      ndig = ch_to_digit(c);

      if(ihave == 3)
      {
        rtop = 10.0 * rtop +(double) ndig;
      }
      else if(ihave == 5)
      {
        rtop = 10.0 * rtop +(double) ndig;
        rbot = 10.0 * rbot;
      }
      else if(ihave == 8)
      {
        jtop = 10 * jtop + ndig;
      }
      else if(ihave == 10)
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if(iterm == 1 || nchar <= lchar + 1)
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if(iterm != 1 && lchar + 1 == nchar)
  {
    lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if(ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7)
  {
    error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if(jtop == 0)
  {
    rexp = 1.0;
  }
  else
  {
    if(jbot == 1)
    {
      rexp = pow(10.0, jsgn * jtop);
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow(10.0, rexp);
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}

//**********************************************************************
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    May 31 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
//**********************************************************************
void timestamp(){
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time(NULL);
  tm = localtime(&now);

  strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}