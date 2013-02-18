#include "STL_MACRO.h"


void matrix2d() {
  int sizex = 20, sizey = 10;
  int **a;
  a = new int*[sizex];              // allows a[0...sizex-1];
  a[0] = new int[sizex * sizey];    // full allocation in one piece

  //instead of doing a[0] = new int[sizey], a[1] = new int[sizey] ...
  for(int x=1; x<sizex; ++x)
    a[x] = a[x-1] /*adres of previous*/ + sizey /*stride of block*/ ;
  //this shifts "a[2]" by "sizeof(typeof(a[1])) * sizey"

  REP(i, sizex)
    REP(j, sizey)
      a[i][j] = 1000 * i + j;   //filling

  REP(i, sizex)
      REP(j, sizey)
        assert(a[i][j] == 1000 * i + j);   //reading, no overwrite
  delete a;
}



void matrix3d() {
  int sizex = 20, sizey = 10, sizez=5;

  int ***a;
  a = new int**[sizex];
  //Normal: a[0...sizex] = new int*[sizey];
  //        a[0...sizex][0...sizey] = new int[sizez]
  //Gives non-contiguous memory [0][1] -> ptr, not adjacent to [1][2]->ptr.
  //
  //Idea -- first declare one large piece of memory, and then adjust pointers.

  a[0]    = new int*[sizey * sizex];    //chop-off later for a[1..sizex]
  a[0][0] = new int [sizez * sizex * sizey];    //true memory
  //Allowed operation on a[0][0][0...sizez-1] ... to "sizex*sizey*sizez-1".

  //Goal of setup: fill a[x][y], for all values of x,y, with appropriate pointers.

  //so far: "a" is a number
  //        "a[0]" is a number
  //        "a[1]" is rubbish
  //        "a[0][0]" is a pointer to first element of z-vector
  //        "a[0][1]" is rubbish

  //form first row, as we know "a[0][0]"
  for(int y=1; y<sizey; ++y)
    a[0][y] = a[0][y-1] + sizez;

  //        each "a[0][i]" now points to first element in z-vector

  //Generalize  a[1] = a[0] + sizey; for all "x"; just poiters to
  // interesting elements.
  for(int x=1; x<sizex; ++x)
    a[x] = a[x-1] + sizey;

  //Fill first elements of each "a[i]"
  for(int x=1; x<sizex; ++x) {
    a[x][0] = &(a[x-1][sizey-1][sizez-1])+1;      //need pointer to z-vector
    for(int y=1; y<sizey; ++y)
      a[x][y] = a[x][y-1] + sizez;
  }

  //Write
  REP(x, sizex)
    REP(y, sizey)
      REP(z, sizez)
        a[x][y][z] = 100000*x + 1000*y + z;
  //Read and check overwrites
  REP(x, sizex)
    REP(y, sizey)
      REP(z, sizez)
        assert(a[x][y][z] ==  100000*x + 1000*y + z);

}


int main(int cc, char **vv) {
  matrix3d();
}
