#include <cstdio>
#include <utility>
using std::pair;
#define MP make_pair




/**
 * 3D Matrix from numerical recipes. 
 * Storage "v" defined in peculiar way. 
 * 
 */
template<class T>
class NRMat3d {
public:
  NRMat3d();
  NRMat3d(int sizex, int sizey, int sizez);
  inline T** operator[](const int i);  //subscripting: pointer to row i
  inline const T* const * operator[](const int i) const;
  inline int dim1() const;
  inline int dim2() const;
  inline int dim3() const;

  //Our addition -- returning pointer to the raw matrix, and x,y size.
  T* RawMatrix2D() {
    return &(v[0][0][0]);
  }

  ~NRMat3d();
private:
  int nn;  //dimension x
  int mm;  //dimension y
  int kk;  //dimension z
  T ***v;
};

template<class T>
NRMat3d<T>::NRMat3d() :
    nn(0), mm(0), kk(0), v(NULL) {
}



/**
 * Verbose version of main constructor.
 */
template<class T>
NRMat3d<T>::NRMat3d(int sizex, int sizey, int sizez) :
    nn(sizex), mm(sizey), kk(sizez) {
  /**
   * Idea: _) allocate space first (T[sizex*sizey*sizez])
   *       _) update pointers
   * Result: memory is contiguous, and addressing via v[1][3][17] works.
   */
  v = new T**[sizex];             //for each "x" there is an [yz] matrix
  v[0]    = new T*[sizex*sizey];  //2d matrix (face), from which vectors(z) grow
  v[0][0] = new T[sizex*sizey*sizez];  //TRUE ALLOCATION of memory

  // Putting all z-vectors (aka v[i][j] together)
  
  //border-case, x=0, ---along y
  for(int j=1; j<sizey; j++)
    v[0][j] = v[0][j-1] + sizez;      //z-vectors, j=1 takes from v[0][0]
  
  for(int i=1; i<sizex; i++) {
    v[i] = v[i-1] + sizey;                //????
    v[i][0] = v[i-1][0] + sizey * sizez;  //border case, y=0, ---along x

    //bulk of the matrix, v[i][j] is pointer to "T[size_z]"
    for(int j=1; j<sizey; j++)
      v[i][j] = v[i][j-1] + sizez;
  }
}

//subscripting: pointer to row i
template<class T>
inline T** NRMat3d<T>::operator[](const int i) {
  return v[i];
}

template<class T>
inline const T* const* NRMat3d<T>::operator[](const int i) const {
  return v[i];
}

template<class T>
inline int NRMat3d<T>::dim1() const {
  return nn;
}

template<class T>
inline int NRMat3d<T>::dim2() const {
  return mm;
}

template<class T>
inline int NRMat3d<T>::dim3() const {
  return kk;
}

template<class T>
NRMat3d<T>::~NRMat3d() {
  if (v != NULL) {
    delete[] (v[0][0]);
    delete[] (v[0]);
    delete[] (v);
  }
}


