/*
 * STL_macros.h
 *
 *  Created on: Dec 15, 2012
 *      Author: fmarecki
 */

#ifndef STL_MACROS_H_
#define STL_MACROS_H_

#include <inttypes.h>
#include <cmath>
#include <cassert>

#include <iostream>
#include <fstream>
#include <sstream>

#include <vector>
#include <list>
#include <set>
#include <map>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <queue>
#include <stack>
#include <bitset>
#include <cstring>
using namespace std;


//Integer types
typedef uint64_t ulint;
typedef int64_t   lint;
typedef uint32_t uint;
typedef uint16_t usint;
typedef uint8_t  uchar;

//containers
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> pii;
typedef vector<pii> vii;

typedef set<int> si;
typedef multiset<int> mi;
typedef list<int> li;
typedef vector<string> vs;

//traversing macros (use google style guide)

#define REP(i,n)  for(int i=0;i<(n);++i)
#define REP1(i,n) for(int i=1;i<(n);++i)
#define FOR(i,l,u) for(int i=l;i<(u);++i)

//use `decltype` in >VC2010
#define TR(c,i) for(__typeof((c).begin()) i = (c).begin(); i != (c).end(); ++i)
#define ALL(c) (c).begin(),(c).end()
#define SS size()

#define SET_MIN(v,nv) if ((v) > (nv)) (v)=(nv);
#define SET_MAX(v,nv) if ((v) < (nv)) (v)=(nv);


//shortcuts
#define PB push_back
#define MP make_pair

#define TOP(s) (*s.begin() )  //no range checking
#define POP(s) s.erase(s.begin()) //no range checking
#define PUSH(s,k)  s.insert( (k) )

//#define PRESENT(c,x) ((c).find(x) != (c).end())
//#define CPRESENT(c,x) (find(ALL(c),x) != (c).end())




const int INF = 1e9;

//Arithmetics

// 1/5 --> 1;   5/5-->1;   6/5-->2
inline int CeilDiv(int a, int b) {
  int r = (a % b)? 1: 0;
  return a / b + r;   ///  !!!! standard stuff == (a+b-1)/b;
}

//Returns (-1)^k
inline int SIGN(int k) {
  return (k % 2)? -1 : 1;
}



inline void EEE() {
  printf("\n");
}

inline void assign_less(int *reg, int val) {
  if (val < *reg) *reg = val;
}


template <typename  T>
void PrintContainer(const T w) {
  printf("{");
  TR(w, i) {
    if (*i != INF)
      printf("%i, ", *i);
    else
      printf("-, ");
  }
  printf("\b\b}\n");
}


#define PrintMatrix(G,N,M) \
  REP(i,N) {                                        \
      printf("[%2i] ",i);                           \
      REP(j,N)                                      \
        printf("%3d", (G[i][j]==INF)?-1:G[i][j]);   \
      EEE();                                        \
    }

#define PrintDMatrix(G,N,M) \
  REP(i,N) {                                        \
      printf("[%2i] ",i);                           \
      REP(j,N)                                      \
        printf("%3.2f ", (G[i][j]==INF)?-1:G[i][j]);   \
      EEE();                                        \
    }




void PrintBit(uint H) {
  REP(i, 32) {
    int no = H;
    int bit = 31-i%32;
    bool val = no & (1<<bit);
    printf("%i%s%s",val,(i+1)%8==0?" ":"",(i+1)%32==0?"| ":"");
  }
  EEE();
}

void PrintBit(ulint H) {
  REP(i, 64) {
    ulint no = H;
    int bit = 63 - i%64;
    bool val = (no & (1UL<<bit));
//    printf("(%i),", bit);
    printf("%i%s%s", val, ((i+1)%8==0)?" ":"",(i+1)%64==0?"| ":"");
  }
//  EEE();
}



uint lrot32(uint x, uint n) {
  return ((x << n) | (x >> (32-n)));
}
uint rrot32(uint x, uint n) {
  return ((x >> n) | (x << (32-n)));
}

//Returns lowest byte of a pointer.
inline uint32_t Lptr(void *ptr) { return (uint64_t)(ptr) & 0xfff;}



#endif /* STL_MACROS_H_ */
