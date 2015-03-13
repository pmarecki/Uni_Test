//
// Created by tttuuu on 3/12/15.
//
#ifndef _BUILDTESTS_MACROS_H_
#define _BUILDTESTS_MACROS_H_

#include <bits/stdc++.h>
#define REP(i,n)  for(int i=0;i<(int)(n);++i)
#define FOR(i,b,n)  for(int i=b;i<(n);++i)
#define ALL(c) (c).begin(),(c).end()
#define SS size()
#define CLR(a,v) memset((a),(v), sizeof a)
#define ST first
#define ND second
template<typename T, typename U> inline void AMIN(T &x, U y) { if(y < x) x = y; }
template<typename T, typename U> inline void AMAX(T &x, U y) { if(x < y) x = y; }
using namespace std;
typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<ll> vl;
typedef vector<vl> vvl;
typedef vector<double> vd;
typedef vector<vd> vvd;
typedef pair<int,int> pii;

template <typename T> void printContainer(T& a) {
  auto it = a.begin();
  cout << "{" << *(it++);
  for(; it!=a.end(); ++it) cout << ", " << (*it);
  cout << "}\n";
}

#endif //_BUILDTESTS_MACROS_H_
