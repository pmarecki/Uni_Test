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
typedef pair<int,int> pii;

#include <Eigen/Dense>


using namespace Eigen;
//using Eigen::Matrix2d;
//using Eigen::VectorXcd;

void halloWorld() {
  Matrix2d m(2,2);
  m(0,0) = 4;
  m(1,0) = 2;
  m(0,1) = 1;
  m(1,1) = 1;
  cout << m << endl;
  cout << m.transpose() << endl;

  VectorXcd eigenval = m.eigenvalues();
  cout << eigenval << endl;
  REP(i, eigenval.size()) {
    cout << eigenval[i].real() << endl;
  }

  m << 4, 2,
       2, 3;


  cout << "--------\n";
  cout << (m.diagonal()) << endl;
  cout << m.norm() << endl;
}


int main() {
//  halloWorld();
  
  MatrixXf m(2,2);
  m(0,0) = 1;
  m(1,0) = 1;
  m(0,1) = 1;
  m(1,1) = 1;

  JacobiSVD<MatrixXf> svd(m, ComputeThinU | ComputeThinV);
  cout << "U=" << svd.matrixU() << endl;
  cout << "V=" << svd.matrixV().transpose() << endl;
  cout << "S=" << svd.singularValues() << endl;

//  MatrixXf check = svd.matrixU() * svd.singularValues() * svd.matrixV();
  MatrixXf check = svd.matrixU() *  svd.matrixU().transpose();
  cout << (check) << endl;

//  Vector2d w(1,2);
//  cout << w << endl;
//  cout << w.asDiagonal() << endl;
}


//include_directories(~/ClionProjects/Libraries/eigen)
