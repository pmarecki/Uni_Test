//
// Created by tttuuu on 3/12/15.
//

#include "EigenInterface.h"

#include <Eigen/Dense>
using namespace Eigen;

vd EigenInterface::eigenVectors(vd mm) {
  Matrix2d m(2,2);
  m(0,0) = mm[0];
  m(1,0) = mm[1];
  m(0,1) = mm[2];
  m(1,1) = mm[3];
  cout << m << endl;
  VectorXcd eigenval = m.eigenvalues();
  cout << "Eigenvalues: "<< endl;
  cout << eigenval << endl;
  cout << "Diagonal:" << endl;
  cout << m.diagonal() << endl;

  Vector2d w = {1,2};
  Matrix2d mu = w.asDiagonal();

  cout << "Diagonal matrix from: " << w << endl;
  cout << mu << "||||\n";

  //Jacobi decomposition: projecting-out the positive part
  JacobiSVD<Matrix2d> svd(m, ComputeFullU | ComputeFullV);
  auto sing = svd.singularValues();
  Matrix2d U = svd.matrixU();
  Vector2d sin = svd.singularValues();
  cout << "Printing the vectof of singular values:\n";
  cout << sin << endl;
  cout << "-------------\n";
  cout << "Printing U:\n";
  cout << svd.matrixU() << endl;
  cout << "Printing V:\n";
  cout << svd.matrixU() << endl;
  cout << "--------------\n";
  EigenSolver<Matrix2d> roz(m);
  cout << "Eigenvalues:\n" << roz.eigenvalues() << endl;
  cout << "Eigenvectors:\n" << roz.eigenvectors() << endl;
//  cout << "Singular values from decomposed m:" << wsing << endl;

  //Type casting:
//  Eigen::MatrixXd d;                       // Matrix of doubles.
//  Eigen::MatrixXf f = d.cast <float> ();   // Matrix of floats.

  return vd();
}

void EigenInterface::diagonalize5x5() {
//  MatrixXd A = MatrixXd::Random(6,6);
//  cout << "Here is a random 6x6 matrix, A:" << endl << A << endl << endl;
//  EigenSolver<MatrixXd> es(A);
//  cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
//  cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;

  MatrixXd m(5,5);
  m << 1, 2, 3, 4, 5,
       3, 1, 1, 1, 1,
       2, 1, 1, 1, 1,
       4, 1, 1, 1, 1,
       5, 1, 1, 1, 1;
  cout << m << endl;
  EigenSolver<MatrixXd> roz(m);
  cout << "Eigenvalues:\n" << roz.eigenvalues() << endl;
  cout << "Eigenvectors:\n" << roz.eigenvectors() << endl;
}


void EigenInterface::testDecompositionMatrices() {
  Matrix2d m;
  m << 0,1,
       1,0;
//  m = m * m;
  EigenSolver<Matrix2d> roz(m);
  cout << "Eigenvalues\n" << roz.eigenvalues() <<endl;
  cout << "Eigenvectors:\n" << roz.eigenvectors() << endl;
  Matrix2cd U = roz.eigenvectors();
  Vector2cd ev = roz.eigenvalues();
  Matrix2cd EV = ev.asDiagonal();
  Matrix2d EVreal = EV.real();
  cout << "EV matrix:\n" << EV << endl;
  //Cut-out negative EV's
  double eps = 0.000001;
  REP(i,2) if (EVreal(i,i) < eps) EVreal(i,i) = eps;

  Matrix2cd M = U * EVreal * U.adjoint();
  cout << "re-transformed:\n" << (M.real()) << endl;

  Matrix2d Mreal = M.real();

//  REP(i,2) EV(i, i) = max(EV(i, i), 0.00001);
//  cout << "positive-part" << (U * EV * U.adjoint()) << endl;



//  cout << "U=\n" << U << endl;
}

























