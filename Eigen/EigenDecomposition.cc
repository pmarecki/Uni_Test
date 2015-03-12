
void EigenInterface::doIt() {
//  MatrixXd A = MatrixXd::Random(6,6);
//  cout << "Here is a random 6x6 matrix, A:" << endl << A << endl << endl;
//  EigenSolver<MatrixXd> es(A);
//  cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
//  cout << "The matrix of eigenvectors, V, is:" << endl << es.eigenvectors() << endl << endl;

  MatrixXd m(5,5);
  m << 1, 2, 3, 4, 5,
       2, 1, 1, 1, 1,
       3, 1, 1, 1, 1,
       4, 1, 1, 1, 1,
       5, 1, 1, 1, 1;
  cout << m << endl;
  EigenSolver<MatrixXd> roz(m);
  cout << "Eigenvalues:\n" << roz.eigenvalues() << endl;
  cout << "Eigenvectors:\n" << roz.eigenvectors() << endl;
}
