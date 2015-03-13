//
// Created by tttuuu on 3/12/15.
//

#ifndef _BUILDTESTS_EIGENINTERFACE_H_
#define _BUILDTESTS_EIGENINTERFACE_H_

#include "macros.h"

class EigenInterface {

public:
  vd eigenVectors(vd matrix2d);
  void diagonalize5x5();

  void testDecompositionMatrices();

};


#endif //_BUILDTESTS_EIGENINTERFACE_H_
