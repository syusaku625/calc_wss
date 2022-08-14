#ifndef _SHAPE_FUNCTION_H_
#define _SHAPE_FUNCTION_H_


//##################################################################################
//
// FEM Base
//
// Copyright (c) 2016 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   shape_function.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */

#include "fem_define.h"
using namespace std;
class ShapeFunction3D{
 public:
  static void C3D6_N(vector<double> &N,const double &L1,const double &L2,const double &L3,const double &g1)
  {
    N[0] = 5e-1 * (1e0 - L1 - L2) * (1e0 - L3);
    N[1] = 5e-1 * L1 * (1e0 - L3);
    N[2] = 5e-1 * L2 * (1e0 - L3);
    N[3] = 5e-1 * (1e0 - L1 - L2) * (1e0 + L3);
    N[4] = 5e-1 * L1 * (1e0 + L3);
    N[5] = 5e-1 * L2 * (1e0 + L3);
  }

  static void C3D6_dNdr(vector<vector<double>> &dNdr,const double &L1,const double &L2,const double &L3,const double &g1)
  {
    dNdr[0][0] = -5e-1 * (1e0 - L3);
    dNdr[0][1] = -5e-1 * (1e0 - L3);
    dNdr[0][2] = -5e-1 * (1.0 - L1 - L2);
    dNdr[1][0] = 5e-1*(1e0-L3);
    dNdr[1][1] = 0e0;
    dNdr[1][2] = -5e-1*L1;
    dNdr[2][0] = 0e0;
    dNdr[2][1] = 5e-1*(1e0-L3);
    dNdr[2][2] = -5e-1*L2;
    dNdr[3][0] = -5e-1*(1e0+L3);
    dNdr[3][1] = -5e-1*(1e0+L3);
    dNdr[3][2] = 5e-1*(1.0-L1-L2);
    dNdr[4][0] = 5e-1*(1e0+L3);
    dNdr[4][1] = 0e0;
    dNdr[4][2] = 5e-1*L1;
    dNdr[5][0] = 0e0;
    dNdr[5][1] = 5e-1*(1e0+L3);
    dNdr[5][2] = 5e-1*L2;
  }
};

#endif //_SHAPE_FUNCTION_H_