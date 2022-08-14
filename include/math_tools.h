#ifndef _MATH_TOOLS_H_
#define _MATH_TOOLS_H_

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
 * @file   math_tools.h
 * @brief  FEMBase Header
 * @author T. Otani
 */

#include "fem_define.h"

class mathTool {
 public:
    mathTool(){};
    ~mathTool(){};
    static double vectorNorm(const int &nump,std::vector<double> &x);
    static double innerProduct(const int &nump,std::vector<double> &x,std::vector<double> &y);
    static void crossProduct(const double (&a)[3],const double (&b)[3],double (&c)[3],double &dc);
    static void calcInverseMatrix_3x3(double (&inv_a)[3][3],const double (&a)[3][3]);
    static double calcDeterminant_3x3(const double (&a)[3][3]);
    static void calcMatrix_x_matrix4(double (&ans)[4][4],const double (&a)[4][4],const double (&b)[4][4]);
};
#endif //_MATH_TOOLS_H_
