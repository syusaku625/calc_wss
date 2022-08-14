#ifndef _FEM_BASE_MATH_TOOL_H_
#define _FEM_BASE_MATH_TOOL_H_


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
 * @file   fem_base_mathTooL.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */
#include<vector>

class FEM_MathTool{
  public:
  static void calc_dudX(double (&dudX)[3][3],std::vector<std::vector<double>> &dNdX,std::vector<std::vector<double>> &u,const int &numOfNodeInElm);
  static void calc_dXdr(double (&dXdr)[3][3],std::vector<std::vector<double>> &dNdr,std::vector<std::vector<double>> &x0,const int &numOfNodeInElm);
  static void calc_dxdr(double (&dxdr)[3][3],std::vector<std::vector<double>> &dNdr,std::vector<std::vector<double>> &x,const int &numOfNodeInElm);
  static void calc_dNdX(std::vector<std::vector<double>> &dNdX,std::vector<std::vector<double>> &dNdr,const double (&dXdr)[3][3],const int &numOfNodeInElm);
  static void calc_dNdx(std::vector<std::vector<double>> &dNdx,std::vector<std::vector<double>> &dNdr,const double (&dxdr)[3][3],const int &numOfNodeInElm);
  static void tensorPushForward_4order(double (&c4)[3][3][3][3],const double (&C4)[3][3][3][3],const double (&F)[3][3],const double J);

};

#endif