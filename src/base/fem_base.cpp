/**
 * @file fem_base.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem_base_mathTool.h"
#include "math_tools.h"

using namespace std;

// #################################################################
/**
 * @brief calc dudX
 * @param [in] stress
 */
void FEM_MathTool::calc_dudX(double (&dudX)[3][3],vector<vector<double>> &dNdX,vector<vector<double>> &u,const int &numOfNodeInElm)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      dudX[i][j] = 0e0;
      for(int p=0;p<numOfNodeInElm;p++) dudX[i][j] += dNdX[p][j] * u[p][i];
    }
  }
}

// #################################################################
/**
 * @brief calc dxdr
 * @param [in] stress
 */
void FEM_MathTool::calc_dxdr(double (&dxdr)[3][3],vector<vector<double>> &dNdr,vector<vector<double>> &x1,const int &numOfNodeInElm)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      dxdr[i][j] = 0e0;
      for(int p=0;p<numOfNodeInElm;p++){
        dxdr[i][j] += dNdr[p][j] * x1[p][i];
      }
    }
  }
}
// #################################################################
/**
 * @brief calc dXdr
 * @param [in] stress
 */
void FEM_MathTool::calc_dXdr(double (&dXdr)[3][3],vector<vector<double>> &dNdr,vector<vector<double>> &x0,const int &numOfNodeInElm)
{
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      dXdr[i][j] = 0e0;
      for(int p=0;p<numOfNodeInElm;p++){
        dXdr[i][j] += dNdr[p][j] * x0[p][i];
      }
    }
  }
}

// #################################################################
/**
 * @brief calc dNdx
 * @param [in] stress
 */
void FEM_MathTool::calc_dNdx(vector<vector<double>> &dNdx,vector<vector<double>> &dNdr,const double (&dxdr)[3][3],const int &numOfNodeInElm)
{
  double drdx[3][3];

  mathTool::calcInverseMatrix_3x3(drdx,dxdr);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      dNdx[p][i] = 0e0;
      for(int j=0;j<3;j++) dNdx[p][i] += dNdr[p][j] * drdx[j][i];
    }
  }
}

// #################################################################
/**
 * @brief calc dNdX
 * @param [in] stress
 */
void FEM_MathTool::calc_dNdX(vector<vector<double>> &dNdX,vector<vector<double>> &dNdr,const double (&dXdr)[3][3],const int &numOfNodeInElm)
{
  double drdX[3][3];

  mathTool::calcInverseMatrix_3x3(drdX,dXdr);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      dNdX[p][i] = 0e0;
      for(int j=0;j<3;j++){
        dNdX[p][i] += dNdr[p][j] * drdX[j][i];
      }
    }
  }
}