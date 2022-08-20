/**
 * @file math_tools.cpp
 * @brief mathTool class
 * @author T. Otani
 */

#include "math_tools.h"
#include <omp.h>
#include <string>
#include <iostream>
#include <cmath>

using namespace std;

// #################################################################
/**
 * @brief calc norm of vector
 */
double mathTool::vectorNorm(const int &nump,vector<double> &x)
{
  double norm=0e0;

  #pragma omp parallel for reduction(+:norm)
  for(int i=0;i<nump;i++){
    norm += fabs(x[i]);
  }
  return norm;
}

// #################################################################
/**
 * @brief calc inner product
 */
double mathTool::innerProduct(const int &nump,vector<double> &x,vector<double> &y)
{
  double dot_p=0e0;

  #pragma omp parallel for reduction(+:dot_p)
  for(int i=0;i<nump;i++){
    dot_p += x[i] * y[i];
  }
  return dot_p;
}

// #################################################################
/**
 * @brief calc cross product
 */
void mathTool::crossProduct(const double (&a)[3],const double (&b)[3],double (&c)[3],double &dc)
{
  dc = 0e0;

  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];

  for(int j=0;j<3;j++) dc += (c[j]*c[j]);

  dc = sqrt(dc);
}

// #################################################################
/**
 * @brief calc determinant
 */
double mathTool::calcDeterminant_3x3(const double (&a)[3][3])
{
  double det  = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2]
              - a[2][0] * a[1][1] * a[0][2] - a[1][0] * a[0][1] * a[2][2] - a[0][0] * a[2][1] * a[1][2];
  return det;
}

// #################################################################
/**
 * @brief calc inverse matrix
 */
void mathTool::calcInverseMatrix_3x3(double (&inv_a)[3][3],const double (&a)[3][3])
{
  double det;

  det = calcDeterminant_3x3(a);

  inv_a[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
  inv_a[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2];
  inv_a[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
  inv_a[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
  inv_a[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
  inv_a[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
  inv_a[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
  inv_a[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
  inv_a[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) inv_a[i][j] = inv_a[i][j] / det;
  }
}

void mathTool::calcMatrix_x_matrix4(double (&ans)[4][4],const double (&a)[4][4],const double (&b)[4][4])
{
  for(int i=0;i<4;i++){
    for(int j=0;j<4;j++){
      ans[i][j] = 0e0;
      for(int k=0;k<4;k++) ans[i][j] += a[i][k] * b[k][j];
    }
  }

}

