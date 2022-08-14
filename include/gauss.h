#ifndef _GAUSS_H_
#define _GAUSS_H_


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
 * @file   gauss.h
 * @brief  FEMBase Header
 * @author T. Otani
 */
#include "fem_define.h"
#include <cstdio>

class Gauss{
public:
	Gauss(){}
	Gauss(const int num){
		switch(num){
			case 1:
	  		point[0] = -0.577350296189626; point[1] =  0.577350296189626;
      	weight[0] = 1e0;							 weight[1] = 1e0;
      	break;
      case 2:
      	point[0] = -0.774596669241483; point[1] = 0e0; point[2] = 0.774596669241483;
      	weight[0] = 0.555555555555555; weight[1] = 0.888888888888888; weight[2] = 0.555555555555555;
      	break;
      case 3:
      	point[0] = -0.861135311594053;  point[1] = -0.339981043584856; point[2] = 0.339981043584856;  point[3] = 0.861135311594053;
      	weight[0] = 0.347854845137454;  weight[1] = 0.652145154862546; weight[2] = 0.652145154862546; weight[3] = 0.347854845137454;
      break;
      default:
      	printf("undefined order is set in the gauss integral\n");
		}
	}
  double point[4],weight[4];
};

//caution! Below is the position of natural coordinates (L0,L1,L2,L3)
class GaussWedge{
public:
	GaussWedge(){}
	GaussWedge(const int num){
    double a,b,c;
	  switch (num){
	    case 1:
        a= 0.166666667;
	      b = 0.577350269;
	      c = 0.666666667;
	  	  point[0][0] = a;  point[0][1] = a;  point[0][2] = -b;
	  	  point[1][0] = c;  point[1][1] = a;  point[1][2] = -b;
	  	  point[2][0] = a;  point[2][1] = c;  point[2][2] = -b;
	  	  point[3][0] = a;  point[3][1] = a;  point[3][2] = b;
	  	  point[4][0] = c;  point[4][1] = a;  point[4][2] = b;
	  	  point[5][0] = a;  point[5][1] = c;  point[5][2] = b;
        weight[0] = 1e0;
        break;
      default:
        	printf("undefined order is set in the gauss integral\n");
	  }
	}
  double point[6][3],weight[5];
};

class GaussTetra{
public:
	GaussTetra(){}
	GaussTetra(const int num){
    double a;
		switch(num){
			case 1:
	  		point[0][0] = 2.5e-1; point[0][1] = 2.5e-1; point[0][2] = 2.5e-1; point[0][3] = 2.5e-1;
      	weight[0] = 1e0;
      	break;
      case 2:
        a=0.13819660;
	  		point[0][0] = 1e0-2e0*a; point[0][1] = a; point[0][2] = a; point[0][3] = a;
	  		point[1][0] = a; point[1][1] = 1e0-2e0*a; point[1][2] = a; point[1][3] = a;
	  		point[2][0] = a; point[2][1] = a; point[2][2] = 1e0-2e0*a; point[2][3] = a;
	  		point[3][0] = a; point[3][1] = a; point[3][2] = a; point[3][3] = 1e0-2e0*a;
      	weight[0] = 2.5e-1;weight[1] = 2.5e-1;weight[2] = 2.5e-1;weight[3] = 2.5e-1;
      	break;
      case 3:
	  		point[0][0] = 2.5e-1; point[0][1] = 2.5e-1; point[0][2] = 2.5e-1; point[0][3] = 2.5e-1;
	  		point[1][0] = 5e-1;    point[1][1] = 1e0/6e0; point[1][2] = 1e0/6e0; point[1][3] = 1e0/6e0;
	  		point[2][0] = 1e0/6e0; point[2][1] = 5e-1;    point[2][2] = 1e0/6e0; point[2][3] = 1e0/6e0;
	  		point[3][0] = 1e0/6e0; point[3][1] = 1e0/6e0; point[3][2] = 5e-1;    point[3][3] = 1e0/6e0;
	  		point[4][0] = 1e0/6e0; point[4][1] = 1e0/6e0; point[4][2] = 1e0/6e0; point[4][3] = 5e-1;
      	weight[0] = -8e-1;weight[1] = 9e0/20e0;weight[2] = 9e0/20e0;weight[3] = 9e0/20e0;weight[4] = 9e0/20e0;
      break;
      default:
      	printf("undefined order is set in the gauss integral\n");
		}
	}
  double point[5][4],weight[5];
};

class GaussTriangle{
public:
  GaussTriangle(){}
  GaussTriangle(const int num){
    switch(num){
      case 1:
        point[0][0] = 1e0/3e0; point[0][1] = 1e0/3e0; point[0][2] = 1e0/3e0;
        weight[0] = 1e0;
        break;
      case 2:
        point[0][0] = 0e0; point[0][1] = 5e-1; point[0][2] = 5e-1;
        point[1][0] = 5e-1; point[1][1] = 0e0;point[1][2] = 5e-1;
        point[2][0] = 5e-1; point[2][1] = 5e-1; point[2][2] = 0e0;
        weight[0] = 1e0/3e0;weight[1] = 1e0/3e0;weight[2] = 1e0/3e0;
        break;
      case 3:
        point[0][0] = 1e0/3e0;  point[0][1] = 1e0/3e0;  point[0][2] = 1e0/3e0;
        point[1][0] = 11e0/15e0;point[1][1] = 2e0/15e0; point[1][2] = 2e0/15e0;
        point[2][0] = 2e0/15e0; point[2][1] = 11e0/15e0;point[2][2] = 2e0/15e0;
        point[3][0] = 2e0/15e0; point[3][1] = 2e0/15e0; point[3][2] = 2e0/15e0;
        weight[0] = -27e0/48e0;weight[1] = 25e0/48e0;weight[2] = 25e0/48e0;weight[3] = 25e0/48e0;
      break;
      default:
        printf("undefined order is set in the gauss integral\n");
    }
  }
  double point[4][3],weight[4];
};

#endif //_GAUSS_H_
