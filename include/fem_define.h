#ifndef _FEM_DEFINE_H_
#define _FEM_DEFINE_H_

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
 * @file   fem_define.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */
#include <string>
#include <vector>
#include "vtkCellType.h"

template<typename T>
using VECTOR1D = std::vector<T>;
template<typename T>
using VECTOR2D = std::vector<std::vector<T>>;
template<typename T>
using VECTOR3D = std::vector<std::vector<std::vector<T>>>;
template<typename T>
using VECTOR4D = std::vector<std::vector<std::vector<std::vector<T>>>>;
template<typename T>
using VECTOR5D = std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>;
template<typename T>
using VECTOR6D = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>>;
template<typename T>
using VECTOR7D = std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<std::vector<T>>>>>>>;

#define FEM_VERS "1.0"
constexpr double GRAVITY = 9.80665; // (m/s2)
constexpr double PI = 3.1415926535897932384626;

// general
#define ON          1
#define OFF         0

//data encode
#define INT         0
#define DOUBLE      1
#define ASCII       0
#define BINARY      1

//Constitutive law
#define LinearElasticMaterial 0
#define StVenantMaterial 1
#define NeoHookean 2

// Mesh type
#define C2D4        0
#define C3D4        1
#define C3D8        2

// IO file format
#define ORIGINAL    0
#define HDF5        1

#define WALL 0
#define INLET 1
#define OUTLET 2

#define XAXIS 0
#define YAXIS 1
#define ZAXIS 2

#define VELOCITY 0
#define PRESSURE 1

class ElementType{
 public:
  VTKCellType meshType;
  int materialType;
  int numOfGaussPoint;
  VECTOR1D<int> node;
};


#endif // _FB_DEFINE_H_
