/*------------------------------------------------------------------------*/
/*  Copyright 2025 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/RotatingWallAuxFunction.h>
#include <Realm.h>
#include <SolutionOptions.h>

// basic c++
#include <algorithm>
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

RotatingWallAuxFunction::RotatingWallAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> theParams,
  Realm &realm) :
  AuxFunction(beginPos, endPos)
{
  // extract the params; always assume 3D since even in 2D we need a 3D omega
  if (theParams.size() != 6 )
    throw std::runtime_error("Rotating wall requires 6 parameter size");

  // set the data: for example, the following provides for an omega (rad/s) of 6 at a centroid of 0.0,0
  /*
       wall_user_data:
        user_function_name:
         velocity: rotating_wall
        user_function_parameters:
         velocity: [0.0, 0.0, 6.0, 0.0, 0.0, 0.0]   
  */
  for ( int i = 0; i < 3; ++i ) {
    omega_[i] = theParams[i];
    centroid_[i] = theParams[i+3];
  }
}

RotatingWallAuxFunction::~RotatingWallAuxFunction()
{
  // Nothing required 
}

void
RotatingWallAuxFunction::setup(const double time)
{
  // Nothing required
}

void
RotatingWallAuxFunction::do_evaluate(
  const double *coords,
  const double /*time*/,
  const unsigned /*spatialDimension*/,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  // temp data structures: radius, r and velocity, u
  double r[3] = {};
  double u[3] = {};
  
  for(unsigned p=0; p < numPoints; ++p) {

    // define radius vector; TBD
    
    // compute wall rotating velocity; TBD
    cross_product(r,u);
    
    // assign u to field pointer (velocity)
    for ( unsigned i = 0; i < fieldSize; ++i )
      fieldPtr[i] = u[i];

    // increment field pointer by the size
    fieldPtr += fieldSize;
    coords += fieldSize;
  }
}

void
RotatingWallAuxFunction::cross_product(double *r, double *u) const
{
  // TBD
}

} // namespace nalu
} // namespace Sierra
