/*------------------------------------------------------------------------*/
/*  Copyright 2025 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef RotatingWallAuxFunction_h
#define RotatingWallAuxFunction_h

#include <AuxFunction.h>

#include <string>
#include <vector>

namespace sierra{
namespace nalu{

class Realm;

class RotatingWallAuxFunction : public AuxFunction
{
public:

  RotatingWallAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    std::vector<double> theParams,
    Realm &realm);

  virtual ~RotatingWallAuxFunction();
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;

  void setup(const double time);
  void cross_product(double *c, double *u) const;

private:
  double omega_[3] = {};
  double centroid_[3] = {};
};

} // namespace nalu
} // namespace Sierra

#endif
