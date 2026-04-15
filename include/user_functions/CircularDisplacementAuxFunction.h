/*------------------------------------------------------------------------*/
/*  Copyright 2026 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef CircularDisplacementAuxFunction_h
#define CircularDisplacementAuxFunction_h

#include <AuxFunction.h>

#include <vector>

namespace sierra{
namespace nalu{

class CircularHelper
{
 public:
  CircularHelper(
    const double height, const double sigmaNorm, const double sigmaTan,
    const double normXhat, const double normYhat, const double tanXhat, const double tanYhat,
    const double centerX, const double centerY, const double scaling);
  virtual ~CircularHelper() {};

  double increment(const double cX, const double cY);
  
 private:

  // internal parameters
  const double height_;
  const double sigmaSqNorm_;
  const double sigmaSqTan_;
  const double normXhat_;
  const double normYhat_;
  const double tanXhat_;
  const double tanYhat_;
  const double centerX_;
  const double centerY_;
  const double scaling_;
};

class CircularDisplacementAuxFunction : public AuxFunction
{
public:
  
  CircularDisplacementAuxFunction(
    const unsigned beginPos,
    const unsigned endPos,
    std::vector<double> theParams);

  virtual ~CircularDisplacementAuxFunction();
  
  virtual void do_evaluate(
    const double * coords,
    const double time,
    const unsigned spatialDimension,
    const unsigned numPoints,
    double * fieldPtr,
    const unsigned fieldSize,
    const unsigned beginPos,
    const unsigned endPos) const;
  
private:

  double timeBlending_;
  int totalOuter_;
  double heightOuter_;
  double diamOuter_;
  double sigmaNormOuter_;
  double sigmaTanOuter_;
  double declinationOuter_;
  int totalInner_;
  double heightInner_;
  double diamInner_;
  double sigmaNormInner_;
  double sigmaTanInner_;
  double declinationInner_;
  int totalInnerInner_;
  double heightInnerInner_;
  double diamInnerInner_;
  double sigmaNormInnerInner_;
  double sigmaTanInnerInner_;
  double declinationInnerInner_;
  double originX_;
  double originY_;
  double timeOffset_;
  const double pi_;
  const double Rxx_;
  const double Rxy_;
  const double Ryx_;
  const double Ryy_;

  std::vector<CircularHelper *> cHelpVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
