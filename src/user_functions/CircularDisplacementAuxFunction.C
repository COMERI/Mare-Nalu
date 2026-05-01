/*------------------------------------------------------------------------*/
/*  Copyright 2026 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <user_functions/CircularDisplacementAuxFunction.h>
#include <algorithm>
#include <PecletFunction.h>

#include <NaluEnv.h>

// basic c++
#include <cmath>
#include <vector>
#include <stdexcept>

namespace sierra{
namespace nalu{

CircularDisplacementAuxFunction::CircularDisplacementAuxFunction(
  const unsigned beginPos,
  const unsigned endPos,
  const std::vector<double> theParams) :
  AuxFunction(beginPos, endPos),
  timeBlending_(1.0),
  totalOuter_(23),
  heightOuter_(5.0),
  diamOuter_(100.0),
  sigmaNormOuter_(5.0),
  sigmaTanOuter_(3.0),
  declinationOuter_(0.0),
  totalInner_(23),
  heightInner_(3.0),
  diamInner_(60.0),
  sigmaNormInner_(3.0),
  sigmaTanInner_(2.0),
  declinationInner_(0.0),
  totalInnerInner_(23),
  heightInnerInner_(3.0),
  diamInnerInner_(60.0),
  sigmaNormInnerInner_(3.0),
  sigmaTanInnerInner_(2.0),
  declinationInnerInner_(0.0),  
  icElevation_(0.0),
  icWidth_(0.0),
  icTrans_(0.0),
  originX_(0.0),
  originY_(0.0),
  timeOffset_(0.0),
  pi_(acos(-1.0)),
  Rxx_(std::cos(pi_/2.0)),
  Rxy_(-std::sin(pi_/2.0)),
  Ryx_( std::sin(pi_/2.0)),
  Ryy_(std::cos(pi_/2.0))
{

  // parse the required parameters
  if ( theParams.size() != 25 ) {
    NaluEnv::self().naluOutputP0() << "Error: circular user function requires 25 parameters" << std::endl;
    throw std::runtime_error("CircularDisplacementAuxFunction::Error");
  }

  // allow for linear ramp up in time
  timeBlending_ = theParams[0];
  
  // outer....
  totalOuter_ = theParams[1];
  heightOuter_ = theParams[2];
  diamOuter_ = theParams[3];
  sigmaNormOuter_ = theParams[4];
  sigmaTanOuter_ = theParams[5];
  
  // inner....
  totalInner_ = theParams[6];
  heightInner_ = theParams[7];
  diamInner_ = theParams[8];
  sigmaNormInner_ = theParams[9];
  sigmaTanInner_ = theParams[10];
  
  // inner, inner....
  totalInnerInner_ = theParams[11];
  heightInnerInner_ = theParams[12];
  diamInnerInner_ = theParams[13];
  sigmaNormInnerInner_ = theParams[14];
  sigmaTanInnerInner_ = theParams[15];
  
  // origin of structure
  originX_ = theParams[16];
  originY_ = theParams[17];
  
  // offset
  timeOffset_ = theParams[18];
  
  // declinations (specified in degree)
  declinationOuter_ = theParams[19];
  declinationInner_ = theParams[20];
  declinationInnerInner_ = theParams[21];

  // initial elevation
  icElevation_ = theParams[22];
  icWidth_ = theParams[23];
  icTrans_ = theParams[24];
  
  // define radius
  const double radOuter = diamOuter_*0.5;
  const double radInner = diamInner_*0.5;
  const double radInnerInner = diamInnerInner_*0.5;
  
  //=========================================
  // create the vector of outer structures
  //=========================================
  const double dThetaOuter = 2.0*pi_/(double)totalOuter_;
  double thetaOuter = declinationOuter_*pi_/180.0;
  for ( int k = 0; k < totalOuter_; ++k ) {
    
    // on the unit circle, diameter 1.0
    double cX = std::cos(thetaOuter);
    double cY = std::sin(thetaOuter);
    double mag = std::sqrt(cX*cX + cY*cY);
    
    // unit normal
    double normXhat = cX/mag;
    double normYhat = cY/mag;
    
    // now rotate 90 degrees for a tangential unit vector 
    double tanXhat = Rxx_*normXhat + Rxy_*normYhat;
    double tanYhat = Ryx_*normXhat + Ryy_*normYhat;
    
    // center on the scaled system at the specified origin
    double centerX = originX_ + cX*radOuter;
    double centerY = originY_ + cY*radOuter;
    
    // create the helper
    CircularHelper *cHelp = new CircularHelper(heightOuter_, sigmaNormOuter_, sigmaTanOuter_,
					       normXhat, normYhat, tanXhat, tanYhat,
					       centerX, centerY, 1.0);
    cHelpVec_.push_back(cHelp);
    
    // incremenet theta
    thetaOuter += dThetaOuter;
  }
  
  //=========================================
  // create the vector of inner structures
  //=========================================
  const double dThetaInner = 2.0*pi_/(double)totalInner_;
  double thetaInner = declinationInner_*pi_/180.0;
  for ( int k = 0; k < totalInner_; ++k ) {
    
    // on the unit circle, diameter 1.0
    double cX = std::cos(thetaInner);
    double cY = std::sin(thetaInner);
    double mag = std::sqrt(cX*cX + cY*cY);
    
    // unit normal
    double normXhat = cX/mag;
    double normYhat = cY/mag;
    
    // now rotate 90 degrees for a tangential unit vector 
    double tanXhat = Rxx_*normXhat + Rxy_*normYhat;
    double tanYhat = Ryx_*normXhat + Ryy_*normYhat;
    
    // center on the translated system
    double centerX = originX_ + cX*radInner;
    double centerY = originY_ + cY*radInner;
    
    // create the helper
    CircularHelper *cHelp = new CircularHelper(heightInner_, sigmaNormInner_, sigmaTanInner_,
					       normXhat, normYhat, tanXhat, tanYhat,
					       centerX, centerY, 1.0);
    cHelpVec_.push_back(cHelp);
    
    // incremenet theta
    thetaInner += dThetaInner;    
  }

  //=========================================
  // create the vector of inner (inner) structures
  //=========================================
  const double dThetaInnerInner = 2.0*pi_/(double)totalInnerInner_;
  double thetaInnerInner = declinationInnerInner_*pi_/180.0;
  for ( int k = 0; k < totalInnerInner_; ++k ) {
    
    // on the unit circle, diameter 1.0
    double cX = std::cos(thetaInnerInner);
    double cY = std::sin(thetaInnerInner);
    double mag = std::sqrt(cX*cX + cY*cY);
    
    // unit normal
    double normXhat = cX/mag;
    double normYhat = cY/mag;
    
    // now rotate 90 degrees for a tangential unit vector 
    double tanXhat = Rxx_*normXhat + Rxy_*normYhat;
    double tanYhat = Ryx_*normXhat + Ryy_*normYhat;
    
    // center on the translated system
    double centerX = originX_ + cX*radInnerInner;
    double centerY = originY_ + cY*radInnerInner;
    
    // create the helper
    CircularHelper *cHelp = new CircularHelper(heightInnerInner_, sigmaNormInnerInner_, sigmaTanInnerInner_,
					       normXhat, normYhat, tanXhat, tanYhat,
					       centerX, centerY, -1.0);
    cHelpVec_.push_back(cHelp);
    
    // incremenet theta
    thetaInnerInner += dThetaInnerInner;    
  }

  // provide review
  NaluEnv::self().naluOutputP0() << "timeBlending_:          " << timeBlending_ << std::endl;
  NaluEnv::self().naluOutputP0() << "totalOuter_:            " << totalOuter_ << std::endl;
  NaluEnv::self().naluOutputP0() << "heightOuter_:           " << heightOuter_ << std::endl;
  NaluEnv::self().naluOutputP0() << "diamOuter_:             " << diamOuter_ << std::endl;
  NaluEnv::self().naluOutputP0() << "sigmaNormOuter_:        " << sigmaNormOuter_ << std::endl;
  NaluEnv::self().naluOutputP0() << "sigmaTanOuter_:         " << sigmaTanOuter_ << std::endl;
  NaluEnv::self().naluOutputP0() << "declinationOuter_:      " << declinationOuter_ << std::endl;
  NaluEnv::self().naluOutputP0() << "totalInner_:            " << totalInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "heightInner_:           " << heightInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "diamInner_:             " << diamInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "sigmaNormInner_:        " << sigmaNormInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "sigmaTanInner_:         " << sigmaTanInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "declinationInner_:      " << declinationInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "totalInnerInner_:       " << totalInnerInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "heightInnerInner_:      " << heightInnerInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "diamInnerInner_:        " << diamInnerInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "sigmaNormInnerInner_:   " << sigmaNormInnerInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "sigmaTanInnerInner_:    " << sigmaTanInnerInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "declinationInnerInner_: " << declinationInnerInner_ << std::endl;
  NaluEnv::self().naluOutputP0() << "icElevation_:           " << icElevation_ << std::endl;
  NaluEnv::self().naluOutputP0() << "icWidth_:               " << icWidth_ << std::endl;
  NaluEnv::self().naluOutputP0() << "icTrans_:               " << icTrans_ << std::endl;
  NaluEnv::self().naluOutputP0() << "originX_:               " << originX_ << std::endl;
  NaluEnv::self().naluOutputP0() << "originY_:               " << originY_ << std::endl;
  NaluEnv::self().naluOutputP0() << "timeOffset_:            " << timeOffset_ << std::endl;
}

CircularDisplacementAuxFunction::~CircularDisplacementAuxFunction()
{
  // delete each circular structure helper
  for ( size_t k = 0; k < cHelpVec_.size(); ++k )
    delete cHelpVec_[k];
}

void
CircularDisplacementAuxFunction::do_evaluate(
  const double *coords,
  const double time,
  const unsigned spatialDimension,
  const unsigned numPoints,
  double * fieldPtr,
  const unsigned fieldSize,
  const unsigned /*beginPos*/,
  const unsigned /*endPos*/) const
{
  const double fac = std::min((time-timeOffset_)/timeBlending_, 1.0);
  TanhFunction<double> tanhFunction(icTrans_, icWidth_);
      
  for(unsigned p=0; p < numPoints; ++p) {
    
    const double cX = coords[0];
    const double cY = coords[1];

    // initial elevation; use a tanh?
    const double R = std::sqrt((cX-originX_)*(cX-originX_) + (cY-originY_)*(cY-originY_));
    const double icBlend = 1.0-tanhFunction.execute(R);

    // loop over all circular structures defined and increment z; increment has scaling (+/-)
    double zSum = icElevation_*icBlend;
    for ( size_t k = 0; k < cHelpVec_.size(); ++k ) {
      zSum += cHelpVec_[k]->increment(cX, cY);
    }
    
    fieldPtr[2] = zSum*fac;
    fieldPtr += fieldSize;
    coords += spatialDimension;
  }
}

CircularHelper::CircularHelper(
  const double height, const double sigmaNorm, const double sigmaTan,
  const double normXhat, const double normYhat, const double tanXhat, const double tanYhat,
  const double centerX, const double centerY, const double scaling) :
  height_(height),
  sigmaSqNorm_(sigmaNorm*sigmaNorm),
  sigmaSqTan_(sigmaTan*sigmaTan),
  normXhat_(normXhat),
  normYhat_(normYhat),
  tanXhat_(tanXhat),
  tanYhat_(tanYhat),
  centerX_(centerX),
  centerY_(centerY),
  scaling_(scaling)
{
  // nothing else
}
  
double
CircularHelper::increment(const double cX, const double cY) 
{
  const double pNorm = normXhat_*(cX - centerX_) + normYhat_*(cY - centerY_);
  const double pTan = tanXhat_*(cX - centerX_) + tanYhat_*(cY - centerY_);
  const double dz = scaling_*height_*std::exp(-(pNorm*pNorm/(2.0*sigmaSqNorm_) + pTan*pTan/(2.0*sigmaSqTan_)));
  
  return dz;
}
  
} // namespace nalu
} // namespace Sierra
