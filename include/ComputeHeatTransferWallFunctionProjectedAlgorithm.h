/*------------------------------------------------------------------------*/
/*  Copyright 2026 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeHeatTransferWallFunctionProjectedAlgorithm_h
#define ComputeHeatTransferWallFunctionProjectedAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Ghosting;
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;
class PointInfo;

class ComputeHeatTransferWallFunctionProjectedAlgorithm : public Algorithm
{
public:

  ComputeHeatTransferWallFunctionProjectedAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool &useShifted,
    const double sigmaT,
    std::map<std::string, std::vector<std::vector<PointInfo *> > > &pointInfoMap,
    stk::mesh::Ghosting *wallFunctionGhosting);
  virtual ~ComputeHeatTransferWallFunctionProjectedAlgorithm() {}
  virtual void execute();

  const bool useShifted_;
  std::map<std::string, std::vector<std::vector<PointInfo *> > > &pointInfoMap_;
  stk::mesh::Ghosting *wallFunctionGhosting_;

  const double yplusCrit_;
  const double elog_;
  const double kappa_;
  const double sigmaT_;

  ScalarFieldType *temperature_;
  ScalarFieldType *wallTemperature_;
  ScalarFieldType *density_;
  ScalarFieldType *viscosity_;
  ScalarFieldType *specificHeat_;
  ScalarFieldType *thermalCond_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *wallFrictionVelocityBip_;
  GenericFieldType *wallNormalDistanceBip_;
  // assembled
  ScalarFieldType *assembledWallArea_;
  ScalarFieldType *referenceTemperature_;
  ScalarFieldType *heatTransferCoefficient_;
  ScalarFieldType *normalHeatFlux_;

  // data structure to parallel communicate nodal data to ghosted elements
  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
