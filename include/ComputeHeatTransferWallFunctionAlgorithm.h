/*------------------------------------------------------------------------*/
/*  Copyright 2026 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeHeatTransferWallFunctionAlgorithm_h
#define ComputeHeatTransferWallFunctionAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;

class ComputeHeatTransferWallFunctionAlgorithm : public Algorithm
{
public:

  ComputeHeatTransferWallFunctionAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    const bool &useShifted,
    const double sigmaT);
  virtual ~ComputeHeatTransferWallFunctionAlgorithm() {}
  virtual void execute();

  const bool useShifted_;
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
};

} // namespace nalu
} // namespace Sierra

#endif
