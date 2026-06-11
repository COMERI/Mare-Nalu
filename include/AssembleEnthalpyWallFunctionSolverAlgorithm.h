/*------------------------------------------------------------------------*/
/*  Copyright 2026 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleEnthalpyWallFunctionSolverAlgorithm_h
#define AssembleEnthalpyWallFunctionSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class Realm;

class AssembleEnthalpyWallFunctionSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleEnthalpyWallFunctionSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    const bool &useShifted,
    const double sigmaT);
  virtual ~AssembleEnthalpyWallFunctionSolverAlgorithm() {}
  virtual void initialize_connectivity();
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
};

} // namespace nalu
} // namespace Sierra

#endif
