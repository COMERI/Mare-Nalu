/*------------------------------------------------------------------------*/
/*  Copyright 2025 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleRadTransEdgeSolverAlgorithm_h
#define AssembleRadTransEdgeSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class RadiativeTransportEquationSystem;
class Realm;

class AssembleRadTransEdgeSolverAlgorithm : public SolverAlgorithm
{
public:
  
  AssembleRadTransEdgeSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    RadiativeTransportEquationSystem *radEqSystem);
  virtual ~AssembleRadTransEdgeSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

  const RadiativeTransportEquationSystem *radEqSystem_;

  ScalarFieldType *intensity_;
  VectorFieldType *edgeAreaVec_;
  VectorFieldType *coordinates_;
  ScalarFieldType *absorption_;
  ScalarFieldType *scattering_;
  ScalarFieldType *scalarFlux_;
  ScalarFieldType *radiationSource_;
  ScalarFieldType *dualNodalVolume_;
  
};

} // namespace nalu
} // namespace Sierra

#endif
