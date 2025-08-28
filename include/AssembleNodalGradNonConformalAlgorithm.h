/*------------------------------------------------------------------------*/
/*  Copyright 2025 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef AssembleNodalGradNonConformalAlgorithm_h
#define AssembleNodalGradNonConformalAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class AssembleNodalGradNonConformalAlgorithm : public Algorithm
{
public:

  AssembleNodalGradNonConformalAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    ScalarFieldType *scalarQ,
    VectorFieldType *dqdx);

  ~AssembleNodalGradNonConformalAlgorithm();

  void execute();

  ScalarFieldType *scalarQ_;
  VectorFieldType *dqdx_;
  
  ScalarFieldType *dualNodalVolume_;
  GenericFieldType *exposedAreaVec_;

  std::vector< const stk::mesh::FieldBase *> ghostFieldVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
