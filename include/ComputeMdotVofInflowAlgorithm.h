/*------------------------------------------------------------------------*/
/*  Copyright 2025 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMdotVofInflowAlgorithm_h
#define ComputeMdotVofInflowAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMdotVofInflowAlgorithm : public Algorithm
{
public:

  ComputeMdotVofInflowAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    bool useShifted);
  ~ComputeMdotVofInflowAlgorithm();

  void execute();

  const bool useShifted_;

  VectorFieldType *velocityBC_;
  ScalarFieldType *densityBC_;
  GenericFieldType *exposedAreaVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
