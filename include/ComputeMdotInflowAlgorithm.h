/*------------------------------------------------------------------------*/
/*  Copyright 2025 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeMdotInflowAlgorithm_h
#define ComputeMdotInflowAlgorithm_h

#include<Algorithm.h>
#include<FieldTypeDef.h>

// stk
#include <stk_mesh/base/Part.hpp>

namespace sierra{
namespace nalu{

class Realm;

class ComputeMdotInflowAlgorithm : public Algorithm
{
public:

  ComputeMdotInflowAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    bool useShifted);
  ~ComputeMdotInflowAlgorithm();

  void execute();

  const bool useShifted_;

  VectorFieldType *velocityBC_;
  ScalarFieldType *densityBC_;
  GenericFieldType *exposedAreaVec_;
};

} // namespace nalu
} // namespace Sierra

#endif
