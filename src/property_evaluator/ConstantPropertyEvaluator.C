/*------------------------------------------------------------------------*/
/*  Copyright 2025 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#include <property_evaluator/ConstantPropertyEvaluator.h>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// ConstantPropertyEvaluator
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
ConstantPropertyEvaluator::ConstantPropertyEvaluator(const double & value)
  : value_(value)
{
  // nothing else
}

//--------------------------------------------------------------------------
//-------- destructor ------------------------------------------------------
//--------------------------------------------------------------------------
ConstantPropertyEvaluator::~ConstantPropertyEvaluator()
{
  // nothing
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
double
ConstantPropertyEvaluator::execute(
  double *indVarList,
  stk::mesh::Entity /*node*/)
{
  return value_;
}

} // namespace nalu
} // namespace Sierra


