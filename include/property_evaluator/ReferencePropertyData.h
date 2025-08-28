/*------------------------------------------------------------------------*/
/*  Copyright 2025 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ReferencePropertyData_h
#define ReferencePropertyData_h

#include <Enums.h>

#include <vector>

namespace sierra{
namespace nalu{

class ReferencePropertyData {
public:
  
  ReferencePropertyData();
  ~ReferencePropertyData();
  
  std::string speciesName_;
  double mw_;
  double massFraction_;
  double stoichiometry_;
  double primaryMassFraction_;
  double secondaryMassFraction_;

};

} // namespace nalu
} // namespace Sierra

#endif
