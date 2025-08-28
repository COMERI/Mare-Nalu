/*------------------------------------------------------------------------*/
/*  Copyright 2025 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "kernel/VolumeOfFluidScsAdvElemKernel.h"
#include "AlgTraits.h"
#include "master_element/MasterElement.h"
#include "SolutionOptions.h"

// template and scratch space
#include "BuildTemplates.h"
#include "ScratchViews.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra {
namespace nalu {

template<typename AlgTraits>
VolumeOfFluidScsAdvElemKernel<AlgTraits>::VolumeOfFluidScsAdvElemKernel(
  const stk::mesh::BulkData& bulkData,
  const SolutionOptions& solnOpts,
  ScalarFieldType* vof,
  ElemDataRequests& dataPreReqs)
  : Kernel(),
    lrscv_(sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_)->adjacentNodes())
{
  // save off fields
  const stk::mesh::MetaData& metaData = bulkData.mesh_meta_data();
  
  vofNp1_ = &(vof->field_of_state(stk::mesh::StateNP1));
  volumeFlowRate_ = metaData.get_field<double>(stk::topology::ELEMENT_RANK, "volume_flow_rate_scs");

  MasterElement *meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(AlgTraits::topo_);

  // compute shape function
  get_scs_shape_fn_data<AlgTraits>([&](double* ptr){meSCS->shape_fcn(ptr);}, v_shape_function_);

  // add master elements
  dataPreReqs.add_cvfem_surface_me(meSCS);

  // fields and data
  dataPreReqs.add_gathered_nodal_field(*vofNp1_, 1);
  dataPreReqs.add_element_field(*volumeFlowRate_, AlgTraits::numScsIp_);
}

template<typename AlgTraits>
VolumeOfFluidScsAdvElemKernel<AlgTraits>::~VolumeOfFluidScsAdvElemKernel()
{}

template<typename AlgTraits>
void
VolumeOfFluidScsAdvElemKernel<AlgTraits>::execute(
  SharedMemView<DoubleType **>& lhs,
  SharedMemView<DoubleType *>&rhs,
  ScratchViews<DoubleType>& scratchViews)
{
  SharedMemView<DoubleType*>& v_vofNp1 = scratchViews.get_scratch_view_1D(*vofNp1_);
  SharedMemView<DoubleType*>& v_vdot = scratchViews.get_scratch_view_1D(*volumeFlowRate_);

  //============================================
  // SCS contribution; vdot*alpha*nj*dS
  //============================================
  for ( int ip = 0; ip < AlgTraits::numScsIp_; ++ip ) {
    // left and right nodes for this ip
    const int il = lrscv_[2*ip];
    const int ir = lrscv_[2*ip+1];

    // save off vdot
    const DoubleType vdot = v_vdot(ip);

    // advection and diffusion
    DoubleType vofAdv = 0.0;
    for ( int ic = 0; ic < AlgTraits::nodesPerElement_; ++ic ) {

      // advection
      const DoubleType lhsfacAdv = v_shape_function_(ip,ic)*vdot;
      vofAdv += lhsfacAdv*v_vofNp1(ic);

      // lhs; il then ir
      lhs(il,ic) += lhsfacAdv;
      lhs(ir,ic) -= lhsfacAdv;
    }

    // rhs; il then ir
    rhs(il) -= vofAdv;
    rhs(ir) += vofAdv;
  }
}

INSTANTIATE_KERNEL(VolumeOfFluidScsAdvElemKernel);

}  // nalu
}  // sierra
