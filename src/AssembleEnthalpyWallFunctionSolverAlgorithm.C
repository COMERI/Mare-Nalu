/*------------------------------------------------------------------------*/
/*  Copyright 2026 COMERI.                                                */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Mare-Nalu     */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


// nalu
#include <AssembleEnthalpyWallFunctionSolverAlgorithm.h>
#include <SolverAlgorithm.h>
#include <EquationSystem.h>
#include <LinearSystem.h>
#include <FieldTypeDef.h>
#include <Realm.h>
#include <master_element/MasterElement.h>

// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

// basic c++
#include <cmath>

namespace sierra{
namespace nalu{

//==========================================================================
// Class Definition
//==========================================================================
// AssembleEnthalpyWallFunctionSolverAlgorithm - elem wall function
//==========================================================================
//--------------------------------------------------------------------------
//-------- constructor -----------------------------------------------------
//--------------------------------------------------------------------------
AssembleEnthalpyWallFunctionSolverAlgorithm::AssembleEnthalpyWallFunctionSolverAlgorithm(
  Realm &realm,
  stk::mesh::Part *part,
  EquationSystem *eqSystem,
  const bool &useShifted,
  const double sigmaT)
  : SolverAlgorithm(realm, part, eqSystem),
    useShifted_(useShifted),
    yplusCrit_(11.63),
    elog_(9.8),
    kappa_(realm.get_turb_model_constant(TM_kappa)),
    sigmaT_(sigmaT)
{
  // save off fields
  stk::mesh::MetaData & meta_data = realm_.meta_data();
  temperature_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "temperature");
  wallTemperature_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "temperature_bc");
  density_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "density");
  viscosity_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "viscosity");
  specificHeat_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "specific_heat");
  thermalCond_ = meta_data.get_field<double>(stk::topology::NODE_RANK, "thermal_conductivity");
  exposedAreaVec_ = meta_data.get_field<double>(meta_data.side_rank(), "exposed_area_vector");
  wallFrictionVelocityBip_ = meta_data.get_field<double>(meta_data.side_rank(), "wall_friction_velocity_bip");
  wallNormalDistanceBip_ = meta_data.get_field<double>(meta_data.side_rank(), "wall_normal_distance_bip");
}

//--------------------------------------------------------------------------
//-------- initialize_connectivity -----------------------------------------
//--------------------------------------------------------------------------
void
AssembleEnthalpyWallFunctionSolverAlgorithm::initialize_connectivity()
{
  eqSystem_->linsys_->buildFaceToNodeGraph(partVec_);
}

//--------------------------------------------------------------------------
//-------- execute ---------------------------------------------------------
//--------------------------------------------------------------------------
void
AssembleEnthalpyWallFunctionSolverAlgorithm::execute()
{

  stk::mesh::BulkData & bulk_data = realm_.bulk_data();
  stk::mesh::MetaData & meta_data = realm_.meta_data();

  const int nDim = meta_data.spatial_dimension();

  // space for LHS/RHS; nodesPerFace*nodesPerFace and nodesPerFace
  std::vector<double> lhs;
  std::vector<double> rhs;
  std::vector<int> scratchIds;
  std::vector<double> scratchVals;
  std::vector<stk::mesh::Entity> connected_nodes;

  // nodal fields to gather
  std::vector<double> ws_temperature;
  std::vector<double> ws_wall_temperature;
  std::vector<double> ws_density;
  std::vector<double> ws_viscosity;
  std::vector<double> ws_specific_heat;
  std::vector<double> ws_thermal_cond;

  // master element
  std::vector<double> ws_shape_function;

  // deal with state
  ScalarFieldType &densityNp1 = density_->field_of_state(stk::mesh::StateNP1);

  // define some common selectors
  stk::mesh::Selector s_locally_owned_union = meta_data.locally_owned_part()
    &stk::mesh::selectUnion(partVec_);

  stk::mesh::BucketVector const& face_buckets =
    realm_.get_buckets( meta_data.side_rank(), s_locally_owned_union );
  for ( stk::mesh::BucketVector::const_iterator ib = face_buckets.begin();
        ib != face_buckets.end() ; ++ib ) {
    stk::mesh::Bucket & b = **ib ;

    // face master element
    MasterElement *meFC = sierra::nalu::MasterElementRepo::get_surface_master_element(b.topology());
    const int nodesPerFace = meFC->nodesPerElement_;
    const int numScsBip = meFC->numIntPoints_;

    // mapping from ip to nodes for this ordinal; face perspective (use with face_node_relations)
    const int *ipNodeMap = meFC->ipNodeMap();

    // resize some things; matrix related
    const int lhsSize = nodesPerFace*nodesPerFace;
    const int rhsSize = nodesPerFace;
    lhs.resize(lhsSize);
    rhs.resize(rhsSize);
    scratchIds.resize(rhsSize);
    scratchVals.resize(rhsSize);
    connected_nodes.resize(nodesPerFace);

    // algorithm related; element
    ws_density.resize(nodesPerFace);
    ws_viscosity.resize(nodesPerFace);
    ws_specific_heat.resize(nodesPerFace);
    ws_shape_function.resize(numScsBip*nodesPerFace);

    // pointers
    double *p_lhs = &lhs[0];
    double *p_rhs = &rhs[0];
    double *p_temperature = &ws_temperature[0];
    double *p_wall_temperature = &ws_wall_temperature[0];
    double *p_density = &ws_density[0];
    double *p_viscosity = &ws_viscosity[0];
    double *p_specific_heat = &ws_specific_heat[0];
    double *p_thermal_cond = &ws_thermal_cond[0];
    double *p_shape_function = &ws_shape_function[0];

    // shape functions
    if ( useShifted_ )
      meFC->shifted_shape_fcn(&p_shape_function[0]);
    else
      meFC->shape_fcn(&p_shape_function[0]);

    const stk::mesh::Bucket::size_type length   = b.size();

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      // zero lhs/rhs
      for ( int p = 0; p < lhsSize; ++p )
        p_lhs[p] = 0.0;
      for ( int p = 0; p < rhsSize; ++p )
        p_rhs[p] = 0.0;

      // get face
      stk::mesh::Entity face = b[k];

      //======================================
      // gather nodal data off of face
      //======================================
      stk::mesh::Entity const * face_node_rels = bulk_data.begin_nodes(face);
      for ( int ni = 0; ni < nodesPerFace; ++ni ) {
        stk::mesh::Entity node = face_node_rels[ni];
        connected_nodes[ni] = node;

        // gather scalars
        p_temperature[ni]   = *stk::mesh::field_data(*temperature_, node);
        p_wall_temperature[ni] = *stk::mesh::field_data(*wallTemperature_, node);
        p_density[ni]       = *stk::mesh::field_data(densityNp1, node);
        p_viscosity[ni]     = *stk::mesh::field_data(*viscosity_, node);
        p_specific_heat[ni] = *stk::mesh::field_data(*specificHeat_, node);
        p_thermal_cond[ni]  = *stk::mesh::field_data(*thermalCond_, node);
      }

      // pointer to face data
      const double * areaVec = stk::mesh::field_data(*exposedAreaVec_, face);
      const double *wallNormalDistanceBip = stk::mesh::field_data(*wallNormalDistanceBip_, face);
      const double *wallFrictionVelocityBip = stk::mesh::field_data(*wallFrictionVelocityBip_, face);

      // loop over face nodes
      for ( int ip = 0; ip < numScsBip; ++ip ) {

        const int ipNdim = ip*nDim;
        const int ipNpf = ip*nodesPerFace;

        const int nn = ipNodeMap[ip];

        // zero out vector quantities; squeeze in aMag
        double aMag = 0.0;
        for ( int j = 0; j < nDim; ++j ) {
          const double axj = areaVec[ipNdim+j];
          aMag += axj*axj;
        }
        aMag = std::sqrt(aMag);

        // interpolate to bip
        double tBip = 0.0;
        double walltBip = 0.0;
        double rhoBip = 0.0;
        double muBip = 0.0;
        double cpBip = 0.0;
	double kBip = 0.0;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_shape_function[ipNpf+ic];
          tBip += r*p_temperature[ic];
          walltBip += r*p_wall_temperature[ic];
          rhoBip += r*p_density[ic];
          muBip += r*p_viscosity[ic];
          cpBip += r*p_specific_heat[ic];
          kBip += r*p_thermal_cond[ic];
        }

        // extract bip data
        const double yp = wallNormalDistanceBip[ip];
        const double utau= wallFrictionVelocityBip[ip];

        // determine LOW quantities
        const double yplus = rhoBip*yp*utau/muBip;
        const double uplus = 1.0/kappa_*std::log(elog_*yplus);

	// account for the disparate thermal boundary layer
	const double lamPr = cpBip*muBip/kBip;
	const double prRatio = lamPr/sigmaT_;
	
	const double Pf = 9.24*(std::pow(prRatio,0.75)-1.0)*(1.0+0.28*std::exp(-0.007*prRatio));
	const double Tplus = sigmaT_*(uplus+Pf);
	
        double lambda = kBip/cpBip*aMag;
        if ( yplus > yplusCrit_)
          lambda = rhoBip*cpBip*utau/Tplus*aMag;

	const double hflux = lambda*(tBip-walltBip);
	p_rhs[nn] -= hflux;
	
	// sensitivities
        const int rowR = nn*nodesPerFace;
        const double lhsFac = lambda/cpBip*aMag;
        for ( int ic = 0; ic < nodesPerFace; ++ic ) {
          const double r = p_shape_function[ipNpf+ic];
          p_lhs[rowR+ic] += r*lhsFac;
        }
      }
      
      apply_coeff(connected_nodes, scratchIds, scratchVals, rhs, lhs, __FILE__);

    }
  }
}


} // namespace nalu
} // namespace Sierra
