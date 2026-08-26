#include <bout/immersed_boundary.hxx>

#include <array>
#include <regex>
#include <sstream>

//IB_TODO: Note, cant currently restart code from dump files. Need to dump IB info too I believe?
#include <bout/field.hxx>
#include <bout/globals.hxx>
#include <bout/mesh.hxx>
#include <bout/output.hxx>
using bout::globals::mesh;

ImmersedBoundary::ImmersedBoundary() {
  AUTO_TRACE();

  ASSERT0(mesh->get(bndry_mask,        "in_mask") == 0);
  ASSERT0(mesh->get(ghost_ids,        "ghost_id") == 0);
  ASSERT0(mesh->get(num_ghosts,             "ng") == 0);
  ASSERT0(mesh->get(image_inds,     "image_inds") == 0);
  ASSERT0(mesh->get(image_points,    "image_pts") == 0);
  ASSERT0(mesh->get(bndry_points,    "bndry_pts") == 0);
  ASSERT0(mesh->get(normals,           "normals") == 0);
  ASSERT0(mesh->get(norm_dist,       "norm_dist") == 0);
  ASSERT0(mesh->get(num_weights,            "nw") == 0);
  ASSERT0(mesh->get(ghost_use_bc, "ghost_use_bc") == 0);
  ASSERT0(mesh->get(ghost_use_gs, "ghost_use_gs") == 0);
  ASSERT0(mesh->get(weights_dir,   "weights_dir") == 0);
  ASSERT0(mesh->get(weights_neu,   "weights_neu") == 0);
  ASSERT0(mesh->get(R3,                      "R") == 0);
  ASSERT0(mesh->get(Z3,                      "Z") == 0);
  ASSERT0(mesh->get(bound_ids,        "bound_id") == 0);
  ASSERT0(mesh->get(bound_counts,  "bound_count") == 0);
  ASSERT0(mesh->get(num_bounds,             "nb") == 0);
  ASSERT0(mesh->get(mid_pts,           "mid_pts") == 0);
  ASSERT0(mesh->get(bnorms,             "bnorms") == 0);
  ASSERT0(mesh->get(bd_area,           "bd_area") == 0);
  ASSERT0(mesh->get(s1,                     "s1") == 0);
  ASSERT0(mesh->get(s2,                     "s2") == 0);
  ASSERT0(mesh->get(bweights,         "bweights") == 0);
  ASSERT0(mesh->get(bbase_inds,      "base_inds") == 0);
  ASSERT0(mesh->get(x_face_frac,    "face_fac_x") == 0);
  ASSERT0(mesh->get(z_face_frac,    "face_fac_z") == 0);
  ASSERT0(mesh->get(fx_grad_offset, "fx_grad_offset") == 0);
  ASSERT0(mesh->get(fz_grad_offset, "fz_grad_offset") == 0);

  if (num_weights <= 0 or num_ghosts <= 0 or num_bounds <= 0) {
    throw BoutException("Invalid number of ghost cells or weights or cut cells.");
  }

  LoadGhostPolyStencil();

  //Set up new boundary and plasma regions.
  Region<Ind3D>::RegionIndices xbndry_indices;
  Region<Ind3D>::RegionIndices zbndry_indices;
  Region<Ind3D>::RegionIndices nobndry_indices;
  Region<Ind3D>::RegionIndices nobndrycut_indices;
  Region<Ind3D>::RegionIndices bndry_indices;
  Region<Ind3D>::RegionIndices gst_indices;
  BOUT_FOR(i, bndry_mask.getRegion("RGN_ALL")) {
    //First handle guard cells correctly.
    if (i.y() < mesh->ystart   || i.y() > mesh->yend) {continue;}
    //Note, need to handle xstart-1 guard further below as special case for dagp xbndry.
    if (i.x() < mesh->xstart-1 || i.x() > mesh->xend) {
      //Add on all boundary ghost values which are ok for now. IB_TODO: Clean up logic to just necessary ghosts for dagp operators at the moment?
      if (IsGhost(i)) {
        //Ignore any ghost which can't interpolate for its image point due to lack of information.
        const auto gid = get_as<int>(ghost_ids[i]);
        const auto img_indx = GetGridInd(image_inds(gid, 0));
        if (IsBadInterpForMPI(img_indx)) {continue;}
        gst_indices.push_back(i);
      }
      continue;
    }

    if (IsInside(i)) {
      nobndry_indices.push_back(i);
      if (IsCutCell(i)) {nobndrycut_indices.push_back(i);}
    }
    else {
      if (IsGhost(i)) {
        //Need to handle MPI left boundary carefully. Sometimes unneeded ghost cells have far away images so ignore them.
        const auto gid = get_as<int>(ghost_ids[i]);
        const auto img_indx = GetGridInd(image_inds(gid, 0));
        if (i.x() == mesh->xstart-1 && IsBadInterpForMPI(img_indx)) {continue;}
        gst_indices.push_back(i);
        //Note, dont if/else here because need both indices in both regions sometimes.
        if (IsInside(i.xp())) {xbndry_indices.push_back(i);}
        if (IsInside(i.zp())) {zbndry_indices.push_back(i);}
      }
      bndry_indices.push_back(i);
    }
  }

  //IB_TODO: Only need upsert to replace RGN_NOBNDRY if considered worthwhile.
  //FV operators work on +x/+z faces only so subset of ghosts here.
  mesh->upsertRegion3D("RGN_dagp_fv_xbndry",   Region<Ind3D>(xbndry_indices));
  mesh->upsertRegion3D("RGN_dagp_fv_zbndry",   Region<Ind3D>(zbndry_indices));
  //Ghost indices.
  mesh->upsertRegion3D("RGN_IMM_BNDRY_GST",    Region<Ind3D>(gst_indices));
  //Plasma indices.
  mesh->upsertRegion3D("RGN_NO_IMM_BNDRY",     Region<Ind3D>(nobndry_indices));
  mesh->upsertRegion3D("RGN_NO_IMM_BNDRY_CUT", Region<Ind3D>(nobndrycut_indices));
  //Non-plasma indices (including ghosts).
  mesh->upsertRegion3D("RGN_IMM_BNDRY",        Region<Ind3D>(bndry_indices));

  //Check for mpi out of bounds issues with images for ghost points in given processor idx.
  const auto proc_idx = mesh->getXProcIndex();
  BOUT_FOR(i, mesh->getRegion("RGN_IMM_BNDRY_GST")) {
    const auto gid = get_as<int>(ghost_ids[i]);
    CheckInterpOkWithMPI(GetGridInd(image_inds(gid, 0)), gid, proc_idx,
                      "Image point is out of MPI bounds for ghost id ");
  }

  //Check for mpi out of bounds issues with interpolating at normal points perp to cut-cell boundary in given processor.
  BOUT_FOR(i, mesh->getRegion("RGN_NO_IMM_BNDRY_CUT")) {
    const int first_bid = get_as<int>(bound_ids[i]);
    const int bound_count = get_as<int>(bound_counts[i]);
    if (first_bid < 0 || bound_count <= 0 || first_bid + bound_count > num_bounds) {
      throw BoutException("Invalid cut-cell boundary metadata.");
    }
    for (int n = 0; n < bound_count; ++n) {
      const int bid = first_bid + n;
      CheckInterpOkWithMPI(GetGridInd(bbase_inds(bid, 0)), bid, proc_idx,
                        "Boundary normal point is out of MPI bounds for cut cell ");
      CheckInterpOkWithMPI(GetGridInd(bbase_inds(bid, 2)), bid, proc_idx,
                        "Boundary normal point is out of MPI bounds for cut cell ");
    }
  }
}

void ImmersedBoundary::LoadGhostPolyStencil() {
  use_ghost_poly =
      Options::root()["mesh"]["ib_use_ghost_poly"]
          .doc("Use precomputed higher-order immersed-boundary ghost reconstruction "
               "from ghost_poly_* grid metadata when available.")
          .withDefault<bool>(false);
  if (use_ghost_poly) {
      ASSERT0(mesh->get(ghost_poly_valid_dir,     "ghost_poly_valid_dir") == 0);
      ASSERT0(mesh->get(ghost_poly_valid_neu,     "ghost_poly_valid_neu") == 0);
      ASSERT0(mesh->get(ghost_poly_stencil_count, "ghost_poly_stencil_count") == 0);
      ASSERT0(mesh->get(ghost_poly_stencil_i,     "ghost_poly_stencil_i") == 0);
      ASSERT0(mesh->get(ghost_poly_stencil_z,     "ghost_poly_stencil_z") == 0);
      ASSERT0(mesh->get(ghost_poly_weights_dir,   "ghost_poly_weights_dir") == 0);
      ASSERT0(mesh->get(ghost_poly_weights_neu,   "ghost_poly_weights_neu") == 0);
      ASSERT0(mesh->get(ghost_poly_bc_weight_dir, "ghost_poly_bc_weight_dir") == 0);
      ASSERT0(mesh->get(ghost_poly_bc_weight_neu, "ghost_poly_bc_weight_neu") == 0);
  }
}

bool ImmersedBoundary::IsGhost(const Ind3D& ind) const {
  return get_as<int>(ghost_ids[ind]) >= 0;
}

bool ImmersedBoundary::IsUsableGhost(const Ind3D& ind) const {
  if (!IsGhost(ind)) {
    return false;
  }

  //IB_TODO: Need to check for now because might be a bad ghost (at extremes)
  //which was ignored in loops above for setting regions.
  //Though MPI GS loop should fix that if all data shared.
  const auto gid = get_as<int>(ghost_ids[ind]);
  return !IsBadInterpForMPI(GetGridInd(image_inds(gid, 0)));
}

bool ImmersedBoundary::IsInside(const Ind3D& ind) const {
  return get_as<bool>(bndry_mask[ind]);
}

bool ImmersedBoundary::IsCutCell(const Ind3D& ind) const {
  return get_as<int>(bound_ids[ind]) >= 0;
}

float ImmersedBoundary::xFaceFrac(const Ind3D& ind) const {
  return x_face_frac[ind];
}

float ImmersedBoundary::zFaceFrac(const Ind3D& ind) const {
  return z_face_frac[ind];
}

int ImmersedBoundary::xFaceGradOffset(const Ind3D& ind) const {
  return get_as<int>(fx_grad_offset[ind]);
}

int ImmersedBoundary::zFaceGradOffset(const Ind3D& ind) const {
  return get_as<int>(fz_grad_offset[ind]);
}

bool ImmersedBoundary::IsBadInterpForMPI(const int global_indx) const {
  const auto indx = mesh->getLocalXIndex(global_indx);
  const bool interpBad = indx < 0 || indx + 1 >= mesh->LocalNx;
  return interpBad;
}

void ImmersedBoundary::CheckInterpOkWithMPI(const int global_indx, const int cell_id,
                            const int proc_idx, const std::string& description) const {
    //Convert global to local index for point of interest and check in bounds of processor.
    //Note, only need to check x for now since z not parallelized.
    if (IsBadInterpForMPI(global_indx)) {
      std::stringstream strstm;
      strstm << description << cell_id << " on proc " << proc_idx
        << "." << " Increase # of guard cells in x for easy fix. " << std::endl;
      throw BoutException(strstm.str());
    }
}

ImmersedBoundary::BC_Info ImmersedBoundary::ReadBC(const std::string& bc_info) const {
  //IB_TODO, the way this is done should probably pull from the boundary factory code and be shared.
  //Match input BC expression. Allows for neumann/dirichlet, scientific notation and various brackets.
  static const std::regex bc_value_regex(
    R"(\s*(neumann|dirichlet)\s*[\(\[\{\<]\s*([-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)\s*[\)\]\}\>]\s*)"
  );
  static const std::unordered_map<std::string, BoundCond> str_to_bc {
    { "neumann",  BoundCond::NEUMANN },
    { "dirichlet", BoundCond::DIRICHLET }
  };

  std::smatch match;
  if (std::regex_match(bc_info, match, bc_value_regex)) {
    const std::string bc_type = match[1];
    const BoutReal value = std::stod(match[2]);

    auto it = str_to_bc.find(bc_type);
    if (it == str_to_bc.end()) {
      throw BoutException("Unknown BC type for immersed boundary: '" + bc_type + "'");
    }

    const auto& bc_info = BC_Info(it->second, value);
    return bc_info;
  } else {
    throw BoutException("Unknown BC input for immersed boundary: '" + bc_info + "'");
  }
}

//IB_TODO: This can take a field and check its name rather than taking a name...
//IB_TODO: Also this can probably be merged with FieldSetup to just check this always first...and continue if not needed.
bool ImmersedBoundary::CheckFieldSetUp(const std::string& name) {
    //Returns true if the field was already seen; inserts it if not.
    //IB_TODO: When calling, should throw error if fine when expected to not exist? Shouldnt be calling all over the place, just once for a field really...
    auto [_, inserted] = bc_map.try_emplace(name);
    return !inserted;
 }

void ImmersedBoundary::FieldSetup(Field3D& f) {
  const auto& bc_type = Options::root()[f.name][bc_key]
        .doc("Boundary condition to use at immersed boundary wall.")
        .withDefault("neumann(0.0)"); //Default to no flux conditions.

  //Default region to plasma cells.
  //IB_TODO: Don't set region because sort of complex at the moment...
  //If field3d operator need to include ghosts. If derivative operator, just plasma cells...
  //For now explicitly use plasma region where necessary, else use old RGN_NOBNDRY...
  //f.setRegion("RGN_NO_IMM_BNDRY");
  //ddt(f).setRegion("RGN_NO_IMM_BNDRY");

  //Load BC info and clear out non-plasma cells.
  const auto& bc_info = ReadBC(bc_type);
  bc_map[f.name] = bc_info;
  BOUT_FOR(i, f.getRegion("RGN_IMM_BNDRY")) {
    f[i] = 0.0;
    ddt(f)[i] = 0.0;
  }
}

void ImmersedBoundary::FloorField(Field3D& f, const float val) const {
  //Just floor plasma cells, not ghosts.
  //IB_TODO: Reset boundary after flooring? Do outside afterwards on case-by-case basis currently.
  BOUT_FOR(i, f.getRegion("RGN_NO_IMM_BNDRY")) {
    if (f[i] < val) {
      f[i] = val;
    }
  }
}

//IB_TODO: Interpolate dfdn from dfdn at segment endpoints (bdy/cell intersections).
//IB_TODO: Also need fb and a at boundary approximation.
void ImmersedBoundary::ComputeBoundaryFluxes(const Field3D& a, const Field3D& f, Field3D& result) const {
  //Get bc info from field first.
  const auto bc_type = bc_map.at(f.name).first;
  const auto bc_val  = bc_map.at(f.name).second;

  //Loop over cut-cells w/ centers in plasma and calculate flux through wall boundary.
  BOUT_FOR(i, f.getRegion("RGN_NO_IMM_BNDRY_CUT")) {
    const int first_bid = get_as<int>(bound_ids[i]);
    const int bound_count = get_as<int>(bound_counts[i]);
    if (first_bid < 0 || bound_count <= 0 || first_bid + bound_count > num_bounds) {
      throw BoutException("Invalid cut-cell boundary metadata.");
    }

    for (int n = 0; n < bound_count; ++n) {
      const int b = first_bid + n;

      BoutReal dfdn = 0.0;
      if (bc_type == BoundCond::NEUMANN) {
        dfdn = bc_val; //Note, using outward normal derivative here.
      } else if (bc_type == BoundCond::DIRICHLET) {
        //Perform bilinear interpolation for points along the boundary normal, and 
        //calculate normal derivative at boundary.
        const int ai0_global = GetGridInd(bbase_inds(b, 0));
        const int aj0 = GetGridInd(bbase_inds(b, 1));
        const int bi0_global = GetGridInd(bbase_inds(b, 2));
        const int bj0 = GetGridInd(bbase_inds(b, 3));

        const int ai0 = mesh->getLocalXIndex(ai0_global);
        const int bi0 = mesh->getLocalXIndex(bi0_global);

        const BoutReal wA00 = bweights(b, 0);
        const BoutReal wA01 = bweights(b, 1);
        const BoutReal wA10 = bweights(b, 2);
        const BoutReal wA11 = bweights(b, 3);

        const BoutReal wB00 = bweights(b, 4);
        const BoutReal wB01 = bweights(b, 5);
        const BoutReal wB10 = bweights(b, 6);
        const BoutReal wB11 = bweights(b, 7);

        // Bilinear interpolation using the correct local stencil
        const BoutReal fA = wA00*f(ai0,1,aj0) + wA01*f(ai0,1,aj0+1) + wA10*f(ai0+1,1,aj0) + wA11*f(ai0+1,1,aj0+1);
        const BoutReal fB = wB00*f(bi0,1,bj0) + wB01*f(bi0,1,bj0+1) + wB10*f(bi0+1,1,bj0) + wB11*f(bi0+1,1,bj0+1);

        // Boundary value and one-sided derivative from J&C '98.
        const BoutReal fb = bc_val;

        const BoutReal term1 = (s2/s1) * (fb - fA);
        const BoutReal term2 = (s1/s2) * (fb - fB);

        dfdn = (term1 - term2) / (s2 - s1);
      } else {
         throw BoutException(bc_exception);
      }

      //Flux assumes normal is outward for dfdn.
      const auto bdy_flux = -a[i] * bd_area[b] * dfdn;
      result[i] -= bdy_flux;
    }
  }
}

// Calculate image value from nearby grid points. Note weights are defined as
// w00, w01, w10, w11 with indices (x,z). All other values follow from there.f
BoutReal ImmersedBoundary::GetImageValue(Field3D& f, const int gid,
            const BoutReal bc_val, const BoundCond bc_type) const {
  // Get nearby vals to image from floating point index.
  // Only x mpi-parallelized at the moment.
  const int indx_global = GetGridInd(image_inds(gid,0));
  const int indx = mesh->getLocalXIndex(indx_global);
  const int indz = GetGridInd(image_inds(gid,1));
  
  //IB_TODO: Need indy for ghost cells in y? Same below. Or do everything in 2d?
  //IB_TODO: Use num_weights instead of 4, same for below where 4s used.
  auto node_vals = std::array<BoutReal, 4>{f(indx,1,indz), f(indx,1,indz+1), 
                                           f(indx+1,1,indz), f(indx+1,1,indz+1)};

  if (!(bc_type == BoundCond::DIRICHLET || bc_type == BoundCond::NEUMANN)) {
    throw BoutException(bc_exception);
  }
  const auto& weights = bc_type == BoundCond::DIRICHLET ? weights_dir : weights_neu;

  for (size_t i = 0; i < num_weights; ++i) {
    if (get_as<bool>(ghost_use_bc(gid, i))) {
      node_vals[i] = bc_val;
    }
  }

  BoutReal image_val = 0.0;
  for (size_t i = 0; i < num_weights; ++i) {
    image_val += weights(gid,i)*node_vals[i];
  }

  return image_val;
}

BoutReal ImmersedBoundary::GetGhostValue(const BoutReal image_val, const int gid,
                          const BoutReal bc_val, const BoundCond bc_type) const {
  switch (bc_type) {
    case BoundCond::DIRICHLET:
      return 2.0*bc_val - image_val;
    case BoundCond::NEUMANN:
      //Note, use + here for outward normal convention.
      return image_val + 2.0*norm_dist[gid]*bc_val;
    default:
        throw BoutException(bc_exception);
  }
}

void ImmersedBoundary::SetBoundary(Field3D& f) {
  const auto bc_type = bc_map.at(f.name).first;
  const auto bc_val  = bc_map.at(f.name).second;

  //Multiple ghost node case requires Gauss-Seidel convergence due to coupled initial guesses.
  constexpr int max_gs_iters = 12; //IB_TODO: Make a vector to store old values, reset to 0 at start, and compare to stop sweeping.
  for (int it = 0; it < max_gs_iters; ++it) {
    //NOTE: This is "Gauss–Seidel-like" because we update f in-place as we sweep.
    BOUT_FOR(i, f.getRegion("RGN_IMM_BNDRY_GST")) {
      const auto gid = get_as<int>(ghost_ids[i]);

      if (use_ghost_poly && GhostPolyValid(gid, bc_type) && IsPolynomialGhostStencilLocal(gid)) {
        if (it == 0) {
          f[i] = GetPolynomialGhostValue(f, gid, i.y(), bc_val, bc_type);
        }
        continue;
      }

      //If first iteration or more than 1 ghost (GS iteration to converge).
      //IB_TODO: Make GS iteration MPI compliant.
      if (it == 0 || get_as<bool>(ghost_use_gs[gid])) {
        const auto image_val = GetImageValue(f, gid, bc_val, bc_type);
        const auto ghost_val = GetGhostValue(image_val, gid, bc_val, bc_type);
        f[i] = ghost_val;
      }
    }
  }

  //IB_TODO: Update to only MYG guard cells when no longer axisymmetric and ny != 1.
  BOUT_FOR(i, f.getRegion("RGN_IMM_BNDRY_GST")) {
    //Loop in case MYG > 1.
    for (int yoffset = 1; yoffset <= f.getMesh()->ystart; ++yoffset) {
      f[i.yp(yoffset)] = f[i];
      f[i.ym(yoffset)] = f[i];
    }
  }
}


//Experimental higher order polynomial ghost stencil functions.
//All experimental code GPT-generated and full of overkill checks...Needs major clean up.
bool ImmersedBoundary::GhostPolyValid(const int gid,
                          const BoundCond bc_type) const {
  switch (bc_type) {
    case BoundCond::DIRICHLET:
      return get_as<bool>(ghost_poly_valid_dir[gid]);
    case BoundCond::NEUMANN:
      return get_as<bool>(ghost_poly_valid_neu[gid]);
    default:
      throw BoutException(bc_exception);
  }
}

bool ImmersedBoundary::IsPolynomialGhostStencilLocal(const int gid) const {
  //Many of these checks seem redundant and can be handled in python/grid generation.
  const auto [ngi, nstencil] = ghost_poly_stencil_i.shape();

  const int stencil_count = get_as<int>(ghost_poly_stencil_count[gid]);
  if (stencil_count <= 0 || stencil_count > nstencil) {
    return false;
  }

  for (int n = 0; n < stencil_count; ++n) {
    const int global_x = get_as<int>(ghost_poly_stencil_i(gid, n));
    const int local_x = mesh->getLocalXIndex(global_x);
    const int z = get_as<int>(ghost_poly_stencil_z(gid, n));
    if (local_x < 0 || local_x >= mesh->LocalNx || z < 0 || z >= mesh->LocalNz) {
      return false;
    }
  }

  return true;
}

BoutReal ImmersedBoundary::GetPolynomialGhostValue(Field3D& f, const int gid,
                          const int y, const BoutReal bc_val,
                          const BoundCond bc_type) const {
  const auto& weights =
      bc_type == BoundCond::DIRICHLET ? ghost_poly_weights_dir : ghost_poly_weights_neu;
  const auto& bc_weight =
      bc_type == BoundCond::DIRICHLET ? ghost_poly_bc_weight_dir : ghost_poly_bc_weight_neu;
  const auto [ng, nstencil] = weights.shape();
  if (gid >= ng) {
    throw BoutException("Invalid immersed-boundary polynomial ghost id.");
  }

  const int stencil_count = get_as<int>(ghost_poly_stencil_count[gid]);
  BoutReal ghost_val = bc_weight[gid] * bc_val;
  for (int n = 0; n < stencil_count; ++n) {
    const int global_x = get_as<int>(ghost_poly_stencil_i(gid, n));
    const int local_x = mesh->getLocalXIndex(global_x);
    const int z = get_as<int>(ghost_poly_stencil_z(gid, n));
    ghost_val += weights(gid, n) * f(local_x, y, z);
  }

  return ghost_val;
}