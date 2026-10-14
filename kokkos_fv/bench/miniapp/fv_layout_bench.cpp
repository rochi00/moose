// Standalone Kokkos benchmark of linear FV assembly layouts.
//
// Reads a mesh dump written by the MOOSE test object FVMeshDumper (ElemInfo/FaceInfo geometry)
// and assembles the same steady advection-diffusion-source system as
// kokkos_fv_bench/inputs/advdiff_kokkos.i with several data layouts and assembly algorithms:
//
//   V0 directed_atomic   today's Kokkos-MOOSE scheme: per-(side, cell) padded geometry, one thread
//                        per directed side, coefficients evaluated per contribution, column found
//                        by binary search, atomic accumulation, separate matrix/RHS loops
//   V1 face_atomic       face SoA, one thread per face, coefficients evaluated once per face,
//                        binary search, atomics
//   V2 face_atomic_slots V1 with precomputed CSR slots instead of the search
//   V3 face_gather       face pass writes per-face scratch (diag_E, diag_N, rhs_E, rhs_N) and the
//                        off-diagonals into their slots; per-cell gather over elem_faces sums the
//                        diagonal and RHS (no atomics)
//   V4 halfface_gather   face pass writes half-face scratch (side, cell); per-cell gather writes
//                        the whole matrix row and RHS (no atomics)
//
// Green-Gauss cell gradient variants (--grad-variants), on a smooth analytic cell field:
//
//   G0 side_gather       per-cell gather over today's padded per-(side, cell) geometry
//   G1 face_gather       per-cell gather over the face SoA through elem_faces (no atomics)
//   G2 face_atomic       one thread per face, atomic adds into both cells, then a volume pass
//   G3 face_atomic_vol   G2 with 1/V folded into the atomic adds (single face pass)
//   G0C                  G0 with the padded side arrays cell-fastest (coalesced on GPUs)
//   G1C                  G1 with elem_faces cell-fastest
//   G1CS                 G1C with face normals and centroids as separate x, y, z arrays (SoA)
//
// V2F is V2 with the diffusion and advection face passes fused into one kernel.
// --launch-bench times kernel launches with a small functor (passed as a kernel argument on CUDA)
// and a 6 KB functor (the size of Kokkos-MOOSE linear FV functors, which CUDA launches through
// constant memory), with and without a fence after every launch.
//
// Face values follow the CPU LinearFVGradientManager: linear interpolation with FaceInfo::gC() on
// internal faces, the Dirichlet value on Dirichlet faces, the cell value elsewhere. Every
// gradient variant is checked against a serial host face loop in the original cell numbering.
//
// Cell numbering (native | morton | rcm | scrambled) and face order (native | owner) are options.
// Every variant is checked against V0 entrywise; V0 is checked against the MOOSE-assembled
// matrix when its PETSc binary dump is given. All variants use the analytic coefficients of the
// MOOSE input evaluated inline (no parsed-function interpreter).

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <iterator>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <string>
#include <vector>

using ExecSpace = Kokkos::DefaultExecutionSpace;
using MemSpace = ExecSpace::memory_space;
template <typename T>
using View1 = Kokkos::View<T *, MemSpace>;
template <typename T>
using View2L = Kokkos::View<T **, Kokkos::LayoutLeft, MemSpace>; // (side, cell): side fastest
template <typename T>
using View2R = Kokkos::View<T **, Kokkos::LayoutRight, MemSpace>; // (side, cell): cell fastest

struct R3
{
  double x, y, z;
};

// ------------------------------------------------------------------------------------------------
// Problem definition (matches inputs/advdiff_kokkos.i)
// ------------------------------------------------------------------------------------------------
namespace phys
{
constexpr double vx = 1.0, vy = 0.5, vz = 0.25;
KOKKOS_INLINE_FUNCTION double
k(double x, double, double)
{
  return 1.0 + 0.5 * x;
}
KOKKOS_INLINE_FUNCTION double
s(double x, double y, double z)
{
  return 1.0 + x * y * z;
}
KOKKOS_INLINE_FUNCTION double
g(double x, double y, double z)
{
  return 1.0 + x + y + z;
}
// GeneratedMesh 3D boundary ids: back 0, bottom 1, right 2, top 3, left 4, front 5
KOKKOS_INLINE_FUNCTION bool
isDirichlet(int bid)
{
  return bid == 0 || bid == 1 || bid == 4;
}
// Smooth cell field for the gradient variants; its Dirichlet faces use phi at the face centroid
KOKKOS_INLINE_FUNCTION double
phi(double x, double y, double z)
{
  return 1.0 + x * y + z * z + Kokkos::sin(3.0 * x);
}
}

// ------------------------------------------------------------------------------------------------
// Host mesh as dumped by FVMeshDumper
// ------------------------------------------------------------------------------------------------
struct HostMesh
{
  int64_t nc = 0, nf = 0;
  std::vector<int64_t> elem_id, dof;
  std::vector<double> vol;
  std::vector<R3> cen;
  std::vector<int64_t> fe, fn;
  std::vector<int32_t> se, sn, bid;
  std::vector<double> area, gc;
  std::vector<R3> nrm, fc, dcn;
};

template <typename T>
T
rd(std::ifstream & in)
{
  T v;
  in.read(reinterpret_cast<char *>(&v), sizeof(T));
  return v;
}

R3
rd3(std::ifstream & in)
{
  R3 p;
  p.x = rd<double>(in);
  p.y = rd<double>(in);
  p.z = rd<double>(in);
  return p;
}

HostMesh
readMesh(const std::string & path)
{
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("cannot open " + path);
  char magic[4];
  in.read(magic, 4);
  if (std::strncmp(magic, "FVD1", 4))
    throw std::runtime_error("bad magic in " + path);
  HostMesh m;
  m.nc = rd<int64_t>(in);
  m.nf = rd<int64_t>(in);
  m.elem_id.resize(m.nc);
  m.dof.resize(m.nc);
  m.vol.resize(m.nc);
  m.cen.resize(m.nc);
  for (int64_t i = 0; i < m.nc; ++i)
  {
    m.elem_id[i] = rd<int64_t>(in);
    m.dof[i] = rd<int64_t>(in);
    m.vol[i] = rd<double>(in);
    m.cen[i] = rd3(in);
  }
  m.fe.resize(m.nf);
  m.fn.resize(m.nf);
  m.se.resize(m.nf);
  m.sn.resize(m.nf);
  m.area.resize(m.nf);
  m.nrm.resize(m.nf);
  m.fc.resize(m.nf);
  m.dcn.resize(m.nf);
  m.gc.resize(m.nf);
  m.bid.resize(m.nf);
  for (int64_t f = 0; f < m.nf; ++f)
  {
    m.fe[f] = rd<int64_t>(in);
    m.fn[f] = rd<int64_t>(in);
    m.se[f] = rd<int32_t>(in);
    m.sn[f] = rd<int32_t>(in);
    m.area[f] = rd<double>(in);
    m.nrm[f] = rd3(in);
    m.fc[f] = rd3(in);
    m.dcn[f] = rd3(in);
    m.gc[f] = rd<double>(in);
    m.bid[f] = rd<int32_t>(in);
  }
  if (!in)
    throw std::runtime_error("truncated mesh file " + path);
  return m;
}

// ------------------------------------------------------------------------------------------------
// Orderings
// ------------------------------------------------------------------------------------------------
// Returns perm with perm[old_cell] = new_cell
std::vector<int>
cellOrdering(const HostMesh & m, const std::string & kind)
{
  const int n = m.nc;
  std::vector<int> perm(n);
  if (kind == "native")
  {
    std::iota(perm.begin(), perm.end(), 0);
    return perm;
  }
  std::vector<int> order(n); // order[new] = old
  std::iota(order.begin(), order.end(), 0);
  if (kind == "scrambled")
  {
    std::mt19937 rng(12345);
    std::shuffle(order.begin(), order.end(), rng);
  }
  else if (kind == "morton")
  {
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (auto & c : m.cen)
    {
      const double p[3] = {c.x, c.y, c.z};
      for (int d = 0; d < 3; ++d)
        lo[d] = std::min(lo[d], p[d]), hi[d] = std::max(hi[d], p[d]);
    }
    auto spread = [](uint64_t v)
    {
      v &= 0x1fffff;
      v = (v | v << 32) & 0x1f00000000ffffULL;
      v = (v | v << 16) & 0x1f0000ff0000ffULL;
      v = (v | v << 8) & 0x100f00f00f00f00fULL;
      v = (v | v << 4) & 0x10c30c30c30c30c3ULL;
      v = (v | v << 2) & 0x1249249249249249ULL;
      return v;
    };
    std::vector<uint64_t> key(n);
    for (int i = 0; i < n; ++i)
    {
      const double p[3] = {m.cen[i].x, m.cen[i].y, m.cen[i].z};
      uint64_t q[3];
      for (int d = 0; d < 3; ++d)
      {
        const double t = hi[d] > lo[d] ? (p[d] - lo[d]) / (hi[d] - lo[d]) : 0.0;
        q[d] = uint64_t(std::min(t, 1.0) * ((1u << 21) - 1));
      }
      key[i] = spread(q[0]) | spread(q[1]) << 1 | spread(q[2]) << 2;
    }
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return key[a] < key[b]; });
  }
  else if (kind == "rcm")
  {
    std::vector<std::vector<int>> adj(n);
    for (int64_t f = 0; f < m.nf; ++f)
      if (m.fn[f] >= 0)
      {
        adj[m.fe[f]].push_back(m.fn[f]);
        adj[m.fn[f]].push_back(m.fe[f]);
      }
    std::vector<char> seen(n, 0);
    order.clear();
    std::vector<int> by_degree(n);
    std::iota(by_degree.begin(), by_degree.end(), 0);
    std::stable_sort(by_degree.begin(),
                     by_degree.end(),
                     [&](int a, int b) { return adj[a].size() < adj[b].size(); });
    for (int start : by_degree)
    {
      if (seen[start])
        continue;
      std::queue<int> q;
      q.push(start);
      seen[start] = 1;
      while (!q.empty())
      {
        const int c = q.front();
        q.pop();
        order.push_back(c);
        auto nb = adj[c];
        std::sort(
            nb.begin(), nb.end(), [&](int a, int b) { return adj[a].size() < adj[b].size(); });
        for (int x : nb)
          if (!seen[x])
            seen[x] = 1, q.push(x);
      }
    }
    std::reverse(order.begin(), order.end());
  }
  else
    throw std::runtime_error("unknown cell ordering " + kind);
  for (int i = 0; i < n; ++i)
    perm[order[i]] = i;
  return perm;
}

// ------------------------------------------------------------------------------------------------
// Device data for all variants
// ------------------------------------------------------------------------------------------------
struct Data
{
  int nc = 0, nf = 0, n_internal = 0, max_sides = 0;
  // cells (renumbered)
  View1<double> vol;
  View1<R3> cen;
  // faces (SoA): internal faces first, then boundary faces
  View1<int> fe, fn, fse, fsn, fbid;
  View1<double> farea, fdcn_mag, fdcf_mag;
  View1<R3> fnrm, fc;
  // cell -> faces
  View2L<int> elem_faces;
  // CSR (rows = cells, sorted columns)
  View1<int> row_ptr, col;
  int nnz = 0;
  // precomputed slots
  View1<int> diag_slot, slot_EN, slot_NE;
  View2L<int> hf_slot; // (side, cell) -> CSR slot of that side's off-diagonal, -1 if boundary
  // per-(side, cell) padded "directed side" geometry for V0
  View2L<double> ds_area, ds_dcn_mag, ds_dcf_mag;
  View2L<R3> ds_nrm, ds_fc;
  View2L<int> ds_nb, ds_bid;
  View1<int> v0_int_cell, v0_int_side, v0_bnd_cell, v0_bnd_side;
  // outputs and scratch
  View1<double> val, rhs;
  View1<double> fs_dE, fs_dN, fs_rE, fs_rN; // per-face scratch (V3)
  View2L<double> hf_d, hf_o, hf_r;          // half-face scratch (V4)
  // gradients: FaceInfo::gC() per face, the side's own interpolation weight per (side, cell) for
  // G0, the cell field, and the gradient with the components of a cell together (3 * c + d)
  View1<double> fgc;
  View2L<double> ds_gc;
  View1<double> phi, grad;
  // cell-fastest copies for G0C and G1C, and SoA face geometry for G1CS
  View2R<double> dsr_area, dsr_gc;
  View2R<R3> dsr_nrm, dsr_fc;
  View2R<int> dsr_nb, dsr_bid, elem_faces_r;
  View1<double> fnx, fny, fnz, fcx, fcy, fcz;
};

// Host-only information (kept out of Data, which is captured by value in every kernel)
struct HostInfo
{
  std::vector<int> new_to_old_cell;
  size_t bytes_face_soa = 0, bytes_side_padded = 0, bytes_elem_faces = 0, bytes_slots = 0,
         bytes_face_scratch = 0, bytes_hf_scratch = 0, bytes_grad_face = 0, bytes_grad_side = 0;
};

template <typename V, typename T>
V
toDevice(const std::string & name, const std::vector<T> & h)
{
  V d(name, h.size());
  auto hm = Kokkos::create_mirror_view(d);
  for (size_t i = 0; i < h.size(); ++i)
    hm(i) = h[i];
  Kokkos::deep_copy(d, hm);
  return d;
}

template <typename T>
View2L<T>
toDevice2(const std::string & name, const std::vector<T> & h, int n0, int n1)
{
  View2L<T> d(name, n0, n1);
  auto hm = Kokkos::create_mirror_view(d);
  for (int j = 0; j < n1; ++j)
    for (int i = 0; i < n0; ++i)
      hm(i, j) = h[i + size_t(n0) * j];
  Kokkos::deep_copy(d, hm);
  return d;
}

template <typename T>
View2R<T>
toDevice2R(const std::string & name, const std::vector<T> & h, int n0, int n1)
{
  View2R<T> d(name, n0, n1);
  auto hm = Kokkos::create_mirror_view(d);
  for (int j = 0; j < n1; ++j)
    for (int i = 0; i < n0; ++i)
      hm(i, j) = h[i + size_t(n0) * j];
  Kokkos::deep_copy(d, hm);
  return d;
}

Data
build(const HostMesh & m,
      const std::string & cell_order,
      const std::string & face_order,
      HostInfo & H)
{
  Data D;
  const auto perm = cellOrdering(m, cell_order);
  D.nc = m.nc;
  D.nf = m.nf;
  H.new_to_old_cell.resize(m.nc);
  for (int o = 0; o < m.nc; ++o)
    H.new_to_old_cell[perm[o]] = o;

  std::vector<double> vol(m.nc);
  std::vector<R3> cen(m.nc);
  for (int o = 0; o < m.nc; ++o)
    vol[perm[o]] = m.vol[o], cen[perm[o]] = m.cen[o];

  // face order: internal first, then boundary grouped by boundary id
  std::vector<int> forder(m.nf);
  std::iota(forder.begin(), forder.end(), 0);
  auto cat = [&](int f) { return m.fn[f] >= 0 ? -1 : m.bid[f]; };
  if (face_order == "native")
    std::stable_sort(forder.begin(), forder.end(), [&](int a, int b) { return cat(a) < cat(b); });
  else if (face_order == "owner")
    std::stable_sort(forder.begin(),
                     forder.end(),
                     [&](int a, int b)
                     {
                       auto key = [&](int f)
                       {
                         const int e = perm[m.fe[f]];
                         const int n = m.fn[f] >= 0 ? perm[m.fn[f]] : e;
                         return std::make_tuple(cat(f), std::min(e, n), std::max(e, n));
                       };
                       return key(a) < key(b);
                     });
  else
    throw std::runtime_error("unknown face ordering " + face_order);

  std::vector<int> fe(m.nf), fn(m.nf), fse(m.nf), fsn(m.nf), fbid(m.nf);
  std::vector<double> farea(m.nf), fdcn(m.nf), fdcf(m.nf), fgc(m.nf);
  std::vector<R3> fnrm(m.nf), fc(m.nf);
  int max_side = 0;
  for (int i = 0; i < m.nf; ++i)
  {
    const int f = forder[i];
    fe[i] = perm[m.fe[f]];
    fn[i] = m.fn[f] >= 0 ? perm[m.fn[f]] : -1;
    fse[i] = m.se[f];
    fsn[i] = m.sn[f];
    fbid[i] = m.bid[f];
    farea[i] = m.area[f];
    fgc[i] = m.gc[f];
    fnrm[i] = m.nrm[f];
    fc[i] = m.fc[f];
    const R3 d = m.dcn[f];
    const double mag = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
    fdcn[i] = mag; // |C_N - C_E| for internal faces
    fdcf[i] = mag; // |x_f - C_E| for boundary faces (FaceInfo stores that in dCN)
    max_side = std::max({max_side, m.se[f], m.sn[f]});
    if (m.fn[f] >= 0)
      D.n_internal = i + 1;
  }
  D.max_sides = max_side + 1;
  const int S = D.max_sides, NC = m.nc;

  // element -> faces and padded directed-side geometry
  std::vector<int> ef(size_t(S) * NC, -1), ds_nb(size_t(S) * NC, -1), ds_bid(size_t(S) * NC, -1);
  std::vector<double> ds_area(size_t(S) * NC, 0), ds_dcn(size_t(S) * NC, 0),
      ds_dcf(size_t(S) * NC, 0), ds_gc(size_t(S) * NC, 0);
  std::vector<R3> ds_nrm(size_t(S) * NC, R3{0, 0, 0}), ds_fc(size_t(S) * NC, R3{0, 0, 0});
  std::vector<int> v0_ic, v0_is, v0_bc, v0_bs;
  for (int f = 0; f < m.nf; ++f)
  {
    const size_t hE = fse[f] + size_t(S) * fe[f];
    ef[hE] = f;
    ds_area[hE] = farea[f];
    ds_gc[hE] = fgc[f];
    ds_nrm[hE] = fnrm[f];
    ds_fc[hE] = fc[f];
    if (fn[f] >= 0)
    {
      const size_t hN = fsn[f] + size_t(S) * fn[f];
      ef[hN] = f;
      ds_area[hN] = farea[f];
      ds_gc[hN] = 1 - fgc[f];
      ds_nrm[hN] = R3{-fnrm[f].x, -fnrm[f].y, -fnrm[f].z};
      ds_fc[hN] = fc[f];
      ds_nb[hE] = fn[f];
      ds_nb[hN] = fe[f];
      ds_dcn[hE] = ds_dcn[hN] = fdcn[f];
    }
    else
    {
      ds_bid[hE] = fbid[f];
      ds_dcf[hE] = fdcf[f];
    }
  }
  // V0 worklists sorted by (cell, side), as in KokkosLinearFVFluxKernel::initialSetup
  for (int c = 0; c < NC; ++c)
    for (int s = 0; s < S; ++s)
    {
      const size_t h = s + size_t(S) * c;
      if (ef[h] < 0)
        continue;
      if (ds_nb[h] >= 0)
        v0_ic.push_back(c), v0_is.push_back(s);
      else
        v0_bc.push_back(c), v0_bs.push_back(s);
    }

  // CSR with sorted columns
  std::vector<std::vector<int>> cols(NC);
  for (int c = 0; c < NC; ++c)
    cols[c].push_back(c);
  for (int f = 0; f < D.n_internal; ++f)
    cols[fe[f]].push_back(fn[f]), cols[fn[f]].push_back(fe[f]);
  std::vector<int> row_ptr(NC + 1, 0), col;
  for (int c = 0; c < NC; ++c)
  {
    std::sort(cols[c].begin(), cols[c].end());
    row_ptr[c + 1] = row_ptr[c] + cols[c].size();
    col.insert(col.end(), cols[c].begin(), cols[c].end());
  }
  D.nnz = col.size();
  auto slotOf = [&](int r, int c)
  {
    auto b = col.begin() + row_ptr[r], e = col.begin() + row_ptr[r + 1];
    return int(std::lower_bound(b, e, c) - col.begin());
  };
  std::vector<int> diag_slot(NC), slot_EN(m.nf, -1), slot_NE(m.nf, -1), hf_slot(size_t(S) * NC, -1);
  for (int c = 0; c < NC; ++c)
    diag_slot[c] = slotOf(c, c);
  for (int f = 0; f < D.n_internal; ++f)
  {
    slot_EN[f] = slotOf(fe[f], fn[f]);
    slot_NE[f] = slotOf(fn[f], fe[f]);
    hf_slot[fse[f] + size_t(S) * fe[f]] = slot_EN[f];
    hf_slot[fsn[f] + size_t(S) * fn[f]] = slot_NE[f];
  }

  D.vol = toDevice<View1<double>>("vol", vol);
  D.cen = toDevice<View1<R3>>("cen", cen);
  D.fe = toDevice<View1<int>>("fe", fe);
  D.fn = toDevice<View1<int>>("fn", fn);
  D.fse = toDevice<View1<int>>("fse", fse);
  D.fsn = toDevice<View1<int>>("fsn", fsn);
  D.fbid = toDevice<View1<int>>("fbid", fbid);
  D.farea = toDevice<View1<double>>("farea", farea);
  D.fdcn_mag = toDevice<View1<double>>("fdcn", fdcn);
  D.fdcf_mag = toDevice<View1<double>>("fdcf", fdcf);
  D.fnrm = toDevice<View1<R3>>("fnrm", fnrm);
  D.fc = toDevice<View1<R3>>("fc", fc);
  D.elem_faces = toDevice2<int>("elem_faces", ef, S, NC);
  D.row_ptr = toDevice<View1<int>>("row_ptr", row_ptr);
  D.col = toDevice<View1<int>>("col", col);
  D.diag_slot = toDevice<View1<int>>("diag_slot", diag_slot);
  D.slot_EN = toDevice<View1<int>>("slot_EN", slot_EN);
  D.slot_NE = toDevice<View1<int>>("slot_NE", slot_NE);
  D.hf_slot = toDevice2<int>("hf_slot", hf_slot, S, NC);
  D.ds_area = toDevice2<double>("ds_area", ds_area, S, NC);
  D.ds_dcn_mag = toDevice2<double>("ds_dcn", ds_dcn, S, NC);
  D.ds_dcf_mag = toDevice2<double>("ds_dcf", ds_dcf, S, NC);
  D.ds_nrm = toDevice2<R3>("ds_nrm", ds_nrm, S, NC);
  D.ds_fc = toDevice2<R3>("ds_fc", ds_fc, S, NC);
  D.ds_nb = toDevice2<int>("ds_nb", ds_nb, S, NC);
  D.ds_bid = toDevice2<int>("ds_bid", ds_bid, S, NC);
  D.v0_int_cell = toDevice<View1<int>>("v0_ic", v0_ic);
  D.v0_int_side = toDevice<View1<int>>("v0_is", v0_is);
  D.v0_bnd_cell = toDevice<View1<int>>("v0_bc", v0_bc);
  D.v0_bnd_side = toDevice<View1<int>>("v0_bs", v0_bs);
  D.val = View1<double>("val", D.nnz);
  D.rhs = View1<double>("rhs", NC);
  D.fs_dE = View1<double>("fs_dE", m.nf);
  D.fs_dN = View1<double>("fs_dN", m.nf);
  D.fs_rE = View1<double>("fs_rE", m.nf);
  D.fs_rN = View1<double>("fs_rN", m.nf);
  D.hf_d = View2L<double>("hf_d", S, NC);
  D.hf_o = View2L<double>("hf_o", S, NC);
  D.hf_r = View2L<double>("hf_r", S, NC);
  D.fgc = toDevice<View1<double>>("fgc", fgc);
  D.ds_gc = toDevice2<double>("ds_gc", ds_gc, S, NC);
  {
    std::vector<double> phi(NC);
    for (int c = 0; c < NC; ++c)
      phi[c] = phys::phi(cen[c].x, cen[c].y, cen[c].z);
    D.phi = toDevice<View1<double>>("phi", phi);
  }
  D.grad = View1<double>("grad", 3 * size_t(NC));
  D.dsr_area = toDevice2R<double>("dsr_area", ds_area, S, NC);
  D.dsr_gc = toDevice2R<double>("dsr_gc", ds_gc, S, NC);
  D.dsr_nrm = toDevice2R<R3>("dsr_nrm", ds_nrm, S, NC);
  D.dsr_fc = toDevice2R<R3>("dsr_fc", ds_fc, S, NC);
  D.dsr_nb = toDevice2R<int>("dsr_nb", ds_nb, S, NC);
  D.dsr_bid = toDevice2R<int>("dsr_bid", ds_bid, S, NC);
  D.elem_faces_r = toDevice2R<int>("elem_faces_r", ef, S, NC);
  {
    std::vector<double> nx(m.nf), ny(m.nf), nz(m.nf), cx(m.nf), cy(m.nf), cz(m.nf);
    for (int f = 0; f < m.nf; ++f)
    {
      nx[f] = fnrm[f].x, ny[f] = fnrm[f].y, nz[f] = fnrm[f].z;
      cx[f] = fc[f].x, cy[f] = fc[f].y, cz[f] = fc[f].z;
    }
    D.fnx = toDevice<View1<double>>("fnx", nx);
    D.fny = toDevice<View1<double>>("fny", ny);
    D.fnz = toDevice<View1<double>>("fnz", nz);
    D.fcx = toDevice<View1<double>>("fcx", cx);
    D.fcy = toDevice<View1<double>>("fcy", cy);
    D.fcz = toDevice<View1<double>>("fcz", cz);
  }

  const size_t F = m.nf, C = NC, SC = size_t(S) * NC;
  // face SoA used by V1-V4: fe, fn, fse, fsn, fbid, farea, fdcn, fdcf, fnrm, fc
  H.bytes_face_soa = F * (5 * sizeof(int) + 3 * sizeof(double) + 2 * sizeof(R3));
  // padded directed-side geometry used by V0 (as in today's Kokkos mesh) + worklists
  H.bytes_side_padded = SC * (3 * sizeof(double) + 2 * sizeof(R3) + 2 * sizeof(int)) +
                        (v0_ic.size() + v0_bc.size()) * 2 * sizeof(int);
  H.bytes_elem_faces = SC * sizeof(int);
  H.bytes_slots = (C + 2 * F) * sizeof(int) + SC * sizeof(int);
  H.bytes_face_scratch = 4 * F * sizeof(double);
  H.bytes_hf_scratch = 3 * SC * sizeof(double);
  // gradient inputs: G1-G3 read fe, fn, fbid, farea, fgc, fnrm, fc (+ elem_faces for G1);
  // G0 reads ds_nb, ds_bid, ds_area, ds_gc, ds_nrm, ds_fc
  H.bytes_grad_face = F * (3 * sizeof(int) + 2 * sizeof(double) + 2 * sizeof(R3));
  H.bytes_grad_side = SC * (2 * sizeof(int) + 2 * sizeof(double) + 2 * sizeof(R3));
  return D;
}

// ------------------------------------------------------------------------------------------------
// Kernels
// ------------------------------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION int
findCol(const View1<int> & row_ptr, const View1<int> & col, int r, int c)
{
  int lo = row_ptr(r), hi = row_ptr(r + 1);
  while (lo < hi)
  {
    const int mid = (lo + hi) / 2;
    if (col(mid) < c)
      lo = mid + 1;
    else
      hi = mid;
  }
  return lo;
}

KOKKOS_INLINE_FUNCTION double
dot(const R3 & a, double x, double y, double z)
{
  return a.x * x + a.y * y + a.z * z;
}

// Internal-face coefficients of row E (CPU convention; row N is (-a_EN, -a_EE) on (diag, off))
struct FaceCoef
{
  double a_EE, a_EN;
};

KOKKOS_INLINE_FUNCTION FaceCoef
diffusionFace(double kA_over_d)
{
  return {kA_over_d, -kA_over_d};
}

KOKKOS_INLINE_FUNCTION FaceCoef
advectionFace(double flux_A)
{
  return {flux_A > 0 ? flux_A : 0.0, flux_A < 0 ? flux_A : 0.0};
}

using Policy = Kokkos::RangePolicy<ExecSpace>;

// V0: directed sides, per-contribution coefficient evaluation, binary search, atomics
void
assembleV0(Data & D, std::vector<double> & t)
{
  auto d = D;
  Kokkos::Timer tm;
  // BC relation data: (computed inline below as in the relation arrays; cost included in the
  // boundary loops)
  // diffusion internal matrix
  Kokkos::parallel_for(
      "V0_diff_int", Policy(0, d.v0_int_cell.extent(0)), KOKKOS_LAMBDA(int i) {
        const int c = d.v0_int_cell(i), s = d.v0_int_side(i);
        const int n = d.ds_nb(s, c);
        const R3 x = d.ds_fc(s, c);
        // matrix and neighbor contributions each evaluate the coefficient (as today)
        const double a_cc = phys::k(x.x, x.y, x.z) * d.ds_area(s, c) / d.ds_dcn_mag(s, c);
        const double a_cn = -phys::k(x.x, x.y, x.z) * d.ds_area(s, c) / d.ds_dcn_mag(s, c);
        Kokkos::atomic_add(&d.val(findCol(d.row_ptr, d.col, c, c)), a_cc);
        Kokkos::atomic_add(&d.val(findCol(d.row_ptr, d.col, c, n)), a_cn);
      });
  Kokkos::fence();
  t[0] += tm.seconds();
  tm.reset();
  // diffusion boundary matrix + rhs (separate loops in MOOSE; one pass each here)
  Kokkos::parallel_for(
      "V0_diff_bnd_mat", Policy(0, d.v0_bnd_cell.extent(0)), KOKKOS_LAMBDA(int i) {
        const int c = d.v0_bnd_cell(i), s = d.v0_bnd_side(i);
        if (!phys::isDirichlet(d.ds_bid(s, c)))
          return;
        const R3 x = d.ds_fc(s, c);
        const double coef = -1.0 / d.ds_dcf_mag(s, c);
        const double v = -phys::k(x.x, x.y, x.z) * d.ds_area(s, c) * coef;
        Kokkos::atomic_add(&d.val(findCol(d.row_ptr, d.col, c, c)), v);
      });
  Kokkos::parallel_for(
      "V0_diff_bnd_rhs", Policy(0, d.v0_bnd_cell.extent(0)), KOKKOS_LAMBDA(int i) {
        const int c = d.v0_bnd_cell(i), s = d.v0_bnd_side(i);
        if (!phys::isDirichlet(d.ds_bid(s, c)))
          return;
        const R3 x = d.ds_fc(s, c);
        const double src = phys::g(x.x, x.y, x.z) / d.ds_dcf_mag(s, c);
        Kokkos::atomic_add(&d.rhs(c), phys::k(x.x, x.y, x.z) * d.ds_area(s, c) * src);
      });
  Kokkos::fence();
  t[1] += tm.seconds();
  tm.reset();
  // advection internal matrix
  Kokkos::parallel_for(
      "V0_adv_int", Policy(0, d.v0_int_cell.extent(0)), KOKKOS_LAMBDA(int i) {
        const int c = d.v0_int_cell(i), s = d.v0_int_side(i);
        const int n = d.ds_nb(s, c);
        const double f1 = dot(d.ds_nrm(s, c), phys::vx, phys::vy, phys::vz);
        const double a_cc = f1 > 0 ? f1 * d.ds_area(s, c) : 0.0;
        const double f2 = dot(d.ds_nrm(s, c), phys::vx, phys::vy, phys::vz);
        const double a_cn = f2 < 0 ? f2 * d.ds_area(s, c) : 0.0;
        if (a_cc != 0)
          Kokkos::atomic_add(&d.val(findCol(d.row_ptr, d.col, c, c)), a_cc);
        if (a_cn != 0)
          Kokkos::atomic_add(&d.val(findCol(d.row_ptr, d.col, c, n)), a_cn);
      });
  Kokkos::fence();
  t[0] += tm.seconds();
  tm.reset();
  // advection boundary matrix + rhs
  Kokkos::parallel_for(
      "V0_adv_bnd_mat", Policy(0, d.v0_bnd_cell.extent(0)), KOKKOS_LAMBDA(int i) {
        const int c = d.v0_bnd_cell(i), s = d.v0_bnd_side(i);
        const double coef = phys::isDirichlet(d.ds_bid(s, c)) ? 0.0 : 1.0;
        const double v = dot(d.ds_nrm(s, c), phys::vx, phys::vy, phys::vz) * d.ds_area(s, c) * coef;
        if (v != 0)
          Kokkos::atomic_add(&d.val(findCol(d.row_ptr, d.col, c, c)), v);
      });
  Kokkos::parallel_for(
      "V0_adv_bnd_rhs", Policy(0, d.v0_bnd_cell.extent(0)), KOKKOS_LAMBDA(int i) {
        const int c = d.v0_bnd_cell(i), s = d.v0_bnd_side(i);
        if (!phys::isDirichlet(d.ds_bid(s, c)))
          return;
        const R3 x = d.ds_fc(s, c);
        const double v = -dot(d.ds_nrm(s, c), phys::vx, phys::vy, phys::vz) * d.ds_area(s, c) *
                         phys::g(x.x, x.y, x.z);
        Kokkos::atomic_add(&d.rhs(c), v);
      });
  Kokkos::fence();
  t[1] += tm.seconds();
  tm.reset();
  // source
  Kokkos::parallel_for(
      "V0_src", Policy(0, d.nc), KOKKOS_LAMBDA(int c) {
        const R3 x = d.cen(c);
        d.rhs(c) += phys::s(x.x, x.y, x.z) * d.vol(c);
      });
  Kokkos::fence();
  t[2] += tm.seconds();
}

// Shared face computation for V1-V4: all coefficients of one face, evaluated once
struct FaceResult
{
  double dE, oEN, dN, oNE, rE, rN;
};

// kernel: 0 = diffusion, 1 = advection
template <int Kernel>
KOKKOS_INLINE_FUNCTION FaceResult
faceTerms(const Data & d, int f)
{
  FaceResult r{0, 0, 0, 0, 0, 0};
  const double A = d.farea(f);
  if (f < d.n_internal)
  {
    FaceCoef a;
    if (Kernel == 0)
    {
      const R3 x = d.fc(f);
      a = diffusionFace(phys::k(x.x, x.y, x.z) * A / d.fdcn_mag(f));
    }
    else
      a = advectionFace(dot(d.fnrm(f), phys::vx, phys::vy, phys::vz) * A);
    r.dE = a.a_EE;
    r.oEN = a.a_EN;
    r.dN = -a.a_EN;
    r.oNE = -a.a_EE;
  }
  else
  {
    const bool dir = phys::isDirichlet(d.fbid(f));
    const R3 x = d.fc(f);
    if (Kernel == 0)
    {
      if (dir)
      {
        const double kA = phys::k(x.x, x.y, x.z) * A;
        r.dE = kA / d.fdcf_mag(f);
        r.rE = kA * phys::g(x.x, x.y, x.z) / d.fdcf_mag(f);
      }
    }
    else
    {
      const double fluxA = dot(d.fnrm(f), phys::vx, phys::vy, phys::vz) * A;
      if (dir)
        r.rE = -fluxA * phys::g(x.x, x.y, x.z);
      else
        r.dE = fluxA;
    }
  }
  return r;
}

// V1/V2: face pass with atomics (search or precomputed slots)
template <bool Slots>
void
assembleFaceAtomic(Data & D, std::vector<double> & t)
{
  auto d = D;
  Kokkos::Timer tm;
  auto pass = [&](auto kernel_tag)
  {
    constexpr int K = decltype(kernel_tag)::value;
    Kokkos::parallel_for(
        "face_atomic", Policy(0, d.nf), KOKKOS_LAMBDA(int f) {
          const auto r = faceTerms<K>(d, f);
          const int E = d.fe(f);
          const int sEE = Slots ? d.diag_slot(E) : findCol(d.row_ptr, d.col, E, E);
          if (r.dE != 0)
            Kokkos::atomic_add(&d.val(sEE), r.dE);
          if (r.rE != 0)
            Kokkos::atomic_add(&d.rhs(E), r.rE);
          if (f < d.n_internal)
          {
            const int N = d.fn(f);
            const int sEN = Slots ? d.slot_EN(f) : findCol(d.row_ptr, d.col, E, N);
            const int sNN = Slots ? d.diag_slot(N) : findCol(d.row_ptr, d.col, N, N);
            const int sNE = Slots ? d.slot_NE(f) : findCol(d.row_ptr, d.col, N, E);
            if (r.oEN != 0)
              Kokkos::atomic_add(&d.val(sEN), r.oEN);
            if (r.dN != 0)
              Kokkos::atomic_add(&d.val(sNN), r.dN);
            if (r.oNE != 0)
              Kokkos::atomic_add(&d.val(sNE), r.oNE);
          }
        });
  };
  pass(std::integral_constant<int, 0>{});
  pass(std::integral_constant<int, 1>{});
  Kokkos::fence();
  t[0] += tm.seconds();
  tm.reset();
  Kokkos::parallel_for(
      "src", Policy(0, d.nc), KOKKOS_LAMBDA(int c) {
        const R3 x = d.cen(c);
        d.rhs(c) += phys::s(x.x, x.y, x.z) * d.vol(c);
      });
  Kokkos::fence();
  t[2] += tm.seconds();
}

// V2F: V2 with the diffusion and advection face passes fused into one kernel
void
assembleFaceAtomicFused(Data & D, std::vector<double> & t)
{
  auto d = D;
  Kokkos::Timer tm;
  Kokkos::parallel_for(
      "face_atomic_fused", Policy(0, d.nf), KOKKOS_LAMBDA(int f) {
        const auto a = faceTerms<0>(d, f);
        const auto b = faceTerms<1>(d, f);
        const int E = d.fe(f);
        const double dE = a.dE + b.dE, rE = a.rE + b.rE;
        if (dE != 0)
          Kokkos::atomic_add(&d.val(d.diag_slot(E)), dE);
        if (rE != 0)
          Kokkos::atomic_add(&d.rhs(E), rE);
        if (f < d.n_internal)
        {
          const int N = d.fn(f);
          const double oEN = a.oEN + b.oEN, dN = a.dN + b.dN, oNE = a.oNE + b.oNE;
          if (oEN != 0)
            Kokkos::atomic_add(&d.val(d.slot_EN(f)), oEN);
          if (dN != 0)
            Kokkos::atomic_add(&d.val(d.diag_slot(N)), dN);
          if (oNE != 0)
            Kokkos::atomic_add(&d.val(d.slot_NE(f)), oNE);
        }
      });
  Kokkos::fence();
  t[0] += tm.seconds();
  tm.reset();
  Kokkos::parallel_for(
      "src", Policy(0, d.nc), KOKKOS_LAMBDA(int c) {
        const R3 x = d.cen(c);
        d.rhs(c) += phys::s(x.x, x.y, x.z) * d.vol(c);
      });
  Kokkos::fence();
  t[2] += tm.seconds();
}

// V3: face pass into per-face scratch + direct off-diagonal slot writes; per-cell gather
void
assembleFaceGather(Data & D, std::vector<double> & t)
{
  auto d = D;
  Kokkos::Timer tm;
  Kokkos::parallel_for(
      "V3_zero", Policy(0, d.nf), KOKKOS_LAMBDA(int f) {
        d.fs_dE(f) = d.fs_dN(f) = d.fs_rE(f) = d.fs_rN(f) = 0;
      });
  auto pass = [&](auto kernel_tag)
  {
    constexpr int K = decltype(kernel_tag)::value;
    Kokkos::parallel_for(
        "V3_face", Policy(0, d.nf), KOKKOS_LAMBDA(int f) {
          const auto r = faceTerms<K>(d, f);
          d.fs_dE(f) += r.dE;
          d.fs_rE(f) += r.rE;
          if (f < d.n_internal)
          {
            d.fs_dN(f) += r.dN;
            d.fs_rN(f) += r.rN;
            d.val(d.slot_EN(f)) += r.oEN; // each slot belongs to exactly one face
            d.val(d.slot_NE(f)) += r.oNE;
          }
        });
  };
  pass(std::integral_constant<int, 0>{});
  pass(std::integral_constant<int, 1>{});
  Kokkos::fence();
  t[0] += tm.seconds();
  tm.reset();
  const int S = d.max_sides;
  Kokkos::parallel_for(
      "V3_gather", Policy(0, d.nc), KOKKOS_LAMBDA(int c) {
        double diag = 0, b = 0;
        for (int s = 0; s < S; ++s)
        {
          const int f = d.elem_faces(s, c);
          if (f < 0)
            continue;
          const bool is_e = d.fe(f) == c;
          diag += is_e ? d.fs_dE(f) : d.fs_dN(f);
          b += is_e ? d.fs_rE(f) : d.fs_rN(f);
        }
        const R3 x = d.cen(c);
        b += phys::s(x.x, x.y, x.z) * d.vol(c);
        d.val(d.diag_slot(c)) = diag;
        d.rhs(c) = b;
      });
  Kokkos::fence();
  t[2] += tm.seconds();
}

// V4: face pass into half-face scratch; per-cell gather writes the whole row
void
assembleHalfFaceGather(Data & D, std::vector<double> & t)
{
  auto d = D;
  Kokkos::Timer tm;
  Kokkos::deep_copy(d.hf_d, 0.0);
  Kokkos::deep_copy(d.hf_o, 0.0);
  Kokkos::deep_copy(d.hf_r, 0.0);
  auto pass = [&](auto kernel_tag)
  {
    constexpr int K = decltype(kernel_tag)::value;
    Kokkos::parallel_for(
        "V4_face", Policy(0, d.nf), KOKKOS_LAMBDA(int f) {
          const auto r = faceTerms<K>(d, f);
          const int E = d.fe(f), sE = d.fse(f);
          d.hf_d(sE, E) += r.dE;
          d.hf_r(sE, E) += r.rE;
          if (f < d.n_internal)
          {
            const int N = d.fn(f), sN = d.fsn(f);
            d.hf_o(sE, E) += r.oEN;
            d.hf_d(sN, N) += r.dN;
            d.hf_o(sN, N) += r.oNE;
            d.hf_r(sN, N) += r.rN;
          }
        });
  };
  pass(std::integral_constant<int, 0>{});
  pass(std::integral_constant<int, 1>{});
  Kokkos::fence();
  t[0] += tm.seconds();
  tm.reset();
  const int S = d.max_sides;
  Kokkos::parallel_for(
      "V4_gather", Policy(0, d.nc), KOKKOS_LAMBDA(int c) {
        double diag = 0, b = 0;
        for (int s = 0; s < S; ++s)
        {
          diag += d.hf_d(s, c);
          b += d.hf_r(s, c);
          const int k = d.hf_slot(s, c);
          if (k >= 0)
            d.val(k) = d.hf_o(s, c);
        }
        const R3 x = d.cen(c);
        b += phys::s(x.x, x.y, x.z) * d.vol(c);
        d.val(d.diag_slot(c)) = diag;
        d.rhs(c) = b;
      });
  Kokkos::fence();
  t[2] += tm.seconds();
}

// ------------------------------------------------------------------------------------------------
// Green-Gauss cell gradients
// ------------------------------------------------------------------------------------------------
// Value of the field on a boundary face: the Dirichlet value, else the cell value
KOKKOS_INLINE_FUNCTION double
boundaryFaceValue(int bid, const R3 & x, double phi_c)
{
  return phys::isDirichlet(bid) ? phys::phi(x.x, x.y, x.z) : phi_c;
}

// G0/G0C: per-cell gather over the padded per-(side, cell) geometry (side- or cell-fastest)
template <typename VI, typename VD, typename VR3>
void
gradientSideGather(Data & D, VI ds_nb, VI ds_bid, VD ds_area, VD ds_gc, VR3 ds_nrm, VR3 ds_fc)
{
  auto d = D;
  const int S = d.max_sides;
  Kokkos::parallel_for(
      "G0_side_gather", Policy(0, d.nc), KOKKOS_LAMBDA(int c) {
        const double phi_c = d.phi(c);
        double gx = 0, gy = 0, gz = 0;
        for (int s = 0; s < S; ++s)
        {
          const int nb = ds_nb(s, c);
          double pf;
          if (nb >= 0)
          {
            const double w = ds_gc(s, c);
            pf = w * phi_c + (1 - w) * d.phi(nb);
          }
          else
          {
            const int bid = ds_bid(s, c);
            if (bid < 0)
              continue; // padding: this cell has fewer sides
            pf = boundaryFaceValue(bid, ds_fc(s, c), phi_c);
          }
          const double w = pf * ds_area(s, c);
          const R3 n = ds_nrm(s, c);
          gx += w * n.x, gy += w * n.y, gz += w * n.z;
        }
        const double inv = 1 / d.vol(c);
        d.grad(3 * c) = gx * inv, d.grad(3 * c + 1) = gy * inv, d.grad(3 * c + 2) = gz * inv;
      });
  Kokkos::fence();
}

// G1/G1C/G1CS: per-cell gather over the face arrays through elem_faces (side- or cell-fastest),
// with face normals and centroids as R3 (AoS) or as separate x, y, z arrays (SoA)
template <bool SoA, typename VEF>
void
gradientFaceGather(Data & D, VEF elem_faces)
{
  auto d = D;
  const int S = d.max_sides;
  Kokkos::parallel_for(
      "G1_face_gather", Policy(0, d.nc), KOKKOS_LAMBDA(int c) {
        const double phi_c = d.phi(c);
        double gx = 0, gy = 0, gz = 0;
        for (int s = 0; s < S; ++s)
        {
          const int f = elem_faces(s, c);
          if (f < 0)
            continue;
          double pf;
          if (f < d.n_internal)
          {
            const double w = d.fgc(f);
            pf = w * d.phi(d.fe(f)) + (1 - w) * d.phi(d.fn(f));
          }
          else
          {
            const R3 x = SoA ? R3{d.fcx(f), d.fcy(f), d.fcz(f)} : d.fc(f);
            pf = boundaryFaceValue(d.fbid(f), x, phi_c);
          }
          const double w = (d.fe(f) == c ? pf : -pf) * d.farea(f);
          const R3 n = SoA ? R3{d.fnx(f), d.fny(f), d.fnz(f)} : d.fnrm(f);
          gx += w * n.x, gy += w * n.y, gz += w * n.z;
        }
        const double inv = 1 / d.vol(c);
        d.grad(3 * c) = gx * inv, d.grad(3 * c + 1) = gy * inv, d.grad(3 * c + 2) = gz * inv;
      });
  Kokkos::fence();
}

// G2 (FoldVolume = false): face pass with atomics, then a volume pass
// G3 (FoldVolume = true): 1/V folded into the atomic adds, single face pass
template <bool FoldVolume>
void
gradientFaceAtomic(Data & D)
{
  auto d = D;
  Kokkos::deep_copy(d.grad, 0.0);
  Kokkos::parallel_for(
      "G2_G3_face_atomic", Policy(0, d.nf), KOKKOS_LAMBDA(int f) {
        const int E = d.fe(f);
        const double phi_E = d.phi(E);
        double pf;
        if (f < d.n_internal)
        {
          const double w = d.fgc(f);
          pf = w * phi_E + (1 - w) * d.phi(d.fn(f));
        }
        else
          pf = boundaryFaceValue(d.fbid(f), d.fc(f), phi_E);
        const R3 n = d.fnrm(f);
        const double w = pf * d.farea(f);
        const double sE = FoldVolume ? w / d.vol(E) : w;
        Kokkos::atomic_add(&d.grad(3 * E), sE * n.x);
        Kokkos::atomic_add(&d.grad(3 * E + 1), sE * n.y);
        Kokkos::atomic_add(&d.grad(3 * E + 2), sE * n.z);
        if (f < d.n_internal)
        {
          const int N = d.fn(f);
          const double sN = FoldVolume ? w / d.vol(N) : w;
          Kokkos::atomic_add(&d.grad(3 * N), -sN * n.x);
          Kokkos::atomic_add(&d.grad(3 * N + 1), -sN * n.y);
          Kokkos::atomic_add(&d.grad(3 * N + 2), -sN * n.z);
        }
      });
  if (!FoldVolume)
    Kokkos::parallel_for(
        "G2_volume", Policy(0, d.nc), KOKKOS_LAMBDA(int c) {
          const double inv = 1 / d.vol(c);
          d.grad(3 * c) *= inv, d.grad(3 * c + 1) *= inv, d.grad(3 * c + 2) *= inv;
        });
  Kokkos::fence();
}

// Serial host reference in the original cell numbering, as LinearFVGradientManager computes it
std::vector<double>
referenceGradient(const HostMesh & m)
{
  std::vector<double> phi(m.nc), g(3 * size_t(m.nc), 0.0);
  for (int64_t c = 0; c < m.nc; ++c)
    phi[c] = phys::phi(m.cen[c].x, m.cen[c].y, m.cen[c].z);
  for (int64_t f = 0; f < m.nf; ++f)
  {
    const int64_t E = m.fe[f];
    double pf;
    if (m.fn[f] >= 0)
      pf = m.gc[f] * phi[E] + (1 - m.gc[f]) * phi[m.fn[f]];
    else
      pf = boundaryFaceValue(m.bid[f], m.fc[f], phi[E]);
    const double w = pf * m.area[f];
    const double c[3] = {w * m.nrm[f].x, w * m.nrm[f].y, w * m.nrm[f].z};
    for (int d = 0; d < 3; ++d)
    {
      g[3 * E + d] += c[d];
      if (m.fn[f] >= 0)
        g[3 * m.fn[f] + d] -= c[d];
    }
  }
  for (int64_t c = 0; c < m.nc; ++c)
    for (int d = 0; d < 3; ++d)
      g[3 * c + d] /= m.vol[c];
  return g;
}

// Max |difference| relative to the largest reference gradient entry
double
gradientDiff(const std::vector<double> & ref, const Data & D, const HostInfo & H)
{
  auto g = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D.grad);
  double scale = 0, worst = 0;
  for (double x : ref)
    scale = std::max(scale, std::abs(x));
  for (int c = 0; c < D.nc; ++c)
    for (int d = 0; d < 3; ++d)
      worst = std::max(worst, std::abs(g(3 * c + d) - ref[3 * size_t(H.new_to_old_cell[c]) + d]));
  return worst / scale;
}

uint64_t
hashGradient(const Data & D)
{
  auto g = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D.grad);
  uint64_t h = 1469598103934665603ULL;
  auto c = reinterpret_cast<const unsigned char *>(g.data());
  for (size_t i = 0; i < g.size() * sizeof(double); ++i)
    h = (h ^ c[i]) * 1099511628211ULL;
  return h;
}

// ------------------------------------------------------------------------------------------------
// Comparison helpers (host), in the original cell numbering
// ------------------------------------------------------------------------------------------------
struct HostSystem
{
  std::map<std::pair<int64_t, int64_t>, double> A; // (row, col) in original cell indices
  std::vector<double> b;
};

HostSystem
toHost(const Data & D, const HostInfo & H)
{
  auto rp = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D.row_ptr);
  auto cl = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D.col);
  auto v = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D.val);
  auto r = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D.rhs);
  HostSystem h;
  h.b.resize(D.nc);
  for (int c = 0; c < D.nc; ++c)
  {
    const int oc = H.new_to_old_cell[c];
    h.b[oc] = r(c);
    for (int k = rp(c); k < rp(c + 1); ++k)
      if (v(k) != 0)
        h.A[{oc, H.new_to_old_cell[cl(k)]}] += v(k);
  }
  return h;
}

uint64_t
hashSystem(const Data & D)
{
  auto v = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D.val);
  auto r = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), D.rhs);
  uint64_t h = 1469598103934665603ULL;
  auto mix = [&](const void * p, size_t n)
  {
    auto c = static_cast<const unsigned char *>(p);
    for (size_t i = 0; i < n; ++i)
      h = (h ^ c[i]) * 1099511628211ULL;
  };
  mix(v.data(), v.size() * sizeof(double));
  mix(r.data(), r.size() * sizeof(double));
  return h;
}

double
maxRelDiff(const HostSystem & a, const HostSystem & b)
{
  double scale = 0, worst = 0;
  for (auto & kv : a.A)
    scale = std::max(scale, std::abs(kv.second));
  std::set<std::pair<int64_t, int64_t>> keys;
  for (auto & kv : a.A)
    keys.insert(kv.first);
  for (auto & kv : b.A)
    keys.insert(kv.first);
  for (auto & k : keys)
  {
    const double x = a.A.count(k) ? a.A.at(k) : 0.0, y = b.A.count(k) ? b.A.at(k) : 0.0;
    worst = std::max(worst, std::abs(x - y) / std::max({std::abs(x), std::abs(y), 1e-14 * scale}));
  }
  double bs = 0;
  for (double x : a.b)
    bs = std::max(bs, std::abs(x));
  for (size_t i = 0; i < a.b.size(); ++i)
    worst = std::max(worst,
                     std::abs(a.b[i] - b.b[i]) /
                         std::max({std::abs(a.b[i]), std::abs(b.b[i]), 1e-14 * bs}));
  return worst;
}

// PETSc binary (big-endian) reader for the MOOSE-assembled reference, in dof numbering
template <typename T>
T
be(const unsigned char * p)
{
  unsigned char tmp[sizeof(T)];
  for (size_t i = 0; i < sizeof(T); ++i)
    tmp[i] = p[sizeof(T) - 1 - i];
  T v;
  std::memcpy(&v, tmp, sizeof(T));
  return v;
}

bool
readPetscSystem(const std::string & mat,
                const std::string & vec,
                const HostMesh & m,
                HostSystem & out)
{
  auto slurp = [](const std::string & p)
  {
    std::ifstream in(p, std::ios::binary);
    return std::vector<unsigned char>((std::istreambuf_iterator<char>(in)), {});
  };
  auto A = slurp(mat), B = slurp(vec);
  if (A.empty() || B.empty())
    return false;
  const bool i64 = be<int32_t>(A.data()) != 1211216;
  const size_t is = i64 ? 8 : 4;
  auto geti = [&](const unsigned char * p)
  { return i64 ? be<int64_t>(p) : int64_t(be<int32_t>(p)); };
  const int64_t M = geti(A.data() + is), nz = geti(A.data() + 3 * is);
  const unsigned char * p = A.data() + 4 * is;
  std::vector<int64_t> rl(M), cols(nz);
  for (int64_t i = 0; i < M; ++i, p += is)
    rl[i] = geti(p);
  for (int64_t i = 0; i < nz; ++i, p += is)
    cols[i] = geti(p);
  std::map<int64_t, int64_t> dof_to_cell;
  for (int64_t c = 0; c < m.nc; ++c)
    dof_to_cell[m.dof[c]] = c;
  int64_t k = 0;
  for (int64_t r = 0; r < M; ++r)
    for (int64_t j = 0; j < rl[r]; ++j, ++k, p += 8)
    {
      const double v = be<double>(p);
      if (v != 0)
        out.A[{dof_to_cell.at(r), dof_to_cell.at(cols[k])}] += v;
    }
  const unsigned char * q = B.data() + 2 * is;
  out.b.assign(m.nc, 0.0);
  for (int64_t r = 0; r < M; ++r, q += 8)
    out.b[dof_to_cell.at(r)] = be<double>(q);
  return true;
}

// ------------------------------------------------------------------------------------------------
// Launch cost
// ------------------------------------------------------------------------------------------------
template <int Bytes>
struct PaddedFunctor
{
  View1<double> out;
  char pad[Bytes];
  KOKKOS_INLINE_FUNCTION void operator()(int i) const { out(i) += pad[i % Bytes]; }
};

template <int Bytes>
void
launchBench(int n, int launches)
{
  PaddedFunctor<Bytes> f;
  f.out = View1<double>("out", n);
  for (int i = 0; i < Bytes; ++i)
    f.pad[i] = 0;
  for (const bool fence_each : {true, false})
  {
    for (int w = 0; w < 10; ++w)
      Kokkos::parallel_for("launch", Policy(0, n), f);
    Kokkos::fence();
    Kokkos::Timer tm;
    for (int l = 0; l < launches; ++l)
    {
      Kokkos::parallel_for("launch", Policy(0, n), f);
      if (fence_each)
        Kokkos::fence();
    }
    Kokkos::fence();
    std::printf("L,%zu,%d,%d,%.3f\n", sizeof(f), n, int(fence_each), tm.seconds() / launches * 1e6);
  }
}

// ------------------------------------------------------------------------------------------------
int
main(int argc, char ** argv)
{
  Kokkos::initialize(argc, argv);
  int rc = 0;
  {
    std::string mesh_file, ref_mat, ref_vec, cell_order = "native", face_order = "native";
    std::string variants = "V0,V1,V2,V2F,V3,V4", grad_variants = "";
    int reps = 20;
    bool launch_bench = false;
    for (int i = 1; i < argc; ++i)
    {
      std::string a = argv[i];
      auto next = [&]() { return std::string(argv[++i]); };
      if (a == "--mesh")
        mesh_file = next();
      else if (a == "--ref-mat")
        ref_mat = next();
      else if (a == "--ref-rhs")
        ref_vec = next();
      else if (a == "--cell-order")
        cell_order = next();
      else if (a == "--face-order")
        face_order = next();
      else if (a == "--variants")
        variants = next() == "none" ? std::string() : std::string(argv[i]);
      else if (a == "--grad-variants")
        grad_variants = next();
      else if (a == "--reps")
        reps = std::stoi(next());
      else if (a == "--launch-bench")
        launch_bench = true;
    }
    if (launch_bench)
    {
      std::printf("# launch cost, threads=%d\nlaunch,functor_bytes,n,fence_each,us_per_launch\n",
                  ExecSpace().concurrency());
      for (const int n : {1, 100000, 2000000})
      {
        launchBench<64>(n, 2000);
        launchBench<6144>(n, 2000);
      }
      // No views are alive in this scope yet
      Kokkos::finalize();
      return 0;
    }
    if (mesh_file.empty())
    {
      std::fprintf(stderr,
                   "usage: fv_layout_bench --mesh FILE [--ref-mat A.bin --ref-rhs b.bin] "
                   "[--cell-order native|morton|rcm|scrambled] [--face-order native|owner] "
                   "[--variants V0,V1,...|none] [--grad-variants G0,G1,G2,G3] [--reps N] "
                   "[--kokkos-num-threads N]\n");
      Kokkos::finalize();
      return 2;
    }

    const auto m = readMesh(mesh_file);
    HostInfo H;
    auto D = build(m, cell_order, face_order, H);

    std::printf("# mesh=%s cells=%lld faces=%lld internal=%d max_sides=%d nnz=%d threads=%d "
                "cell_order=%s face_order=%s\n",
                mesh_file.c_str(),
                (long long)m.nc,
                (long long)m.nf,
                D.n_internal,
                D.max_sides,
                D.nnz,
                ExecSpace().concurrency(),
                cell_order.c_str(),
                face_order.c_str());
    std::printf("# bytes/cell: face_soa=%.1f side_padded(V0)=%.1f elem_faces=%.1f slots=%.1f "
                "face_scratch(V3)=%.1f halfface_scratch(V4)=%.1f\n",
                double(H.bytes_face_soa) / m.nc,
                double(H.bytes_side_padded) / m.nc,
                double(H.bytes_elem_faces) / m.nc,
                double(H.bytes_slots) / m.nc,
                double(H.bytes_face_scratch) / m.nc,
                double(H.bytes_hf_scratch) / m.nc);

    auto run = [&](const std::string & v, std::vector<double> & t)
    {
      Kokkos::deep_copy(D.val, 0.0);
      Kokkos::deep_copy(D.rhs, 0.0);
      Kokkos::fence();
      if (v == "V0")
        assembleV0(D, t);
      else if (v == "V1")
        assembleFaceAtomic<false>(D, t);
      else if (v == "V2")
        assembleFaceAtomic<true>(D, t);
      else if (v == "V2F")
        assembleFaceAtomicFused(D, t);
      else if (v == "V3")
        assembleFaceGather(D, t);
      else if (v == "V4")
        assembleHalfFaceGather(D, t);
      else
        throw std::runtime_error("unknown variant " + v);
    };

    // Gradient variants
    if (!grad_variants.empty())
    {
      std::printf("# gradient bytes/cell: face_soa=%.1f (+ elem_faces %.1f) side_padded=%.1f "
                  "grad=24.0\n",
                  double(H.bytes_grad_face) / m.nc,
                  double(H.bytes_elem_faces) / m.nc,
                  double(H.bytes_grad_side) / m.nc);
      const auto gref = referenceGradient(m);
      std::printf("gvariant,max_rel_diff_vs_ref,distinct_hashes,ms,ns_per_cell\n");
      size_t gpos = 0;
      while (gpos < grad_variants.size())
      {
        const auto comma = grad_variants.find(',', gpos);
        const auto v = grad_variants.substr(
            gpos, comma == std::string::npos ? std::string::npos : comma - gpos);
        gpos = comma == std::string::npos ? grad_variants.size() : comma + 1;

        auto run_grad = [&]()
        {
          if (v == "G0")
            gradientSideGather(D, D.ds_nb, D.ds_bid, D.ds_area, D.ds_gc, D.ds_nrm, D.ds_fc);
          else if (v == "G0C")
            gradientSideGather(D, D.dsr_nb, D.dsr_bid, D.dsr_area, D.dsr_gc, D.dsr_nrm, D.dsr_fc);
          else if (v == "G1")
            gradientFaceGather<false>(D, D.elem_faces);
          else if (v == "G1C")
            gradientFaceGather<false>(D, D.elem_faces_r);
          else if (v == "G1CS")
            gradientFaceGather<true>(D, D.elem_faces_r);
          else if (v == "G2")
            gradientFaceAtomic<false>(D);
          else if (v == "G3")
            gradientFaceAtomic<true>(D);
          else
            throw std::runtime_error("unknown gradient variant " + v);
        };
        run_grad();
        run_grad();
        const double diff = gradientDiff(gref, D, H);
        std::vector<double> times;
        std::set<uint64_t> hashes;
        for (int r = 0; r < reps; ++r)
        {
          Kokkos::fence();
          Kokkos::Timer tm;
          run_grad();
          times.push_back(tm.seconds());
          hashes.insert(hashGradient(D));
        }
        std::nth_element(times.begin(), times.begin() + times.size() / 2, times.end());
        const double t = times[times.size() / 2];
        std::printf(
            "%s,%.3e,%zu,%.3f,%.2f\n", v.c_str(), diff, hashes.size(), t * 1e3, t / m.nc * 1e9);
        if (diff > 1e-12)
          rc = 1;
      }
    }

    if (!variants.empty())
    {
      // Reference: V0
      std::vector<double> tmp(3, 0.0);
      run("V0", tmp);
      const auto ref = toHost(D, H);
      if (!ref_mat.empty())
      {
        HostSystem moose;
        if (readPetscSystem(ref_mat, ref_vec, m, moose))
          std::printf("# V0 vs MOOSE Kokkos assembly: max rel diff = %.3e\n",
                      maxRelDiff(moose, ref));
        else
          std::printf("# could not read MOOSE reference\n");
      }

      std::printf("variant,max_rel_diff_vs_V0,distinct_hashes,face_ms,boundary_ms,cell_ms,total_ms,"
                  "ns_per_cell\n");
      size_t pos = 0;
      while (pos < variants.size())
      {
        const auto comma = variants.find(',', pos);
        const auto v =
            variants.substr(pos, comma == std::string::npos ? std::string::npos : comma - pos);
        pos = comma == std::string::npos ? variants.size() : comma + 1;

        std::vector<double> warm(3, 0.0);
        run(v, warm);
        run(v, warm);
        const double diff = maxRelDiff(ref, toHost(D, H));
        std::vector<std::vector<double>> times;
        std::set<uint64_t> hashes;
        for (int r = 0; r < reps; ++r)
        {
          std::vector<double> t(3, 0.0);
          run(v, t);
          times.push_back(t);
          hashes.insert(hashSystem(D));
        }
        auto med = [&](int j)
        {
          std::vector<double> x;
          for (auto & t : times)
            x.push_back(j < 0 ? t[0] + t[1] + t[2] : t[j]);
          std::nth_element(x.begin(), x.begin() + x.size() / 2, x.end());
          return x[x.size() / 2];
        };
        const double tot = med(-1);
        std::printf("%s,%.3e,%zu,%.3f,%.3f,%.3f,%.3f,%.2f\n",
                    v.c_str(),
                    diff,
                    hashes.size(),
                    med(0) * 1e3,
                    med(1) * 1e3,
                    med(2) * 1e3,
                    tot * 1e3,
                    tot / m.nc * 1e9);
        if (diff > 1e-11)
          rc = 1;
      }
    }
  }
  Kokkos::finalize();
  return rc;
}
