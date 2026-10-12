# Kokkos-MOOSE: Design Choices

Survey of the Kokkos work in MOOSE as of `next` at `410cc90c7e` (2026-10). Paths are relative to
the `moose/` repository root; line numbers are approximate. Items marked *(inferred)* are readings
of the code rather than statements made in it.

Primary user documentation: `framework/doc/content/syntax/Kokkos/index.md` and the per-system
pages under `framework/doc/content/syntax/Kokkos*/`.

---

## 1. Overall philosophy

1. **Port MOOSE's object systems, keep libMesh as the FE backend.** Kokkos-MOOSE is not a new
   finite element backend (unlike MFEM-MOOSE). It keeps `MooseMesh`/libMesh `MeshBase`, `DofMap`,
   and PETSc, and provides device-portable versions of selected object systems
   (`modules/doc/content/framework_development/computational_backends.md`).
2. **Retain the MOOSE user experience.** Ported objects mirror the CPU API, parameter names, and
   doc pages; concrete objects are named `Kokkos<CPUName>` (`KokkosDiffusion`, `KokkosDirichletBC`,
   `KokkosHeatConduction`), while framework bases live in `namespace Moose::Kokkos` without the
   prefix (`Moose::Kokkos::Kernel`, `IntegratedBC`, `Material`, ...).
3. **Opt-in and mixable with CPU objects.** Kokkos and CPU objects can appear in the same input
   and assemble into the same PETSc vectors/matrices (test `test/tests/kokkos/fe-fv-interplay`).
4. **Portability through Kokkos, with Kokkos taken from PETSc.** Kokkos is never vendored; the
   build reads `petscconf.h` and inherits PETSc's backend (CUDA, HIP, SYCL, or OpenMP on CPU) and
   GPU architecture (`framework/kokkos.mk`). `--with-kokkos=cpu` runs the same code through the
   OpenMP backend for coverage and debugging.

## 2. Build system and compilation model

1. **Separate compilation via `.K` files.** Only `*.K` sources may contain Kokkos code. They are
   compiled by the GPU compiler (`nvcc --extended-lambda -x cu`, `hipcc`, `icpx -fsycl
   -fno-sycl-rdc`, or the host compiler) into separate shared libraries
   (`libmoose_kokkos-<method>.so`, `lib<app>_kokkos-<method>.so`), then linked into the CPU-built
   executable. An NVCC device-link step produces `lib/dlink.o` (`framework/app.mk`,
   `framework/moose.mk`). `.K` files are unity-built per source directory.
2. **Two preprocessor scopes.** `MOOSE_KOKKOS_ENABLED` (configure-time, visible everywhere) guards
   Kokkos-*related* code in ordinary `.C`/`.h` files; `MOOSE_KOKKOS_SCOPE` (only in `.K`
   compilation) guards Kokkos code in shared headers. `KokkosHeader.h` `#error`s outside
   `MOOSE_KOKKOS_SCOPE`.
3. **Class layouts must be identical across compilers.** `MOOSE_KOKKOS_SCOPE` must never guard
   member variables or virtual functions, otherwise the CPU and GPU compilers see different object
   sizes (`syntax/Kokkos/index.md:442`). Consequently data members are guarded only by
   `MOOSE_KOKKOS_ENABLED`.
4. **One class split across two compilers.** Kokkos member functions of `FEProblemBase`,
   `NonlinearSystemBase`, `LinearSystem`, `AuxiliarySystem` are declared in the normal headers and
   defined in `.K` files (`framework/src/kokkos/problems/KokkosFEProblemBase.K`,
   `framework/src/kokkos/systems/KokkosNonlinearSystemBase.K`, ...).
5. **No relocatable device code (RDC).** Disabled because of PETSc/Kokkos limitations, so every
   device function must be inlineable (defined in a header or the same TU). Enabling RDC is listed
   as high priority; this is the main reason so much code lives in headers.
6. **Configurable index width.** `--with-kokkos-index-size=4|8` sets `MOOSE_KOKKOS_INDEX_TYPE`
   (`uint32_t` default), used for thread IDs and the default `Array` index type. No mixed/single
   precision option: `Real` is MOOSE's `Real`.

## 3. Static polymorphism without virtual functions

1. **No device virtuals; no CRTP class templates.** Virtuals require constructing objects on the
   device to populate vtables, are not supported by every backend, and block inlining. CRTP was
   rejected to spare users class templates (`syntax/Kokkos/index.md`, "Static Polymorphism").
2. **Templated hook methods.** Base classes are ordinary classes; hooks are
   `template <typename Derived> KOKKOS_FUNCTION ... const` member functions (e.g.
   `computeQpResidual<Derived>(i, qp, datum)`). Users "override" by *hiding*; hooks must be public.
   This replaced an earlier "shim method" design (`computeQpResidualShim`, removed in #32512).
3. **Dispatcher registry.** `registerKokkos*` macros (`framework/include/kokkos/base/KokkosDispatcher.h`)
   call `registerMooseObject` and, through a static initializer, register a type-erased
   dispatcher builder in a singleton `DispatcherRegistry` keyed by
   `(typeid(OperationTag), registered type name)`. At run time the object's `type()` string finds
   the builder, which captured the concrete C++ type at registration. The only virtual layer is
   host-side (`DispatcherBase::parallelFor/parallelReduce`).
4. **Operation tags.** Loops are selected by empty tag structs (`ResidualLoop`, `JacobianLoop`,
   `OffDiagJacobianLoop`, `OffDiagJacobianScalarLoop`, `PresetLoop`, `ElementLoop`, `NodeLoop`,
   ...). The functor's `operator()(Tag, tid, const Derived & self)` receives the fully-typed object
   back, so base-class template loops can call `self.template computeQpResidual<Derived>(...)`.
5. **Detecting user-defined hooks at registration.** The registry compares member-function
   pointers against the base default (`&Object::template computeQpJacobian<Object> !=
   Object::template defaultJacobian<Object>()`). Unimplemented Jacobian/off-diagonal loops are never
   launched; default hooks `Kokkos::abort` if reached. Caveat: a misspelled hook silently falls back
   to the default.
6. **Compile-time traits** select which dispatchers get registered: `uses_precompute_hooks`,
   `supports_scalar_jacobian`, and SFINAE `has_split_linear_fv_flux_dispatchers`.
7. **Exception: Functions use real device virtuals.** `Moose::Kokkos::Function` is a type-erased
   wrapper whose device object is placement-new'd inside a one-thread `parallel_for` so the vtable
   holds device pointers (`KokkosFunctorWrapper.h`, `KokkosFunctionWrapper.h`). Because this needs
   RDC, retrieving a function as the abstract type is an error on GPU; GPU code must use
   `getKokkosFunction<ConcreteType>()`.

## 4. Functor copy semantics and host/device data binding

1. **Every Kokkos object is a functor, copied to the device at each dispatch.** A single `const`
   copy is shared by all threads, so members must not hold thread-private state and device-callable
   methods must be `const`.
2. **The copy constructor is the synchronization hook.** Every MOOSE base/interface gained a
   special copy constructor taking `Moose::Kokkos::FunctorCopy`. `Dispatcher` re-copies from the host
   reference on each launch, so each dispatch snapshots current host state (e.g. `ResidualObject`
   rebuilds its tag arrays from `TaggingInterface`; `KokkosResidualObject.K:62-85`).
3. **Value binding instead of reference binding.** Host references are useless on device;
   `ReferenceWrapper<T>` holds a host reference plus a value copy refreshed at dispatch, and
   `Scalar<T>` adds arithmetic. Used for `_t`, `_dt`, postprocessor values
   (`Moose::Kokkos::PostprocessorValue`), and restartable data.
4. **Holder pattern.** `MeshHolder`, `AssemblyHolder`, `SystemHolder`/`FESystemHolder`,
   `MemoryPoolHolder` follow the same "host reference + refreshed device copy" idiom; accessors pick
   one with `KOKKOS_IF_ON_HOST`. On the CPU backend "device" means "inside an OpenMP parallel
   region" (`KOKKOS_IF_ON_HOST` is redefined as `if (!omp_get_level())`, `KokkosHeader.h:54-66`).
5. **Dispatches are synchronous.** Every `parallel_for`/`parallel_reduce` is followed by
   `Kokkos::fence()`.

## 5. Device containers and memory

1. **Custom containers instead of `Kokkos::View`.** `Moose::Kokkos::Array` (up to N-D via variadic
   templates, LEFT layout by default, configurable index type), `JaggedArray` (flat storage + dope
   vectors to avoid nested-array fragmentation), and `Map` (host `std::map` + device FNV-1a hash
   table, populated on host only, not performance-tuned). Device memory comes from `kokkos_malloc`.
2. **Explicit dual host/device storage.** `create()`, `createHost()`, `createDevice()`, plus
   `aliasHost()`/`aliasDevice()` to wrap external pointers (e.g. PETSc arrays). No automatic
   synchronization: explicit `copyToHost()`, `copyToDevice()`, `copyToDeviceNested()`. Device memory
   is uninitialized on allocation.
3. **Shallow copy by default, reference counted.** Copy is shallow (shared `_counter`; last owner
   frees) to make functor dispatch cheap. `ArrayDeepCopy<T>` opts a type into per-entry copy
   construction (used by `MaterialProperty` so clones resync).
4. **Non-default-constructible entries.** `create<false>()` mallocs raw storage and `emplace()`
   constructs entries, with sparse construction tracked.
5. **No dynamic allocation in device code.** Fixed-size local arrays are preferred (e.g.
   `MAX_CACHED_DOF = 30` stack buffers). An opt-in N-way `MemoryPool`
   (`MooseApp::allocateKokkosMemoryPool(size, ways)`, RAII `MemoryChunk`) exists; no framework code
   allocates one.
6. **Compile-time capacity limits.** `MAX_TAG = 30` tag slots, `MAX_CACHED_DOF = 30`,
   `MOOSE_AD_MAX_DOFS_PER_ELEM` for AD, parsed-expression stack depth 10, VPP per-element cache 10
   (`KokkosHeader.h:81-85`).

## 6. Execution model and parallel decomposition

1. **Flat `RangePolicy` everywhere.** No `TeamPolicy`, `MDRangePolicy`, or scratch memory.
   `Thread<>` decomposes a 1-D index into up to 4 strided dimensions (dimension 0 fastest).
2. **Thread mapping.**
   - Kernels: `(local_thread, element)`; off-diagonal adds a `jvar` dimension; scalar off-diagonal
     is `(scalar jvar, element)`.
   - Integrated BCs: boundary element-sides. Nodal BCs/NodalKernels: nodes. AuxKernels: elements or
     nodes. Materials: elements/element-sides. The qp loop is serial within a thread.
3. **Intra-element "local DOF" parallelism.** `num_local_residual_threads` (default 2) and
   `num_local_jacobian_threads` (default 4) split the test-function loop of one element across
   adjacent threads; defaults were "chosen empirically for first-order Lagrange shape functions in
   3D" (`KokkosLocalParallelInterface`). Forced to 1 on the CPU backend.
4. **Batched local accumulation.** Residuals use a stack `local_re[MAX_CACHED_DOF]`, processing
   DOFs in batches across local threads; Jacobians are column-parallel (threads own `j`, batch `i`).
5. **Single MOOSE thread.** Kokkos objects exist only on thread 0 (`n_threads = 1` for warehouses
   and materials); on-node parallelism comes from Kokkos, not MOOSE threads.

## 7. Mesh, assembly, and FE data on device

1. **Contiguous renumbering.** `Moose::Kokkos::Mesh` assigns contiguous per-process element and node
   IDs, ordered by subdomain so each subdomain is a contiguous range; stored as libMesh extra
   integers (`kokkos_contiguous_elem_id`, `kokkos_contiguous_node_id`). `(ElemType, p_level)` pairs
   are enumerated into element-type IDs. Rebuilt lazily on `meshChanged`.
2. **Ghost layer only when needed.** One layer of off-process neighbors is numbered after local
   elements only if side geometry is required (linear FV).
3. **Hybrid precompute / on-the-fly assembly** (`KokkosAssembly`):
   - Precomputed per reference element, indexed by (subdomain, element type, FE type): quadrature
     points/weights, `phi`, `grad_phi`, vector shape functions, and geometric map shape functions,
     obtained by reinitializing libMesh `FE` on a reference element.
   - Cached per physical volume qp: inverse Jacobian, `JxW`, `xyz` (hence no displaced mesh).
   - Computed on the fly for faces: Jacobian, `JxW`, `xyz`, normals.
   - Face qp storage is allocated only for boundaries that consume boundary material properties.
4. **Variable values: precompute volume, compute faces on demand.** `FESystem::reinit` dispatches
   over `(qp, elem, var, tag)` to fill `_qp_solutions` (values, gradients, vector values, curls) for
   active variables/tags only; face values are interpolated per call. Volume terms trade memory
   for speed.
5. **Coordinate systems** (XYZ, RZ with general axis, RSPHERICAL) are folded into `JxW` per
   subdomain (`Assembly::coordTransformFactor`).
6. **Datum objects** are thread-private stack temporaries carrying the element/side/node context:
   `MeshDatum`, `FVDatum` (face centroid, normal, area, volume, `faceDCNMag`), `Datum` (assembly,
   qp offsets, one-qp transform cache), `AssemblyDatum` (`ivar`, `jvar`, DOF counts,
   `do_derivatives`). They replace `_current_elem`, `_qp`, `_i`, `_j`, which cannot be members.
7. **Supported discretizations are deliberately narrow.** FE families LAGRANGE, L2_LAGRANGE,
   MONOMIAL, SCALAR, and LAGRANGE_VEC only. Rejected: displaced mesh, adaptivity, side with
   multiple boundary IDs, non-contiguous nodal DOFs.

## 8. Solution vectors, matrices, and PETSc

1. **PETSc is the memory boundary.** `Moose::Kokkos::Vector` uses `VecGhostGetLocalForm` +
   `VecGetArrayAndMemType`: zero-copy aliasing of device-resident PETSc vectors
   (`-vec_type kokkos`), otherwise a device mirror aliased to host memory with explicit copies.
   Device vec/mat types are not selected automatically; the user sets PETSc options.
2. **Old/older/dot solutions are tags.** `Variable` is (system, variable, tag); time history and
   time derivatives are just different vector tags. Uncoupled variables carry default values.
3. **Per-evaluation sync of only what is needed.** Each residual/Jacobian/UO pass resolves
   variable, tag, and material-property dependencies, then syncs host→device only the active
   variables/tags, and syncs results back afterward.
4. **Off-process residual contributions go through PETSc host assembly.** A separate `_ghost` array
   collects off-process sums, closed with `VecSetValues(..., ADD_VALUES)`. A device-side
   `MPI_Alltoallv` path exists but is compiled out (`// #define DEVICE_ASSEMBLY`,
   `KokkosVector.K:13`).
5. **Matrices: CSR on device, COO into PETSc.** `Moose::Kokkos::Matrix` builds CSR from libMesh's
   *full* sparsity pattern (forced whenever Kokkos residual objects exist,
   `NonlinearSystemBase.C:236-239`), including ghost rows fetched from owners. Assembly uses
   `MatSetPreallocationCOO` once and `MatSetValuesCOO(..., ADD_VALUES)` with host or device values
   depending on PETSc memory type. Device entry lookup is a binary search within the row.
6. **Scalable ghost setup.** Communication lists use sparse `TIMPI::push_parallel_vector_data`
   (no O(p) loops over ranks, #33937).
7. **System allocation by need.** Plain `System` for FV-only/scalar-only, `FESystem` for qp-based,
   and both sharing device memory for mixed FE/FV systems.

## 9. Residual and Jacobian assembly

1. **Atomics, not coloring.** Element contributions use `Kokkos::atomic_add` into every active
   tagged vector/matrix; nodal objects (one thread per node) write directly; linear FV chooses
   `AccumulationMode::{Atomic, NonAtomic}`. There is no mesh coloring.
2. **Kokkos runs first, then CPU.** `computeKokkosResidual`/`computeKokkosJacobian`/
   `computeKokkosResidualAndJacobian` run before the host loops in `NonlinearSystemBase`; nodal BCs
   run in a separate pass. Each stage is timed separately (KokkosCopy, KokkosReinit,
   KokkosMaterial, KokkosKernel, KokkosClose).
3. **Nodal BC race avoidance.** `FESystem::setupNodalBCDofs` builds a per-matrix-tag mask of
   nodal-BC rows; element accumulation skips those rows, and the nodal BC zeros the row and sets
   the diagonal. Nodal BC residuals are *set*, not added.
4. **Hand-coded Jacobians with full coupling.** Off-diagonal `jvar` lists come from the problem's
   coupling matrix (`SMP full=true` gives full blocks); field-to-scalar off-diagonals supported via
   `computeQpOffDiagJacobianScalar`.
5. **Device AD via MetaPhysicL.** `Moose::Kokkos::ADReal = DualNumber<Real,
   KokkosSemiDynamicSparseNumberArray<..., MOOSE_AD_MAX_DOFS_PER_ELEM>>`, derivatives keyed by
   *global* DOF index. Seeds (1, `duDotDu` for dot, 0 for old) are set at dispatch; qp values are
   re-interpolated from DOFs when seeding. AD objects register only `ResidualLoop`, which serves
   residual, Jacobian, or both; derivatives scatter straight into CSR rows. AD exists for scalar
   Kernel, TimeKernel, IntegratedBC, NodalBC only (no AD vector kernels, nodal kernels, or AD
   material properties).
6. **Precompute variants.** `KernelValue`/`KernelGrad`/`IntegratedBCValue` (+ vector versions)
   factor out the test function via `precomputeQpResidual/Jacobian` for performance.
7. **Overridable inner loops.** Derived classes may hide `computeResidualInternal`/
   `computeJacobianInternal` (e.g. `KokkosTimeDerivative` lumping) or AuxKernel
   `computeElementInternal`.

## 10. Object-system specifics

**Kernels / BCs / NodalKernels**
- Hierarchy rooted at `Moose::Kokkos::ResidualObject` (derives from CPU `::ResidualObject` plus the
  holders). `computeOffDiagJacobian(jvar)` is `final` and errors.
- Variable access through callable objects: `_u(datum, qp)`, `_grad_test(datum, i, qp)`,
  coupling through `kokkosCoupled*` APIs (value, gradient, old/older, dot, nodal, DOF, vector tag,
  vector FE, AD variants). No second derivatives.
- Nodal hooks keep a `qp` argument (always 0) for interface uniformity.
- Preset Dirichlet BCs use a `PresetLoop` writing into a `parallel_solution` tag.

**AuxKernels**
- Nodal: one thread per node, direct write. Elemental: per-element local L2 projection solved with
  `Utils::choleskySolve`. Kokkos aux kernels run after all CPU aux kernels. No boundary-restricted
  elemental aux kernels; no nodal aux kernels using material properties.

**Linear FV** (#32492, July 2026; stated as groundwork for application GPU acceleration such as
OpenPronghorn)
- Separate `LinearSystemContributionObject` base with matrix/RHS hooks on an `FVDatum`; designed to
  avoid FE assembly machinery.
- Flux kernels loop over *directed* element sides: each interior face is visited from both sides
  and each visit writes only its own row. Split loops `Internal/Boundary x Matrix/RightHandSide`.
- BCs precompute per-face `BoundaryRelation{coefficient, source}`, read by flux kernels. Mixing
  non-Kokkos linear FV BCs with Kokkos flux kernels is an error.
- Current objects: `KokkosLinearFVDiffusion` (two-point flux, no non-orthogonal correction
  *(inferred)*), `KokkosLinearFVAdvection` (constant velocity, first-order upwind),
  `KokkosLinearFVSource`, functor Dirichlet/Neumann BCs. Not threaded.

**Materials**
- Per-property storage `Array<Array<T, dim+1>>`: outer per subdomain, inner device-only with a
  trailing qp/element/1 index depending on `constant_on`. Rank is a template parameter; sizes are
  runtime and may vary per subdomain.
- All qps of all elements are stored (GPU parallelism prevents per-element storage); materials are
  computed eagerly, one `parallel_for` per material over all its elements, in the order block ->
  face -> neighbor -> boundary-restricted.
- Properties are returned *by value* as shallow clones that resync from a reference property at
  each dispatch (survives reallocation and stateful swaps).
- Stateful: every state is a full device array; `shift()` is an O(1) array swap. Old/older requests
  do not create dependencies.
- Dependency sorting uses the normal resolver but across all blocks at once ("Kokkos computes all
  elements and faces at once regardless of block and boundary").
- On-demand properties allocated only if requested. `constant_on` = NONE/ELEMENT/SUBDOMAIN.
- **Strict CPU/Kokkos separation**: CPU objects cannot read Kokkos properties and vice versa; the
  bridge is `KokkosMaterialRealAux`. Trivially copyable types recommended, not enforced (stateful
  types need `dataStore`/`dataLoad`).

**Functions**
- Concrete-type retrieval (`getKokkosFunction<T>`) is the GPU path; abstract retrieval is CPU-only
  (see 3.7). Numeric defaults create `KokkosConstantFunction`.
- **Parsed expressions without JIT.** PEG parser -> AST -> RPN on host; a device stack machine
  (`RPNEvaluator`, stack depth 10) interprets it. Symbols: constants, postprocessors, variables,
  material properties, functions.

**UserObjects / Postprocessors / VectorPostprocessors / Reporters**
- Two modes chosen automatically: `execute<Derived>` -> `parallel_for`, or
  `reduce<Derived>` + `join` + `init` -> `parallel_reduce` with a runtime-length `Real[]` reduction
  landing directly in a host `View`. All reduced values are `Real`. `join` replaces `threadJoin`;
  MPI reduction happens in `finalize()`. This API intentionally differs from CPU postprocessors.
- VPPs may choose `calculation_mode = atomic|reduction`; reducers suit small dense buffers, atomics
  suit large or sparse ones.
- Results enter the same `ReporterData` as CPU objects.
- Within an execution group Kokkos UOs run before CPU UOs; cross-dependencies between them require
  `execution_order_group` (automatic resolution is planned).

## 11. Integration with the rest of MOOSE

1. **Same input blocks, flag-based routing.** Kokkos objects go in the normal `[Kernels]`, `[BCs]`,
   `[Materials]`, ... blocks. A private parameter `_kokkos_object` (`isKokkosObject()`) makes the
   existing `Add*Action`s route to `FEProblemBase::addKokkos*`. This *reversed* the original design
   of separate `[Kokkos*]` blocks and `AddKokkos*Action`s (`6054977ec3`, Nov 2025).
2. **Separate warehouses.** `_kokkos_materials`, `_kokkos_functions`, Kokkos kernel/BC warehouses,
   and a `TheWarehouse` attribute `AttribKokkos` so host loops filter Kokkos objects out.
3. **Base names reuse CPU ones** (`registerBase("Kernel")`), so syntax/documentation tooling treats
   them as normal objects.
4. **Initialization.** `FEProblemBase::initKokkos()` (only when Kokkos objects exist) builds the
   Kokkos mesh, assembly, systems, and material storages; it rejects displaced meshes and
   adaptivity.
5. **GPU selection per rank.** `device_id = local_rank % num_gpus` on a node-local communicator,
   using the minimum GPU count across ranks; multiple ranks per GPU expected. For CUDA the device
   is also set *before* `MPI_Init` from launcher environment variables, to avoid >100x slowdowns
   with Open MPI smcuda/NVLink (`MooseInit.C`, #33585). `Kokkos::finalize` runs at PETSc finalize.
6. **Restart.** Kokkos material storages are restartable data with custom `dataStore`/`dataLoad`
   (per-property blobs, type-erased per-type store functions); `Array` and `Map` serialize with
   type/rank/extent checks. Objects can declare restartable `Moose::Kokkos::Array` data.
7. **Automatic scaling** works by assembling into the real `PetscMatrix` and using
   `MatGetDiagonal`/`MatGetRowSumAbs` (#33867).
8. **Testing.** Specs use `capabilities = 'kokkos'` and `compute_devices = 'cpu cuda xpu'`; the
   Kokkos backend is fixed at build time, so CPU builds exercise the same tests via OpenMP.

## 12. Evolution (key milestones)

| Date | Change |
|---|---|
| 2025-06 | Initial PR #30835 (originally named "GPU", renamed Kokkos with namespace). |
| 2025-07 | Configure option; CPU-only mode; Kokkos added to the CPU PETSc build. |
| 2025-09 | Merged: Kernels, NodalKernels, BCs, Materials (NVIDIA only). |
| 2025-09/10 | AuxKernels, JaggedArray, dispatcher registry. |
| 2025-11 | Functions; separate Kokkos syntax merged into MOOSE syntax. |
| 2026-01/02 | UserObjects/Postprocessors on Kokkos reducers; on-demand properties; `constant_on`. |
| 2026-03/04 | Shim methods replaced by templated hooks; variadic Array/MaterialProperty dims; parsed objects; VPPs and Reporters; local DOF parallelism; first module object. |
| 2026-04/05 | Device AD (MetaPhysicL); unified `registerKokkosResidualObject`. |
| 2026-06/07 | Linear FV; mixed FE/FV systems; HIP and SYCL builds; vector FE. |
| 2026-09 | FV variables coupleable in FE objects; scalar variables; pre-p-refined meshes; automatic scaling. |

Module usage today: `heat_transfer` (conduction kernels, material, convective BCs) and
`thermal_hydraulics` (radiative heat flux BC). **`navier_stokes` has no Kokkos code.**

## 13. Gaps relevant to a Navier-Stokes port *(inferred from the above)*

- No Kokkos nonlinear FV (`FVKernels`), DG/interface kernels, constraints/mortar, or array
  variables; linear FV is limited to diffusion, constant-velocity first-order upwind advection, and
  source, without non-orthogonal correction or gradient reconstruction.
- No AD material properties; AD limited to scalar kernels/BCs.
- No interior-face or functor materials; CPU functor materials cannot feed Kokkos objects.
- No displaced mesh or adaptivity.
- Off-process residual assembly and many sync points still go through host PETSc.

## 14. Possible defects noticed during the survey

- `framework/src/multiapps/FullSolveMultiApp.C:113` uses `#ifdef KOKKOS_ENABLED`, but
  `AX_PREFIX_CONFIG_H` (`configure.ac:28`) renames the macro to `MOOSE_KOKKOS_ENABLED`, so the
  Kokkos stateful-property check there is likely dead code (not confirmed with a build).
- `MaterialData::getKokkosProperty` (`framework/include/materials/MaterialData.h:466-471`) loops
  `for s <= state` but passes `state` instead of `s` to `addKokkosPropertyHelper`, so intermediate
  state slots are not reserved by that loop. Requesting only an `Older` property may break `shift()`;
  existing tests always request `Old` alongside `Older`.
- `NodalReducer::_unique_node_execute` is read from parameters but never used
  (`framework/include/kokkos/userobjects/KokkosNodalReducer.h:52`).
