# Cleanup of `vtkAMRAmazeReaderInternal.h` / `.cxx` and `vtkAMRAmazeReader.h` / `.cxx`

Cleanup of the low-level Amaze reader internals in
`ParaViewPlugin/src/Reader/vtkAMRAmazeReaderInternal.{h,cxx}` and of the
public reader `vtkAMRAmazeReader.{h,cxx}`.

Net result: **−327 lines** (1198 insertions, 1525 deletions), clean build with
zero compiler warnings, and a latent link bug found and fixed along the way.

## The big find: a lost function

Commit `5aeb71d` ("upgrade to use vtkAMRBaseReader") accidentally deleted the
definition of `vtkAMRAmazeReaderInternal::AxisSymStarSource()` while leaving its
declaration and its call in `BuildStars()`. The shared library still linked only
because Linux `.so` files allow undefined symbols — `nm` confirmed
`U vtkAMRAmazeReaderInternal::AxisSymStarSource(...)`, i.e. a `dlopen` crash
waiting to happen for any file containing axisymmetric stars.

The implementation was restored from git history (before `5aeb71d`), cleaned up,
and is now a defined symbol (`T`) in `libAMAZEReader.so`.

## Header (`vtkAMRAmazeReaderInternal.h`)

- Removed declarations that had no definition anywhere: `PrintSelf`,
  `ReadRectilinearGrid`, `ReadStructuredGrid`, `ReadStructuredGrid2` (only the
  old VisIt plugin needed those; it carries its own copy), plus the unused
  `GetStar` and the `FILE *errs` member.
- Removed unused constants (`adG_NAMELENGTH`, `adG_PATHLENGTH`,
  `adG_BINARY`/`adG_ASCII` macros) and unused includes.
- Modernized C-style `typedef struct X {...} X;` to plain `struct X {...};`.
- Added deleted copy constructor / assignment operator.
- `FileName` is now a `std::string` (was a non-owning raw `char*`).
- `newstars` is now a `std::vector<newstar>` (was a leaked `new[]` array).
- Members now have default initializers (`NumberOfLevels`, `Dimensionality`,
  `AMAZETime`, `AMAZETimeScalor`, `DataScale`, …) — previously several were
  uninitialized garbage if the corresponding HDF5 attributes were missing.
- HDF5 internals (`ReadHDF5MetaData`, `ReadHDF5VariablesMetaData`,
  `ReadHDF5GridsMetaData`, `ReadVar`, `FindLevelAndBlock`, `FindDomainId`,
  `AxisSymStarSource`) moved to `private`; the public API is now exactly what
  `vtkAMRAmazeReader.cxx` uses.
- Dropped the never-implemented `SetFileName` ownership ambiguity: it now simply
  copies the name (and closes any open file when set to null).

## Implementation (`vtkAMRAmazeReaderInternal.cxx`)

- Dropped ~250 lines of commented-out debug blocks, per-star `cerr`/`cout`
  spam, and dead code: `CheckVarSize` (computed values and did nothing with
  them), the `if(0) //GetCellCentered()` branch, the `LatLongTessellation`
  dead branch in `AxisSymStarSource`, the no-op `stars[i].Radius *= 1;`, and
  the empty axis-symmetric `while` loop.
- Replaced manual `new[]/delete[]` (`interactions`, `models`, `stars`,
  `spherStarData`, `axiStarData`) with `std::vector` / stack objects.
- Collapsed duplicated code: the "Cell Edge Length Y/Z" attribute reads, the
  "Polytropic EOS" vs "Ideal Gas EOS" branches, and the repeated HDF5
  string/array compound-type boilerplate (now `CreateStringType()` /
  `CreateDoubleArrayType()` helpers).
- Moved file-local structs and the attribute-iterator callback into an
  anonymous namespace / `static` linkage.
- Fixed real memory leaks:
  - The `vtkSphereSource` and `vtkTransformFilter` per-star leaks (star outputs
    are now `Register(nullptr)`-ed and the filters deleted). This also fixes the
    leak on every `RequestData`, since `BuildStars()` rebuilt stars on each
    update without deleting (or even releasing) the previous ones.
- Fixed error paths in `BuildStars()`: the early `return 0` (no APR section)
  previously left the HDF5 file handle open and HDF5 error reporting disabled.
- Fixed the `AxiSymStarCurrent` compound type being closed *inside* the
  axis-symmetric star loop (invalid handle for the 2nd star); it is now closed
  once after the loop.
- Small robustness fixes without behavior change: guards on failed `H5Dopen`
  of "Grid Info", null check before `*strchr(...) = '\0'`, and initialization of
  `nb_stars` / `Radius` which were previously returned/read uninitialized on
  some paths.
- Modernized VTK calls for VTK 9 (`vtkPoints::Reserve`, `vtkCellArray`
  initializer-list `InsertNextCell`, `vtkNew` for transforms/filters).

## Behavior preserved / notes

- `MakeVariableNames`, `GetSpacing`, `GetAMRGrid`, `ReadVar` and the HDF5
  metadata reads produce identical results as before; the reader
  (`vtkAMRAmazeReader.*`) was not modified and compiles against the new header.
- `CheckVarSize` was removed entirely: it opened datasets, read their extents,
  compared them to the grid size, and did nothing either way.
- `ReadHDF5GridsMetaData(bool shiftedGrid)` lost its parameter: the "shifted
  grid" branch was never called with `true` anywhere (the reader's
  `ShiftedGrid` property was already a no-op) and the branch also closed the
  shared file handle in surprising ways.
- **Follow-up suggestion (not done here):** `vtkAMRAmazeReader`'s
  `LogData` / `ScaleChoice` / `DataScale` / `ShiftedGrid` properties are (and
  always were) never forwarded to the `Internal` object, so the UI checkboxes
  mostly have no effect. Worth a look separately.
- **Validation done:** full `make` of the plugin in `buildDev` succeeds with no
  warnings; `nm` shows no remaining undefined `vtkAMRAmaze*` symbols.
  **Validation not done:** no runtime smoke test, since no A-MAZE dataset path
  is configured in `Testing/Python`. Recommended: open a dataset with stars
  (ideally one with an axisymmetric star) in ParaView before committing.

---

# Cleanup of the public reader (`vtkAMRAmazeReader.h` / `.cxx`)

Same treatment applied to the reader class itself. All properties bound by
`AmazeReader_Server.xml` (`SetLogData`, `SetShiftedGrid`, `SetLengthScale`,
`SetScaleChoice`, `SetEnableCaching`, `SetMaxLevel`, point-array selections,
`GetNumberOfLevels/Blocks`) are untouched, so the ParaView proxy keeps working.

## Header (`vtkAMRAmazeReader.h`)

- All members now have default initializers matching the proxy defaults
  (`LogData{0}`, `LengthScale{1}`, `ShiftedGrid{0}`, `DataScale{1}`,
  `ScaleChoice{3}` = "No-scale", `IsReady{false}`). Previously `LengthScale`,
  `ScaleChoice` and `DataScale` were uninitialized until the UI pushed values
  (the ctor's `DataScaleOn()` call is replaced by the initializer).
- Renamed the cryptic `nbstars` to `NumberOfStars` and documented `ScaleChoice`'s
  default.
- Removed `GetStarsOutput()`: it was unused everywhere, and it was also broken —
  it looked for output port 2 while the reader only has 2 ports (indices 0–1).
- Dropped now-unneeded forward declarations (`vtkOverlappingAMR`).

## Implementation (`vtkAMRAmazeReader.cxx`)

- Removed the commented-out dead code (old `CanReadFile` implementation,
  array-selection hacks, debug `cerr` lines) and the misleading duplicate
  comments.
- Destructor body removed: the base `vtkAMRBaseReader` destructor already frees
  `FileName` and the `Internal` `unique_ptr` cleans itself up; the destructor is
  now `= default`.
- `RequestData`:
  - removed the `if(1) // only MPI task 0 ...` block and its unreachable `else`
    branch (the code was unconditionally executing the "task 0" path);
  - removed the unused `amr->CheckValidity()` call and unused temporaries;
  - **fixed a potential out-of-bounds read**: the per-star `FieldData` loop
    indexed `Internal->newstars[s]` for `s < nbstars` — with old-style
    (pre-2008) star files `newstars` is empty while `nbstars > 0`. The loop is
    now bounded by `newstars.size()`;
  - clearer star-count mismatch message.
- `FillMetaData`: inverted the big `if (nbstars == -1) { ... } else { return 1; }`
  wrapper into an early return; de-nested and clarified the refinement-ratio
  comment.
- `GetAMRGridPointData`: added a null guard on `ReadVisItVar()` — it returns
  `nullptr` for unknown variable names, and the old code passed that straight
  to `AddArray()` (crash).
- `SetFileName`: same logic, simplified branch conditions; still calls
  `Internal->SetFileName(nullptr)` first so the internal reader releases the
  stale path before the buffer is reallocated.
- `CanReadFile`: null-pointer guard (`fname && H5Fis_hdf5(fname) > 0`).
- Includes cleaned: removed unused `<cfloat>` / `<map>`, added explicit
  `<algorithm>`, `<cstring>` and the VTK headers the code actually uses
  (`vtkFieldData.h`, `vtkPolyData.h`, `vtkUniformGrid.h`) instead of relying on
  transitive includes.

## Follow-up suggestion (partially resolved)

`vtkAMRAmazeReader`'s `ScaleChoice` / `ShiftedGrid` / `DataScale` properties
are (and always were) never forwarded to the `Internal` object, so the
corresponding UI checkboxes have no effect on reading. The defaults were kept
consistent with the proxy XML, but making them functional would require
forwarding them in `SetFileName`/`FillMetaData`. (`LogData` was in the same
situation and has since been removed entirely — see the next section.)

---

# Removal of the dead `LogData` feature (2025-11)

`LogData` ("Compute Log10" checkbox: log10 of Density/Temperature/Pressure)
was removed from the whole stack. It was **dead in practice**: the public
reader's `LogData` member was never forwarded to `vtkAMRAmazeReaderInternal`
(the internal member defaulted to 0 and had no public setter path), so the
UI checkbox never changed any output. Removing it therefore cannot change any
result — confirmed by rerunning `TestAMAZEFullFile.py` after the change
(all values identical: 6 levels, 384 blocks, 19,805,712 cells, Density range
`[4.655599e-15, 4.926240e-15]`, 2 stars).

Removed:

- `AmazeReader_Server.xml`: the `LogData` `IntVectorProperty` and its entry in
  the core proxy's `ExposedProperties`.
- `vtkAMRAmazeReader.h`: `vtkSet/Get/BooleanMacro(LogData)` and the member.
- `vtkAMRAmazeReaderInternal.h`: `LogData` member and the `VarNamesToLog` map.
- `vtkAMRAmazeReaderInternal.cxx`:
  - the `VarNamesToLog` seeding in the constructor (now `= default`);
  - the `Log10({}) [{}]` branch of `MakeVariableNames()` (names are simply
    `label [unit]` now);
  - the log10 pass over `dArray` in `ReadVar()`;
  - the `LogData ? log10(...) : ...` ternary on star temperature in
    `BuildStars()`;
  - the log10 post-pass on `Temperature` in `AxisSymStarSource()`.
- `src/Testing/Cxx/TestAMRValidity.cxx`: `reader->LogDataOff()`.
- `src/Testing/Cxx/TestJet3D.cxx`: `reader->LogDataOn()`; the never-compiled
  `#ifdef LOG` variants of the contour/color array names and lookup-table
  range were dropped (keeping the non-LOG branch), and the stale
  `SelectColorArray("Log10(Density)")` (inside the `#ifdef ALL` block) was
  corrected to `"Density [N/cm^3]"`.

After the removal a full `make` in `buildDev` is clean (including both Cxx
tests) and `r.SMProxy.GetProperty("LogData")` returns `None` in pvbatch while
sibling properties (`Level`, ...) still resolve. The `VisIt/` copies of the
reader were left untouched, as elsewhere in this cleanup.

---

# Regression test: `src/Testing/Python/TestAMAZEFullFile.py.in`

The smoke test that validated this cleanup on the real 1.7 GB file
`DATA_M1.CygX-...ML=5-7...R=256.COOL.NT=0000023158...amr5` is now a permanent
regression test, following the existing `TestAMAZEReader.py.in` pattern
(`$<TARGET_FILE:pvAMAZEReader>` plugin path, registered through
`paraview_add_test_python` in `src/Testing/Python/CMakeLists.txt`).

It verifies, with values confirmed identical between the original and cleaned
code (A/B run with `git stash`):

- metadata: 6 levels, 384 blocks, timestep `4057652.8125`;
- array selection: selecting `Density` + `Velocity` loads exactly those two
  arrays per block, `Velocity` with 3 components;
- grid geometry: level-0 first block dims `(31, 31, 133)`, spacing `1.953125e9`;
- full data read: 19,805,712 cells over all 384 blocks, Density range
  `[4.655599e-15, 4.926240e-15]` (relative tolerance 1e-5);
- root `FieldData`: `CenterOfMass`, `Omega`, `Thermodynamics_EOS` present;
- stars port: 2 blocks, at least one with geometry, and a second
  `UpdatePipeline()` (star rebuild) leaves everything unchanged — this is the
  path that used to leak/hang before the cleanup.

Notes:

- The test skips (prints `SKIP`, exit code 0, `SKIP_REGULAR_EXPRESSION` set in
  CTest) when the data file is absent, so it is safe on machines without the
  `/local/data` volume. Override the path with the `AMAZE_FULL_TEST_FILE`
  environment variable.
- Manual run: `pvbatch TestAMAZEFullFile.py` after substituting the plugin
  path, e.g.
  `sed 's|\$<TARGET_FILE:pvAMAZEReader>|<build>/lib/paraview-6.2/plugins/pvAMAZEReader/pvAMAZEReader.so|' TestAMAZEFullFile.py.in > /tmp/t.py && pvbatch /tmp/t.py`
- Python CTests are currently not registered in `buildDev` because
  `PARAVIEW_USE_PYTHON` is not set in that build's cache; configure with
  `-DPARAVIEW_USE_PYTHON=ON` to enable both Python tests.
- Known pre-existing issue seen by the test: the `IsotrInfWind` (wind) star
  block is `None` in the stars output for this file — identical in the original
  code, so it is *not* asserted on; only the accretor star's geometry is
  required. Worth investigating separately (star count vs. built-stars
  mismatch in `BuildStars()`).
