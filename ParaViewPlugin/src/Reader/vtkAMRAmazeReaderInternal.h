/**
 * @class   vtkAMRAmazeReaderInternal
 * @brief   Consists of the low-level Amaze Reader used by the vtkAMRAmazeReader.
 *
 * @sa
 *  vtkAMRAmazeReader
 */

#ifndef vtkAMRAmazeReaderInternal_h
#define vtkAMRAmazeReaderInternal_h

#include "vtkAMRBox.h"
#include "vtkDoubleArray.h"
#include "vtk_hdf5.h"

#include <map>
#include <string>
#include <vector>

class vtkPolyData;
class vtkUniformGrid;

enum class ScaleType : int
{
  pc = 0,
  AU = 1,
  RSun = 2,
  NoScale = 3
};

enum class MapName : int
{
  NoMap = 0,
  Sphere_LogR = 1,
  DCR_Cart2Spheres = 2
};

constexpr int adG_MAXDIM = 3;        // maximum number of dimensions
constexpr int adG_LABELLENGTH = 64;  // maximum length of a label string
constexpr int adG_UNITLENGTH = 32;   // maximum length of a unit string

VTK_ABI_NAMESPACE_BEGIN

struct Thermodynamics_EOS
{
  double Adiabatic_Exponent;
  double Specific_Gasconstant;
  double Scaled_Specific_Gasconstant;
  double Scale_Factor;
  double Mean_Particle_Mass;
};

struct adG_grid_layout
{
  int    grid_nr;                   // number of grid, part of filename
  int    level;                     // level to which grid belongs
  int    dimensions[adG_MAXDIM];    // array with nx, ny (, nz)
  double origin[2 * adG_MAXDIM];    // min (followed by max) grid extent
  int    box_corners[2 * adG_MAXDIM]; // min (followed by max) grid extent in cell indices
};

struct adG_grid
{
  adG_grid_layout layout;
  vtkAMRBox       amrbox;
};

struct level
{
  int    GridsPerLevel;
  int    RefRatio;
  double DXs[3];
};

struct adG_component
{
  int   vec_len;               // vector-length of component
  char  label[adG_LABELLENGTH];// label of data component
  char  unit[adG_UNITLENGTH];  // unit for data component
  float scalefactor;           // scale-factor for component
};

struct AxiSymStarCurrent
{
  double Theta;
  double Radius;
  double Temperature;
  double Luminosity;
  double Epsilon;
  float  Omega[3];
  float  BField[3];
};

struct SpherSymStarCurrent
{
  double Radius;
  double Temperature;
  double Luminosity;
  double Omega[3];
  double BField[3];
};

struct newstar
{
  double StarTime;
  double CompRadiusFrac;
  double Mass;
  char   SpectralType[100];
  double Position[3];
  double Velocity[3];
  char   StarModel[100];
  char   StellarEvolution[100];
  char   StarModelFileName[100];
  char   InteractionModel[100];
  char   IActionEvolution[100];
  char   IActionModelFileName[100];
};

struct ParamToPhysicalMapping
{
  char   label[15]; // label of data component
  double Rmin;
  double Rmax;
  double Tmin;
  double Tmax;
  double Pmin;
  double Pmax;
  double StretchFactorOBSOLETE;
};

struct DCR_Mapping
{
  double Rmin;
  double Rmax;
  char   MapCase[6];
  char   MapLunarity[7];
  int    Dimension;
};

struct Fixed_Binary_Orbit
{
  double CenterofMass_x;
  double CenterofMass_y;
  double CenterofMass_z;
  double Eccentric_Anomaly;
  double Eccentricity_Ecc;
  double Inclination;
  double Initial_Orbital_Angle;
  double Initial_Phase;
  double Initial_Time;
  double Long_half_axis;
  double Long_half_axis1;
  double Long_half_axis2;
  double Omega;
  double Orbital_Period;
  double Orbital_Phase;
  double Orientation_of_Rotation;
  double Periastron_Angle;
  double Radius;
  double Reduced_Mass_Star1;
  double Reduced_Mass_Star2;
  double Star1angle;
  double Star1radius;
  double Star2radius;
  double sqrt_1_Ecc_2;
};

// ----------------------------------------------------------------------------
//                     Class  vtkAMRAmazeReaderInternal (begin)
// ----------------------------------------------------------------------------

class vtkAMRAmazeReaderInternal
{
public:
  vtkAMRAmazeReaderInternal();
  ~vtkAMRAmazeReaderInternal();

  vtkAMRAmazeReaderInternal(const vtkAMRAmazeReaderInternal&) = delete;
  void operator=(const vtkAMRAmazeReaderInternal&) = delete;

  void SetFileName(const char* fileName);

  int ReadMetaData(); // returns nb of stars

  int  GetBlockLevel(int domain) const;
  void GetSpacing(int level, double* spacing);
  vtkUniformGrid* GetAMRGrid(int blockIdx);
  vtkDoubleArray* ReadVisItVar(int domain, const char* varname);

  int BuildStars(); // returns nb of stars

  void MakeVariableNames();
  int  GridsPerLevels(int l) { return this->Levels[l].GridsPerLevel; }

  ScaleType ScaleChoice{ScaleType::NoScale};

  int    Dimensionality{0};
  int    NumberOfLevels{0};
  int    NumberOfComponents{0};
  int    NumberOfGrids{0};
  int    NumberOfSphericallySymmetricStars{0};
  int    NumberOfAxisSymmetricStars{0};
  MapName MappedGrids{MapName::NoMap};
  double AMAZETimeScalor{0.0};
  double AMAZETime{0.0}; // simulation time read from HDF5 file

  std::vector<adG_component> Labels;
  std::vector<adG_grid>      Grids;
  std::vector<level>         Levels;
  std::vector<vtkPolyData*>  Stars;
  std::vector<newstar>       newstars;
  ParamToPhysicalMapping     SphereLogRMappings[2]{};
  DCR_Mapping                DCR_Mappings{};
  Thermodynamics_EOS         thermodynamics_eos{};
  Fixed_Binary_Orbit         FixedBinaryOrbit{}; // initialized to all 0.0

  vtkTypeBool LengthScale{1}; // will automatically scale the grids to real length
  double      LengthScaleFactor{1.0};

protected:
  std::string FileName;
  int         DataScale{1};

  std::map<std::string, std::string> PVlabels;

  hid_t file_id{0};

private:
  void ReadHDF5GridsMetaData();
  void ReadHDF5VariablesMetaData();
  int  ReadHDF5MetaData(); // returns nb of stars

  void FindLevelAndBlock(int domain, int& level, int& block) const;
  int  FindDomainId(int level, int block);
  vtkDoubleArray* ReadVar(int levelId, int block, adG_component& variable);
  vtkPolyData* AxisSymStarSource(newstar* astar, AxiSymStarCurrent* axiStarData,
                                 int AngleResolution);
};

// ----------------------------------------------------------------------------
//                     Class  vtkAMRAmazeReaderInternal ( end )
// ----------------------------------------------------------------------------
VTK_ABI_NAMESPACE_END
#endif /* vtkAMRAmazeReaderInternal_h */
// VTK-HeaderTest-Exclude: vtkAMRAmazeReaderInternal.h
