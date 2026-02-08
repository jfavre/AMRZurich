// .NAME vtkCSCSAMRReader - Reads CSCSAMR files (in development)
// .SECTION Description
// This is an experimental CSCSAMR file reader. 

#ifndef __vtkAMAZEReader_h
#define __vtkAMAZEReader_h

class vtkPolyData;
class vtkUniformGrid;
class vtkRectilinearGrid;
class vtkStructuredGrid;
class vtkDoubleArray;
class vtkDataArray;
enum ScaleOption {pc=0, AU, RSun, NoScale};
enum MapName     {NoMap=0, Sphere_LogR, DCR_Cart2Spheres};
//enum Lunarity    {FULL=0, HALF=1, QUARTER=2}
#include <vector> // Needed for vector ivar
#include <map>
#include <string>

//#ifdef AMAZE_EXPORTS
//#include <vtkhdf5/hdf5.h>
//#else
#include <vtk_hdf5.h>
//#endif

#include <stddef.h>

#include "vtkObject.h"
#include "vtkAMRBox.h"

#define adG_MAXDIM		3	  // maximum number of dimensions
#define adG_NAMELENGTH	400	  // maximum length of a name string
#define adG_PATHLENGTH	1024  // maximum length of a path string
#define adG_LABELLENGTH	64	  // maximum length of a label string
#define adG_UNITLENGTH	32	  // maximum length of a unit string

#define adG_BINARY		0
#define adG_ASCII		1

typedef struct adG_grid
 {
  int       grid_nr;                   // number of grid, part of filename
  int       level;                     // level to which grid belongs
  int       dimensions[adG_MAXDIM];    // array with nx, ny (, nz)
  double    origin[2*adG_MAXDIM];      // min (followed by max) grid extent*/
  int       box_corners[2*adG_MAXDIM]; // min (followed by max) grid extent in cell indices with 
  vtkAMRBox amrbox;
} adG_grid;

typedef struct level
 {
 int          GridsPerLevel;
 int          RefRatio;
 double       DXs[3];
 } level;

typedef struct adG_component
 {	
  int           vec_len;                  // vector-length of component
  char label[adG_LABELLENGTH];   // label of data component
  char unit[adG_UNITLENGTH];     // unit for data component
  float         scalefactor;              // scale-factor for component
 } adG_component;

typedef struct AxiSymStarCurrent
{
  double Theta;
  double Radius;
  double Temperature;
  double Luminosity;
  double Epsilon;
  float  Omega[3];
  float  BField[3];
} AxiSymStarCurrent;

typedef struct SpherSymStarCurrent
{
  double Radius;
  double Temperature;
  double Luminosity;
  double  Omega[3];
  double  BField[3];
} SpherSymStarCurrent;

typedef struct newstar
{
  double StarTime;
  double CompRadiusFrac;
  double Mass;
  char   SpectralType[100];
  double  Position[3];
  double  Velocity[3];
  char   StarModel[100];
  char   StellarEvolution[100];
  char   StarModelFileName[100];
  char   InteractionModel[100];
  char   IActionEvolution[100];
  char   IActionModelFileName[100];
} newstar;

typedef struct ParamToPhysicalMapping
{
  char     label[15];   // label of data component
  double   Rmin;
  double   Rmax;
  double   Tmin;
  double   Tmax;
  double   Pmin;
  double   Pmax;
  double   StretchFactorOBSOLETE;
} ParamToPhysicalMapping;

typedef struct DCR_Mapping
{
  double   Rmin;
  double   Rmax;
  char     MapCase[6];
  char     MapLunarity[7];
  int      Dimension;
} DCR_Mapping;

class vtkAMAZEReader: public vtkObject
{
public:
  vtkTypeMacro(vtkAMAZEReader, vtkObject);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkAMAZEReader* New();

  //void SetFileName( const char * fileName );
  virtual void SetFileName(VTK_FILEPATH const char* filename);
  vtkGetFilePathMacro(FileName);
  
  vtkSetStringMacro(HDF5SaveFileName);
  vtkGetStringMacro(HDF5SaveFileName);

  // Description:
  // Get the number of components (similar to NumberOfComponents
  // in data arrays)
  vtkGetMacro(NumberOfComponents, int);

  vtkGetMacro(NumberOfLevels, int);
  vtkGetMacro(NumberOfGrids, int);

  vtkGetMacro(Dimensionality, int);

  vtkSetMacro(Time, double);
  vtkGetMacro(Time, double);

  vtkSetMacro(MaxLevelWrite, int);
  vtkGetMacro(MaxLevelWrite, int);

  // The range of valid levels values.
  int LevelRange[2];
  vtkGetVector2Macro(LevelRange, int);

  int LevelRead[2];
  vtkSetVector2Macro(LevelRead, int);
  vtkGetVector2Macro(LevelRead, int);

  vtkSetMacro(LogData, int);
  vtkGetMacro(LogData, int);
  vtkBooleanMacro(LogData, int);

  vtkSetMacro(LengthScale, int);
  vtkGetMacro(LengthScale, int);
  vtkBooleanMacro(LengthScale, int);

  vtkSetMacro(LengthScaleFactor, double);
  vtkGetMacro(LengthScaleFactor, double);

  vtkSetMacro(DataScale, int);
  vtkGetMacro(DataScale, int);
  vtkBooleanMacro(DataScale, int)

  vtkSetMacro(CellCentered, int);
  vtkGetMacro(CellCentered, int);
  vtkBooleanMacro(CellCentered, int);

  void  ReadMetaData();

  void FindLevelAndBlock(int domain, int &level, int &block) const;
  int GetBlockLevel(const int domain) const;
  int  FindDomainId(int level, int block);
  void GetSpacing(int level, double *spacing);
  vtkUniformGrid* ReadUniformGrid(int levelId, int domain);
  vtkRectilinearGrid* ReadRectilinearGrid(int domain);
  vtkStructuredGrid* ReadStructuredGrid(int domain);  // for Sphere_LogR
  vtkStructuredGrid* ReadStructuredGrid2(int domain); // for DCR_Cart2Spheres
  vtkDoubleArray* ReadVisItVar(int domain, const char *varname);
  vtkDoubleArray* ReadVar(int levelId, int block, adG_component&);
  void CheckVarSize(int levelId, int block, adG_component &variable);

  int BuildStars();
  vtkPolyData* GetStar(int domain);
  vtkPolyData* AxisSymStarSource(newstar *astar,
                                 struct AxiSymStarCurrent *axiStarData,
                                 int AngleResolution);
  
  ScaleOption            ScaleChoice;
  int MaxLevelRead;
  int MinLevelRead;
  int NumberOfLevels;
  int NumberOfComponents;
  int NumberOfGrids;
  int NumberOfSphericallySymmetricStars;
  int NumberOfAxisSymmetricStars;
  MapName MappedGrids;
  double TimeScalor;
  double Time; // simulation time read from HDF5 file

  std::vector<adG_component> Labels;
  std::vector<adG_grid>      Grids;
  std::vector<level>         Levels;
  std::vector<vtkPolyData*>  Stars;
  ParamToPhysicalMapping     SphereLogRMappings[2];
  DCR_Mapping                DCR_Mappings;
  std::map<std::string, int> VarNamesToLog;
  std::map<std::string, std::string> PVlabels;
  void ReadHDF5GridsMetaData(bool);
  void MakeVariableNames();
  int LengthScale; // will automatically scale the grids to real length
  double LengthScaleFactor;

protected:
  vtkAMAZEReader();
  ~vtkAMAZEReader();
  // The input file's name.
  hid_t file_id;
  //std::string FileName;
  char* FileName = nullptr;
  char* HDF5SaveFileName;
  int LogData; // will automatically calculate log10() for Density, Temperature and Pressure
  int Dimensionality;
  int DataScale;

  int CellCentered;

  int MaxLevelWrite;

  FILE *errs;
  void ReadHDF5MetaData();
  void ReadHDF5VariablesMetaData();

private:
};

#endif
