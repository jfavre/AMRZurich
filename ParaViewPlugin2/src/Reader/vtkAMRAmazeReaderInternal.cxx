
#include "vtkAMRAmazeReaderInternal.h"

#include "vtkAMRBox.h"
#include "vtkCellArray.h"
#include "vtkCharArray.h"

#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkErrorCode.h"
#include "vtkFloatArray.h"
#include "vtkMath.h"
#include "vtkNew.h"

#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkRectilinearGrid.h"
#include "vtkSphereSource.h"
#include "vtkStringArray.h"
#include "vtkStructuredGrid.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkUniformGrid.h"

#include <stdio.h>
#include <math.h>
#include <stddef.h>

//#define SINGLE_OUTPUT_PORT 1
//#define PARALLEL_DEBUG 1

#include <string>
#include <map>
#include <format>
#include <sstream>

#include <algorithm>
#include <vector>
#include <iostream>
//------------------------------------------------------------------------------
VTK_ABI_NAMESPACE_BEGIN

typedef struct interaction
{
  int StarNumber;
  int NTIME;
  int NAbund;
  double CompRadius;
  double MassLoss;
  double VInf;
  double Temp;
} interaction;


//----------------------------------------------------------------------------
static hid_t Create_Interaction_Compound()
{
  hid_t    id;

  id = H5Tcreate(H5T_COMPOUND, sizeof(interaction));

  H5Tinsert(id, "Star Number", HOFFSET(interaction, StarNumber),
            H5T_NATIVE_INT);

  H5Tinsert(id, "NTIME", HOFFSET(interaction, NTIME),
		    H5T_NATIVE_INT);

  H5Tinsert(id, "NAbund", HOFFSET(interaction, NAbund),
            H5T_NATIVE_INT);

  H5Tinsert(id, "CompRadius", 
            HOFFSET(interaction, CompRadius), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "MassLoss", 
            HOFFSET(interaction, MassLoss), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "VInf", 
            HOFFSET(interaction, VInf), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Temp", 
            HOFFSET(interaction, Temp), H5T_NATIVE_DOUBLE);

  return id;
};

typedef struct model
{
  char Model[20];
  char ModelFileName[20];
  char IntActType[20];
  char IntActModel[20];
  char IModelFileName[20];
  int NTime;
  int NTimePos;
} model;

static hid_t Create_StarModel_Compound()
{
  hid_t    id, labelstring;

  labelstring = H5Tcopy(H5T_C_S1);
                H5Tset_size(labelstring, 20);
                H5Tset_strpad(labelstring, H5T_STR_NULLTERM);

  id = H5Tcreate(H5T_COMPOUND, sizeof(model));

  H5Tinsert(id, "Model", HOFFSET(model, Model), labelstring);
  H5Tinsert(id, "ModelFileName", HOFFSET(model, ModelFileName), labelstring);
  H5Tinsert(id, "IntActType", HOFFSET(model, IntActType), labelstring);
  H5Tinsert(id, "IntActModel", HOFFSET(model, IntActModel), labelstring);
  H5Tinsert(id, "I-ModelFileName", HOFFSET(model, IModelFileName), labelstring);
  H5Tinsert(id, "NTime", HOFFSET(model, NTime), H5T_NATIVE_INT);
  H5Tinsert(id, "NTimePos", HOFFSET(model, NTimePos), H5T_NATIVE_INT);

  H5Tclose(labelstring);

  return id;
};

static hid_t Create_NewStar_Compound()
{
  hid_t    id, labelstring, array1_id;
  hsize_t  dim[1];
  labelstring = H5Tcopy(H5T_C_S1);
                H5Tset_size(labelstring, 100);
                H5Tset_strpad(labelstring, H5T_STR_NULLTERM);

  dim[0] = 3;
  array1_id = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dim);

  id = H5Tcreate(H5T_COMPOUND, sizeof(newstar));
  H5Tinsert(id, "StarTime [y]",          HOFFSET(newstar, StarTime), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "CompRadiusFrac",        HOFFSET(newstar, CompRadiusFrac), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Mass [msol]",           HOFFSET(newstar, Mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Spectral Type",         HOFFSET(newstar, SpectralType), labelstring);
  H5Tinsert(id, "Position",              HOFFSET(newstar, Position), array1_id);
  H5Tinsert(id, "Velocity",              HOFFSET(newstar, Velocity), array1_id);
  H5Tinsert(id, "Star Model",            HOFFSET(newstar, StarModel), labelstring);
  H5Tinsert(id, "Stella rEvolution",     HOFFSET(newstar, StellarEvolution), labelstring);
  H5Tinsert(id, "Star Model FileName",   HOFFSET(newstar, StarModelFileName), labelstring);
  H5Tinsert(id, "Interaction Model",     HOFFSET(newstar, InteractionModel), labelstring);
  H5Tinsert(id, "I-ActionEvolution",     HOFFSET(newstar, IActionEvolution), labelstring);
  H5Tinsert(id, "I-ActionModelFileName", HOFFSET(newstar, IActionModelFileName), labelstring);
  H5Tclose(array1_id);
  H5Tclose(labelstring);

  return id;
};

typedef struct star
{
  char   Type[20];
  double Position[3];
  double Velocity[3];
  double Radius;
  double Mass;
  double Temperature;
  double Luminosity;
  double Rotation[3];
  double BField[3];
  char   Interaction[20];
} star;

static hid_t Create_Star_Compound()
{
  hid_t    id, labelstring, array1_id;
  hsize_t  dim[1];
  labelstring = H5Tcopy(H5T_C_S1);
                H5Tset_size(labelstring, 20);
                H5Tset_strpad(labelstring, H5T_STR_NULLTERM);

  dim[0] = 3;
  array1_id = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dim);

  id = H5Tcreate(H5T_COMPOUND, sizeof(star));

  H5Tinsert(id, "Type",        HOFFSET(star, Type), labelstring);
  H5Tinsert(id, "Position",    HOFFSET(star, Position), array1_id);
  H5Tinsert(id, "Velocity",    HOFFSET(star, Velocity), array1_id);
  H5Tinsert(id, "Radius",      HOFFSET(star, Radius), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Mass",        HOFFSET(star, Mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Temperature", HOFFSET(star, Temperature), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Luminosity",  HOFFSET(star, Luminosity), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Rotation",    HOFFSET(star, Rotation), array1_id);
  H5Tinsert(id, "B-Field",     HOFFSET(star, BField), array1_id);
  H5Tinsert(id, "Interaction", HOFFSET(star, Interaction), labelstring);
  H5Tclose(array1_id);
  H5Tclose(labelstring);

  return id;
};

static hid_t Create_AxiSymStar_Compound()
{
  hsize_t  dim[1];
  dim[0] = 3;
  hsize_t array1_id = H5Tarray_create2(H5T_NATIVE_FLOAT, 1, dim);

  hid_t   id = H5Tcreate(H5T_COMPOUND, sizeof(AxiSymStarCurrent));

  H5Tinsert(id, "Theta",       HOFFSET(AxiSymStarCurrent, Theta),       H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Radius",      HOFFSET(AxiSymStarCurrent, Radius),      H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Temperature", HOFFSET(AxiSymStarCurrent, Temperature), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Luminosity",  HOFFSET(AxiSymStarCurrent, Luminosity),  H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Epsilon",     HOFFSET(AxiSymStarCurrent, Epsilon),     H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Omega",       HOFFSET(AxiSymStarCurrent, Omega),       array1_id);
  H5Tinsert(id, "BField",      HOFFSET(AxiSymStarCurrent, BField),      array1_id);
  H5Tclose(array1_id);

  return id;
}

static hid_t Create_SpherSymStar_Compound()
{
  hsize_t  dim[1];
  dim[0] = 3;
  hsize_t array1_id = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dim);

  hid_t   id = H5Tcreate(H5T_COMPOUND, sizeof(SpherSymStarCurrent));

  H5Tinsert(id, "Radius",      HOFFSET(SpherSymStarCurrent, Radius),      H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Temperature", HOFFSET(SpherSymStarCurrent, Temperature), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Luminosity",  HOFFSET(SpherSymStarCurrent, Luminosity),  H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Omega",       HOFFSET(SpherSymStarCurrent, Omega),       array1_id);
  H5Tinsert(id, "BField",      HOFFSET(SpherSymStarCurrent, BField),      array1_id);
  H5Tclose(array1_id);

  return id;
}

vtkStandardNewMacro(vtkAMRAmazeReaderInternal);

vtkAMRAmazeReaderInternal::vtkAMRAmazeReaderInternal()
{
  //this->FileName = NULL;
  this->file_id = 0;
  this->Dimensionality = 0;
  this->NumberOfLevels = 0;
  this->NumberOfComponents = 0;
  this->Labels.clear();
  this->Grids.clear();
  this->Stars.clear();
  this->LogDataOn();
  this->DataScaleOn();
  this->CellCenteredOff();
  this->DebugOff();
  this->MaxLevelWrite = -1;

  this->LevelRead[0] = -1;
  this->LevelRead[1] = -1;

  this->LevelRange[0] = -1;
  this->LevelRange[1] = -1;
  this->LengthScale = true; // GUI will ALWAYS overwrite that value
  this->LengthScaleFactor = 1;
  this->ScaleChoice = NoScale;
  //cerr << "AMAZEConstructor\n";
  this->VarNamesToLog["Density"] = 1;
  this->VarNamesToLog["Pressure"] = 1;
  this->VarNamesToLog["Temperature"] = 1;

  this->NumberOfSphericallySymmetricStars = 0;
  this->NumberOfAxisSymmetricStars = 0;
  this->MappedGrids = NoMap;
}

vtkAMRAmazeReaderInternal::~vtkAMRAmazeReaderInternal()
{
  int i;
  //cerr << "AMAZEDestructor\n";
  this->SetFileName(nullptr);
    
  for(i=0; i < this->Stars.size(); i++)
    {
    (this->Stars[i])->Delete();
    }
  this->Stars.clear();
  
  this->Labels.clear();
  
  this->Grids.clear();
  
  this->VarNamesToLog.clear();
  if(this->file_id)
    {
    H5Fclose(this->file_id);
    this->file_id = 0;
    //cerr << "465: H5Fclose( " << this->FileName << ")\n";
    }
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReaderInternal::ReadMetaData()
{
  int levelId, GridId, i, node_veclen;
  double time, time_scalor;

  this->Levels.clear();
  this->Labels.clear();
  if (!this->FileName || *this->FileName == 0)
    {
    //this->SetErrorCode(vtkErrorCode::NoFileNameError);
    //vtkErrorMacro(<< "Must specify adG file");
    return 0;
    }
  //cout << __LINE__ << ": H5Fopen( " << this->FileName << ")\n";
  this->file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if(this->file_id<0)
    {
    std::cerr << "file could not be opened. Check filename " << endl;
    return 0;
    }
  auto nbstars = this->ReadHDF5MetaData();

  this->Levels.resize(this->NumberOfLevels);

  this->Labels.resize(this->NumberOfComponents);

  this->ReadHDF5VariablesMetaData();

  this->Grids.resize(this->NumberOfGrids);

  this->ReadHDF5GridsMetaData(false);

  this->CheckVarSize(0, 0, this->Labels[0]);

  if(this->file_id)
    {
    H5Fclose(this->file_id);
    this->file_id = 0;
    //cout << __LINE__ << ": H5Fclose( " << this->FileName << ")\n";
    }

///////////////////////////////////////////////////////////////////
  int current_level = -1;
  std::vector<adG_grid> &grid = this->Grids;

  for (i = 0; i < this->NumberOfGrids; i++)
    {
    this->Grids[i].amrbox = vtkAMRBox(grid[i].layout.box_corners, grid[i].layout.box_corners+3);
    if(this->Grids[i].layout.level != current_level)
      { // first grid of a given level
      current_level = this->Grids[i].layout.level;
      this->Levels[current_level].GridsPerLevel = 1;
      }
    else
      {
      this->Levels[current_level].GridsPerLevel++;
      }
    }

  this->LevelRange[0] = 0;
  this->LevelRange[1] = this->NumberOfLevels-1;

  this->MinLevelRead = this->LevelRange[0];
  this->MaxLevelRead = this->LevelRange[1];

  int firstLevel = this->MinLevelRead;
  int lastLevel = this->MaxLevelRead;

  /*
  cout << __LINE__ << "\tDimensionality: " << this->Dimensionality 
                 << "\n\tNumberOfComponents: " << this->NumberOfComponents
                 << "\n\tNumberOfLevels: " << this->NumberOfLevels
                 << "\n\tNumberOfGrids: " << this->NumberOfGrids << endl;
  */
  //grid = this->Grids;
  //int size = grid[0].dimensions[0];

  return nbstars;
}

// deciding on adding the Log10 prefix to the name can only be done after UpdateInformation
// so we make this a separate call
void vtkAMRAmazeReaderInternal::MakeVariableNames()
{
  for(int c=0; c < this->NumberOfComponents; c++)
    {
    std::ostringstream varName;
    if(strlen(Labels[c].unit))
      {

      if(this->VarNamesToLog.find(Labels[c].label) == this->VarNamesToLog.end())
        {
        varName << Labels[c].label << " [" << Labels[c].unit << "]" << ends;
        }
      else
        {
        if(this->LogData)
          varName << "Log10(" << Labels[c].label << ")" << " [" << Labels[c].unit << "]" << ends;
        else
          varName << Labels[c].label << " [" << Labels[c].unit << "]" << ends;
        }
      }
    else
      { // unit-less variables such as Mach Number are never log() anyway, so no need to check for LogData
      varName << Labels[c].label << ends;
      }
    PVlabels[(const char *)Labels[c].label] = varName.str();
    //delete[] varName.str();
    }
}

int vtkAMRAmazeReaderInternal::ReadHDF5MetaData()
{
  hid_t   root_id, apr_root_id, StarsDS, models_root_id;
  hid_t   attr1;
  herr_t  status;
// turn off error reporting
  H5E_auto2_t func;
  void *client_data;
  H5Eget_auto2(H5E_DEFAULT, &func, &client_data);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

  root_id = H5Gopen(this->file_id, "/", H5P_DEFAULT);

  attr1 = H5Aopen_name(root_id, "NumberLevels");
  status = H5Aread(attr1, H5T_NATIVE_INT, &this->NumberOfLevels);
  status = H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "NumberOfComponents");
  status = H5Aread(attr1, H5T_NATIVE_INT, &this->NumberOfComponents);
  status = H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "NumberOfGrids");
  status = H5Aread(attr1, H5T_NATIVE_INT, &this->NumberOfGrids);
  status = H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "Dimensionality");
  status = H5Aread(attr1, H5T_NATIVE_INT, &this->Dimensionality);
  status = H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "Time");
  status = H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->AMAZETime);
  if (status < 0)
    {
    vtkGenericWarningMacro("Failed to open Time Attribute " << endl);
    return 0;
    }
  status = H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "Time Scaling Factor");
  status = H5Aread(attr1, H5T_NATIVE_DOUBLE, &AMAZETimeScalor);
  if (status < 0)
    {
    vtkGenericWarningMacro("Failed to open Time Scaling Factor Attribute " << endl);
    }
  else
    this->AMAZETime /= this->AMAZETimeScalor;
  status = H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "Length Scale Factor");
  if(attr1 >0 )
    {
    status = H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->LengthScaleFactor);
    status = H5Aclose(attr1);
    }
  else
    this->LengthScaleFactor = 1.0;

  switch(this->ScaleChoice)
    {
    case 0: // pc
      this->LengthScaleFactor  *= 3.08567782e18;
      std::cout << "!!!\nUsing PARSEC with Length Scale Factor * by 3.08567782e18 = " << this->LengthScaleFactor << "!!!\n";
    break;
    case 1: // AU
      this->LengthScaleFactor  *= 1.49597870700e13;
      std::cout << "!!!\nUsing AU with Length Scale Factor * by 1.49597870700e13 = " << this->LengthScaleFactor << "!!!\n";
    break;
    case 2: // RSun
      this->LengthScaleFactor  *= 6.96342e10;
      std::cout << "!!!\nUsing RSun with Length Scale Factor * by 6.96342e10 = " << this->LengthScaleFactor << "!!!\n";
    break;
    case 3:
      //std::cout << "!!!\nLength Scale Factor  is untouched!!!\n";
    break;
    }

  if(H5Lexists(root_id, "/Map", H5P_DEFAULT))
    {
    hid_t map_id = H5Gopen(root_id, "/Map", H5P_DEFAULT);
    attr1 = H5Aopen_name(map_id, "Map Type");
    hid_t t1 = H5Aget_type(attr1);
    char map_type[30];
    status = H5Aread(attr1, t1, map_type);
    H5Tclose(t1);

    if(!strncmp(map_type, "Sphere-LogR", 11))
      {
      this->MappedGrids = Sphere_LogR;
      std::cerr << "Using Mapped Grids:" << this->MappedGrids << endl;
      }
    else if(!strncmp(map_type, "DCR_Cart2Spheres", 16))
      {
      this->MappedGrids = DCR_Cart2Spheres;
      std::cerr << "Using Mapped Grids:" << this->MappedGrids << endl;
      }
    status = H5Aclose(attr1);
    status = H5Gclose(map_id);
    }
  
  H5Gclose(root_id);
  
  int nb_stars=0;

  apr_root_id = H5Gopen(this->file_id, "/APR_StellarSystems", H5P_DEFAULT);
  models_root_id = H5Gopen(apr_root_id, "Stars", H5P_DEFAULT);
  hid_t dataspace;
  hsize_t dims_out[1]; // we assume rank = 1

  if((StarsDS = H5Dopen(models_root_id, "Stars", H5P_DEFAULT)) >= 0)
    {
    dataspace = H5Dget_space(StarsDS);
    status = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    nb_stars = dims_out[0];
    H5Sclose(dataspace);
    H5Gclose(StarsDS);
    }

  H5Gclose(models_root_id);
  H5Gclose(apr_root_id);
  H5Eset_auto2(H5E_DEFAULT, func, client_data);
  return nb_stars;
} // end of ReadHDF5MetaData()



VTK_ABI_NAMESPACE_END
