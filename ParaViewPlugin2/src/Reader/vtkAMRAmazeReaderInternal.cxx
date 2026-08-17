
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
#include "vtkPointSet.h"
#include "vtkRectilinearGrid.h"
#include "vtkSphereSource.h"
#include "vtkStringArray.h"
#include "vtkStructuredGrid.h"
#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include "vtkUniformGrid.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

VTK_ABI_NAMESPACE_BEGIN

struct interaction
{
  int StarNumber;
  int NTIME;
  int NAbund;
  double CompRadius;
  double MassLoss;
  double VInf;
  double Temp;
};

struct model
{
  char Model[20];
  char ModelFileName[20];
  char IntActType[20];
  char IntActModel[20];
  char IModelFileName[20];
  int NTime;
  int NTimePos;
};

struct star
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
};


//----------------------------------------------------------------------------
static hid_t Create_Interaction_Compound()
{
  hid_t    id;

  id = H5Tcreate(H5T_COMPOUND, sizeof(interaction));

  H5Tinsert(id, "Star Number", HOFFSET(interaction, StarNumber), H5T_NATIVE_INT);
  H5Tinsert(id, "NTIME", HOFFSET(interaction, NTIME), H5T_NATIVE_INT);
  H5Tinsert(id, "NAbund", HOFFSET(interaction, NAbund), H5T_NATIVE_INT);
  H5Tinsert(id, "CompRadius", HOFFSET(interaction, CompRadius), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "MassLoss", HOFFSET(interaction, MassLoss), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "VInf", HOFFSET(interaction, VInf), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Temp", HOFFSET(interaction, Temp), H5T_NATIVE_DOUBLE);

  return id;
};



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

vtkAMRAmazeReaderInternal::vtkAMRAmazeReaderInternal()
{
  this->Dimensionality = 0;
  this->NumberOfLevels = 0;
  this->NumberOfComponents = 0;
  this->Labels.clear();
  this->Grids.clear();
  this->Stars.clear();
  this->DataScale = 1;

  this->ScaleChoice = ScaleType::NoScale;
  //cerr << "AMAZEConstructor\n";
  this->VarNamesToLog["Density"] = 1;
  this->VarNamesToLog["Pressure"] = 1;
  this->VarNamesToLog["Temperature"] = 1;

  this->NumberOfSphericallySymmetricStars = 0;
  this->NumberOfAxisSymmetricStars = 0;
  this->MappedGrids = MapName::NoMap;
}

vtkAMRAmazeReaderInternal::~vtkAMRAmazeReaderInternal()
{
  int i;
  //cerr << "AMAZEDestructor\n";
  this->SetFileName(nullptr);
    
  for (auto* starPtr : this->Stars)
  {
    if (starPtr)
    {
      starPtr->Delete();
    }
  }
  this->Stars.clear();
  
  this->Labels.clear();
  
  this->Grids.clear();
  
  this->VarNamesToLog.clear();

  if (this->file_id > 0)
  {
    H5Fclose(this->file_id);
    this->file_id = -1;
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

  if(this->file_id > 0)
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

void vtkAMRAmazeReaderInternal::CheckVarSize(int levelId, int block, adG_component &variable)
{
  hid_t level_root_id, grid_root_id, dataset_id, mem_space_id;
  int domain = this->FindDomainId(levelId, block);
  /*cerr << "domain = " << domain
       << ", level = " << levelId
       << ", block = " << block
       << ", varname = " << variable.label
       << endl;*/
  const auto& grid = this->Grids[domain];

  level_root_id = H5Gopen(this->file_id, std::format("/Level {}", levelId).c_str(), H5P_DEFAULT);
  if(level_root_id < 0)
    std::cerr << "bad level_root_id returned\n";
  else
    {
    std::string lname = std::format("Grid {}", grid.layout.grid_nr);
    if(H5Lexists(level_root_id, lname.c_str(), H5P_DEFAULT))
      {
      grid_root_id = H5Gopen(level_root_id, lname.c_str(), H5P_DEFAULT);
      if(grid_root_id < 0)
        std::cerr << "ReadVar(): bad grid_root_id returned\n";

      int nvals = grid.layout.dimensions[0] * grid.layout.dimensions[1] * grid.layout.dimensions[2];
  /*
  cerr << lname << ":" << PVlabels[(const char *)variable.label] << "("<< variable.vec_len << "," << nvals << ")\n";
  cerr << "nvals = " << grid.layout.dimensions[0] << "x"<<
                    grid.layout.dimensions[1]<< "x"<<
                    grid.layout.dimensions[2]<< endl;
*/
      dataset_id = H5Dopen(grid_root_id, (const char *) variable.label, H5P_DEFAULT);
      if(dataset_id < 0)
        {
        std::cerr << "error opening HDF5 dataset for var " << variable.label << endl;
        }
      hid_t space_id;
      hsize_t dims[3], maxdims;
      space_id = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(space_id, dims, NULL );
    
      H5Sclose(space_id);
      H5Dclose(dataset_id);
      H5Gclose(grid_root_id);
      H5Gclose(level_root_id);

      if (nvals == static_cast<int>(dims[0]))
        {
        //this->CellCenteredOff();
        }
      else
        {
        //this->CellCenteredOn();
        }
      return;
      }
    }
}

void vtkAMRAmazeReaderInternal::ReadHDF5GridsMetaData(bool shiftedGrid)
{
  hid_t   root_id, dataset, adG_grid_id, mapping_id, label1, label2, unitstring;
  herr_t  status;
  hid_t   attr1, array0_id, array1_id, array2_id;
  hsize_t  dim[2];

  if(shiftedGrid)
    {
    if(!this->file_id)
      {
      this->file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
      }
    //cout << "reading the shifted grid info\n";
    root_id = H5Gopen(this->file_id, "/", H5P_DEFAULT);
    dataset = H5Dopen(root_id, "Shifted Grid Info", H5P_DEFAULT);
    if(dataset < 0)
      {
      std::cerr << "failed to find shifted grid info. Returning without action\n";
      H5Gclose(root_id);
      return;
      }
    if(this->file_id)
      {
      H5Fclose(this->file_id);
      this->file_id = 0;
      //cerr << "584: H5Fclose( " << this->FileName << ")\n";
      }
    }
  else
    {
    root_id = H5Gopen(this->file_id, "/", H5P_DEFAULT);
    dataset = H5Dopen(root_id, "Grid Info", H5P_DEFAULT);
    }
  adG_grid_id = H5Tcreate (H5T_COMPOUND, sizeof(adG_grid));

  H5Tinsert(adG_grid_id, "grid number", offsetof(adG_grid_layout, grid_nr), H5T_NATIVE_INT);
  H5Tinsert(adG_grid_id, "level", HOFFSET(adG_grid_layout, level), H5T_NATIVE_INT);

  dim[0] = adG_MAXDIM;
  array0_id = H5Tarray_create2(H5T_NATIVE_INT, 1, dim);

  dim[0] = 2 * adG_MAXDIM;
  array1_id = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, dim);

  dim[0] = 2 * adG_MAXDIM;
  array2_id = H5Tarray_create2(H5T_NATIVE_INT, 1, dim);

  H5Tinsert(adG_grid_id, "dimensions", HOFFSET(adG_grid_layout, dimensions),
            array0_id);
  H5Tinsert(adG_grid_id, "origin", HOFFSET(adG_grid_layout, origin), array1_id);
  H5Tinsert(adG_grid_id, "box corners", HOFFSET(adG_grid_layout, box_corners), array2_id);

  H5Tclose(array0_id);
  H5Tclose(array1_id);
  H5Tclose(array2_id);

  status = H5Dread(dataset, adG_grid_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->Grids[0]);
  H5Tclose(adG_grid_id);
  H5Dclose(dataset);

  for (int levelId=0; levelId < this->NumberOfLevels; levelId++)
    {
    std::string name = std::format("Level {}", levelId);

    attr1 = H5Aopen_by_name(root_id, name.c_str(), "Refinement Ratio", H5P_DEFAULT, H5P_DEFAULT);
    status = H5Aread(attr1, H5T_NATIVE_INT, &this->Levels[levelId].RefRatio);
    status = H5Aclose(attr1);

    attr1 = H5Aopen_by_name(root_id, name.c_str(), "Cell Edge Length", H5P_DEFAULT, H5P_DEFAULT);
    status = H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->Levels[levelId].DXs[0]);
    status = H5Aclose(attr1);

    if(H5Aexists_by_name(root_id, name.c_str(), "Cell Edge Length Y", H5P_DEFAULT) > 0)
      {
      attr1 = H5Aopen_by_name(root_id, name.c_str(), "Cell Edge Length Y", H5P_DEFAULT, H5P_DEFAULT);
      if(attr1 >= 0)
        {
        status = H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->Levels[levelId].DXs[1]);
        status = H5Aclose(attr1);
        }
      else
        this->Levels[levelId].DXs[1] = this->Levels[levelId].DXs[0];
      }
    else
      this->Levels[levelId].DXs[1] = this->Levels[levelId].DXs[0];

    if(H5Aexists_by_name(root_id, name.c_str(), "Cell Edge Length Z", H5P_DEFAULT) > 0)
      {
      attr1 = H5Aopen_by_name(root_id, name.c_str(), "Cell Edge Length Z", H5P_DEFAULT, H5P_DEFAULT);
      if(attr1 >= 0)
        {
        status = H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->Levels[levelId].DXs[2]);
        status = H5Aclose(attr1);
        }
      else
        this->Levels[levelId].DXs[2] = this->Levels[levelId].DXs[0];
      }
    else
      this->Levels[levelId].DXs[2] = this->Levels[levelId].DXs[0];

    //H5Gclose(level_root_id);
    //cerr << this->Levels[levelId].RefRatio << " ," << this->Levels[levelId].DXs << endl;
    }
  H5Gclose(root_id);

  if((this->MappedGrids !=  MapName::NoMap) && (root_id = H5Gopen(this->file_id, "/Map", H5P_DEFAULT) )>= 0)
    {
    dataset = H5Dopen(root_id, "Map Parameter", H5P_DEFAULT);
    if(dataset < 0)
      {
      std::cerr << "error opening Map_Parameter\n";
      }
    switch(this->MappedGrids) {
      case MapName::Sphere_LogR:
      label1 = H5Tcopy(H5T_C_S1);
      H5Tset_size(label1, 15);
      H5Tset_strpad(label1, H5T_STR_NULLTERM);

      mapping_id = H5Tcreate(H5T_COMPOUND, sizeof(ParamToPhysicalMapping));

      H5Tinsert(mapping_id, "Space Name", HOFFSET(ParamToPhysicalMapping, label), label1);
      H5Tinsert(mapping_id, "Min Radius", HOFFSET(ParamToPhysicalMapping, Rmin), H5T_NATIVE_DOUBLE);
      H5Tinsert(mapping_id, "Max Radius", HOFFSET(ParamToPhysicalMapping, Rmax), H5T_NATIVE_DOUBLE);
      H5Tinsert(mapping_id, "Min Theta",  HOFFSET(ParamToPhysicalMapping, Tmin), H5T_NATIVE_DOUBLE);
      H5Tinsert(mapping_id, "Max Theta",  HOFFSET(ParamToPhysicalMapping, Tmax), H5T_NATIVE_DOUBLE);
      H5Tinsert(mapping_id, "Min Phi",    HOFFSET(ParamToPhysicalMapping, Pmin), H5T_NATIVE_DOUBLE);
      H5Tinsert(mapping_id, "Max Phi",    HOFFSET(ParamToPhysicalMapping, Pmax), H5T_NATIVE_DOUBLE);
      H5Tinsert(mapping_id, "Stretch Factor", HOFFSET(ParamToPhysicalMapping, StretchFactorOBSOLETE), H5T_NATIVE_DOUBLE);
      H5Tclose(label1);

      if ((status = H5Dread(dataset, mapping_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->SphereLogRMappings)) < 0)
        {
        std::cerr << "error reading Map_Parameter\n";
        }
/*
      cerr << "\nparam =\t" << param[0].label << "\n\t";
      cerr << param[0].Rmin<< "\n\t";
      cerr << param[0].Rmax << "\n\t";
      cerr << param[0].Tmin << "\n\t";
      cerr << param[0].Tmax << "\n\t";
      cerr << param[0].Pmin << "\n\t";
      cerr << param[0].Pmax << "\n\t";
      cerr << param[0].StretchFactorOBSOLETE << "\n\n";
*/
      break;
      
      case MapName::DCR_Cart2Spheres:
      label1 = H5Tcopy(H5T_C_S1);
      H5Tset_size(label1, 6);
      H5Tset_strpad(label1, H5T_STR_NULLTERM);
      label2 = H5Tcopy(H5T_C_S1);
      H5Tset_size(label2, 7);
      H5Tset_strpad(label2, H5T_STR_NULLTERM);

      mapping_id = H5Tcreate(H5T_COMPOUND, sizeof(DCR_Mapping));


      H5Tinsert(mapping_id, "Min Radius Physical Space", HOFFSET(DCR_Mapping, Rmin),        H5T_NATIVE_DOUBLE);
      H5Tinsert(mapping_id, "Max Radius Physical Space", HOFFSET(DCR_Mapping, Rmax),        H5T_NATIVE_DOUBLE);
      H5Tinsert(mapping_id, "Map Case",                  HOFFSET(DCR_Mapping, MapCase),     label1);
      H5Tinsert(mapping_id, "Map Lunarity",              HOFFSET(DCR_Mapping, MapLunarity), label2);
      H5Tinsert(mapping_id, "Dimensionality",            HOFFSET(DCR_Mapping, Dimension),   H5T_NATIVE_INT);
      H5Tclose(label1);
      H5Tclose(label2);

      if ((status = H5Dread(dataset, mapping_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->DCR_Mappings)) < 0)
        {
        std::cerr << "error reading Map_Parameter\n";
        }
      std::cerr << "\nparam =\t"  << "\n\t";
      std::cerr << DCR_Mappings.Rmin<< "\n\t";
      std::cerr << DCR_Mappings.Rmax << "\n\t";
      std::cerr << DCR_Mappings.MapCase << "\n\t";
      std::cerr << DCR_Mappings.MapLunarity << "\n\t";
      std::cerr << DCR_Mappings.Dimension << "\n\n";
      break;
      
      default:
        std::cerr << "error finding an implemented Mapping code\n";
      break;
    }
    H5Tclose(mapping_id); 
    H5Dclose(dataset);
    H5Gclose(root_id);
    }
} // ReadHDF5GridsMetaData()

void vtkAMRAmazeReaderInternal::ReadHDF5VariablesMetaData()
{
  hid_t   root_id, dataset, adG_component_id, labelstring, unitstring;
  herr_t  status;
  
  root_id = H5Gopen(this->file_id, "/", H5P_DEFAULT);
  dataset = H5Dopen(root_id, "Variable Info", H5P_DEFAULT);

  labelstring = H5Tcopy(H5T_C_S1);
                H5Tset_size(labelstring, adG_LABELLENGTH);
                H5Tset_strpad(labelstring, H5T_STR_NULLTERM);
  unitstring =  H5Tcopy(H5T_C_S1);
                H5Tset_size(unitstring, adG_UNITLENGTH-1); 
                H5Tset_strpad(unitstring, H5T_STR_NULLTERM);

  adG_component_id = H5Tcreate(H5T_COMPOUND, sizeof(adG_component));
  
  H5Tinsert(adG_component_id, "vector length", HOFFSET(adG_component, vec_len),
            H5T_NATIVE_INT);
  
  H5Tinsert(adG_component_id, "Variable Name", HOFFSET(adG_component, label),
                    labelstring);
    
  H5Tinsert(adG_component_id, "Variable Unit", HOFFSET(adG_component, unit),
            unitstring);

  H5Tinsert(adG_component_id, "scale factor",
            HOFFSET(adG_component, scalefactor), H5T_NATIVE_FLOAT);

  H5Tclose(labelstring);
  H5Tclose(unitstring);

  status = H5Dread(dataset, adG_component_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Labels[0]);

  H5Tclose(adG_component_id);
  H5Dclose(dataset);
  H5Gclose(root_id);
}


// deciding on adding the Log10 prefix to the name can only be done after UpdateInformation
// so we make this a separate call
void vtkAMRAmazeReaderInternal::MakeVariableNames()
{
  for(int c=0; c < this->NumberOfComponents; c++)
    {
    std::ostringstream varName;
    if(std::strlen(Labels[c].unit) > 0)
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
    case ScaleType::pc:
      this->LengthScaleFactor  *= 3.08567782e18;
      std::cout << "!!!\nUsing PARSEC with Length Scale Factor * by 3.08567782e18 = " << this->LengthScaleFactor << "!!!\n";
    break;
    case ScaleType::AU:
      this->LengthScaleFactor  *= 1.49597870700e13;
      std::cout << "!!!\nUsing AU with Length Scale Factor * by 1.49597870700e13 = " << this->LengthScaleFactor << "!!!\n";
    break;
    case ScaleType::RSun:
      this->LengthScaleFactor  *= 6.96342e10;
      std::cout << "!!!\nUsing RSun with Length Scale Factor * by 6.96342e10 = " << this->LengthScaleFactor << "!!!\n";
    break;
    case ScaleType::NoScale:
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
      this->MappedGrids = MapName::Sphere_LogR;
      //std::cerr << "Using Mapped Grids:" << this->MappedGrids << endl;
      }
    else if(!strncmp(map_type, "DCR_Cart2Spheres", 16))
      {
      this->MappedGrids = MapName::DCR_Cart2Spheres;
      //std::cerr << "Using Mapped Grids:" << this->MappedGrids << endl;
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

// returns the global id for the given block
int vtkAMRAmazeReaderInternal::FindDomainId(int level, int block)
{     
  int domain = 0;
  for(int l=0; l < level; l++)
    domain += this->Levels[l].GridsPerLevel;
  domain += block;
  if (domain >= this->NumberOfGrids)
   
    std::cerr << "level or block too high for this dataset\n";
  return domain;
}

void vtkAMRAmazeReaderInternal::GetSpacing(int level, double *spacing)
{
  if(this->Dimensionality == 3)  // a 3D grid
  {
    if(this->LengthScale)
    {
      spacing[2] = this->Levels[level].DXs[2] / this->LengthScaleFactor;
      spacing[1] = this->Levels[level].DXs[1] / this->LengthScaleFactor;
      spacing[0] = this->Levels[level].DXs[0] / this->LengthScaleFactor;
    }
    else
      {
      spacing[2] = this->Levels[level].DXs[2];
      spacing[1] = this->Levels[level].DXs[1];
      spacing[0] = this->Levels[level].DXs[0];
      }
  }
  else
    {
    if(this->LengthScale)
      {
      spacing[2] = 0.0;
      spacing[1] = this->Levels[level].DXs[1] / this->LengthScaleFactor;
      spacing[0] = this->Levels[level].DXs[0] / this->LengthScaleFactor;
      }
    else
      {
      spacing[2] = 0.0;
      spacing[1] = this->Levels[level].DXs[1];
      spacing[0] = this->Levels[level].DXs[0];
      }
    }
}

int vtkAMRAmazeReaderInternal::GetBlockLevel(const int domain) const
{
  int gid = domain;
  int l=0;

  while(gid >= this->Levels[l].GridsPerLevel)
    {
    gid -= this->Levels[l].GridsPerLevel;
    l++;
    }
  return l;
}

vtkUniformGrid* vtkAMRAmazeReaderInternal::GetAMRGrid(int blockIdx)
{
  int levelId, block_unused;
  // need to retrieve my levelId
  this->FindLevelAndBlock(blockIdx, levelId, block_unused);

  const auto& grid = this->Grids[blockIdx];
  // uniform is uniform, so no need to get all 3 DXs. One is enough

  double dx[3]; // spacing is constant at a given level
  dx[0] = this->Levels[levelId].DXs[0];
  dx[1] = this->Levels[levelId].DXs[1];
  dx[2] = this->Levels[levelId].DXs[2];

  vtkStringArray *sarr = vtkStringArray::New();
  sarr->SetName("GridName");
  sarr->SetNumberOfComponents(1);
  sarr->SetNumberOfTuples(1);
  
  sarr->SetValue(0, std::format("Grid {}", grid.layout.grid_nr));
  vtkUniformGrid* ug = vtkUniformGrid::New();
  ug->Initialize();
  ug->GetFieldData()->AddArray(sarr);
  sarr->Delete();
      
  if(grid.layout.dimensions[2] > 1)  // a 3D grid
    {
    if(this->LengthScale)
       ug->SetSpacing(dx[0]/this->LengthScaleFactor,
                      dx[1]/this->LengthScaleFactor,
                      dx[2]/this->LengthScaleFactor);
    else
      ug->SetSpacing(dx[0], dx[1], dx[2]);
    }
  else
    {
    if(this->LengthScale)
      ug->SetSpacing(dx[0]/this->LengthScaleFactor,
                     dx[1]/this->LengthScaleFactor,
                     0.0);
    else
      ug->SetSpacing(dx[0], dx[1], 0.0);
    }
  if(this->LengthScale)
    ug->SetOrigin(grid.layout.origin[0]/this->LengthScaleFactor,
                  grid.layout.origin[1]/this->LengthScaleFactor,
                  grid.layout.origin[2]/this->LengthScaleFactor);
  else
    ug->SetOrigin(grid.layout.origin[0], grid.layout.origin[1], grid.layout.origin[2]);
  ug->SetDimensions(grid.layout.dimensions[0],
                    grid.layout.dimensions[1],
                    grid.layout.dimensions[2]);

  return ug;
} // GetAMRGrid

// returns the level and block id at that level for the given global id
void  vtkAMRAmazeReaderInternal::FindLevelAndBlock(int domain, int &level, int &block) const
{
  int gid = domain;
  int l=0; 
    
  while(gid >= this->Levels[l].GridsPerLevel)
    {
    gid -= this->Levels[l].GridsPerLevel;
    l++;
    }
  level = l;
  block = gid;
}

vtkDoubleArray* vtkAMRAmazeReaderInternal::ReadVisItVar(int domain, const char *varname)
{
// find which adG_component that is and then go ahead
  int i=0, level, block;
  while(i < this->NumberOfComponents)
    {
    if(!(strcmp(varname, this->Labels[i].label)))
      break;
    i++;
    }
  //cerr << "found label "<< i << ", " << this->Labels[i].label << endl;
  if(i < this->NumberOfComponents) // reached the end without finding it.
    {
    this->FindLevelAndBlock(domain, level, block);
    return this->ReadVar(level, block, this->Labels[i]);
    }
  else
    return nullptr;
}

vtkDoubleArray* vtkAMRAmazeReaderInternal::ReadVar(int levelId, int block, adG_component &variable)
{
  hid_t level_root_id, grid_root_id, dataset_id, mem_space_id;
  int domain = this->FindDomainId(levelId, block);
  /*
  std::cerr << "domain = " << domain
       << ", level = " << levelId
       << ", block = " << block
       << ", varname = " << variable.label
       << endl;
  */
  const auto& grid = this->Grids[domain];

  this->file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  level_root_id = H5Gopen(this->file_id, std::format("/Level {}", levelId).c_str(), H5P_DEFAULT);
  if(level_root_id < 0)
    std::cerr << __LINE__ << ": ReadVar() bad level_root_id returned\n";

  grid_root_id = H5Gopen(level_root_id, std::format("Grid {}", grid.layout.grid_nr).c_str(), H5P_DEFAULT);
  if(grid_root_id < 0)
    std::cerr << __LINE__<< ": ReadVar(): bad grid_root_id returned\n";

// for 2D case with 2D vectors, we create a 3-tuple anyway, fill in the 3rd component with zeroes
// and we must define a hyperslab select to only fill in the first 2 columns.
  vtkDoubleArray*scalars = vtkDoubleArray::New();
  scalars->SetNumberOfComponents(variable.vec_len == 2? 3 : variable.vec_len);
  scalars->SetName((const char*)PVlabels[(const char *)variable.label].c_str());
// default naming. Could be over-written by "Log10()"

  int nvals;
  if(0) //GetCellCentered())
    {
    if(grid.layout.dimensions[2] != 1)
      nvals = (grid.layout.dimensions[0]-1) * (grid.layout.dimensions[1]-1) * (grid.layout.dimensions[2]-1);
    else
      nvals = (grid.layout.dimensions[0]-1) * (grid.layout.dimensions[1]-1);
    }
  else
    nvals = grid.layout.dimensions[0] * grid.layout.dimensions[1] * grid.layout.dimensions[2];
/*
  cerr << lname << ":" << PVlabels[(const char *)variable.label] << "("<< variable.vec_len << "," << nvals << ")\n";
  cerr << "nvals = " << grid.layout.dimensions[0] << "x"<<
                    grid.layout.dimensions[1]<< "x"<<
                    grid.layout.dimensions[2]<< endl;
*/
  hsize_t dims[2], count[2], offset[2];
  dims[0] = nvals;
  dims[1] =  variable.vec_len == 2 ? 3 : variable.vec_len;
  mem_space_id = H5Screate_simple(2,  dims, NULL);

  offset[0] = 0; offset[1] = 0;
  count[0] = nvals; count[1] = variable.vec_len;
  H5Sselect_hyperslab (mem_space_id, H5S_SELECT_SET, offset, NULL, count, NULL);

  scalars->SetNumberOfTuples(nvals);
  void *dataArray = scalars->GetVoidPointer(0);

  dataset_id = H5Dopen(grid_root_id, (const char *) variable.label, H5P_DEFAULT);
  if(dataset_id < 0)
     {
     std::cerr << "error opening HDF5 dataset for var " << variable.label << endl;
     }

  if(H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_space_id, H5S_ALL, H5P_DEFAULT, dataArray) < 0)
  {
    std::cerr << __LINE__ << " :error reading HDF5 dataset for var " << variable.label << endl;
  }

  double *dArray = (double *)dataArray;
  if(variable.vec_len == 2) for(int k=0; k < nvals * 3; k+=3) dArray[k+2] = 0;

  if(this->DataScale == true && variable.scalefactor != 1.0)
    {
    for(int k=0; k < nvals * variable.vec_len; k++)
      {
      dArray[k] /= variable.scalefactor;
      }
    }

  if(this->VarNamesToLog.find((const char *)variable.label) == this->VarNamesToLog.end())
  {
     //std::cout<< variable.label << " is not in the map!"<<endl;
  }
  else if(this->LogData)
  {
    for(int k=0; k < nvals ; k++)
    {
      dArray[k] = log10(dArray[k]);
    }
  }
  H5Sclose(mem_space_id);
  H5Dclose(dataset_id);
  H5Gclose(grid_root_id);
  H5Gclose(level_root_id);
  if(this->file_id)
  {
    H5Fclose(this->file_id);
    this->file_id = 0;
  }
  return scalars;
}

vtkPolyData* vtkAMRAmazeReaderInternal::GetStar(int domain)
{
  return this->Stars[domain];
}

#define THETARES 48
#define PHIRES 24
int vtkAMRAmazeReaderInternal::BuildStars()
{
  hid_t    apr_root_id, dataset1, dataset2, StarsDS, models_root_id;
  hid_t    attr1, attr2, interactions_root_id;
  herr_t   status;
  interaction *interactions = NULL;
  model *models = NULL;

  int  i, j, nb_stars;
  herr_t error;
  this->NumberOfSphericallySymmetricStars = 0;
  this->NumberOfAxisSymmetricStars = 0;
  this->Stars.clear();

// turn off error reporting just in case there are no stars
  H5E_auto2_t func;
  void *client_data;
  H5Eget_auto2(H5E_DEFAULT, &func, &client_data);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  //cerr << "1273: H5Fopen( " << this->FileName << ")\n";
  this->file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  apr_root_id = H5Gopen(this->file_id, "/APR_StellarSystems", H5P_DEFAULT);

  if(apr_root_id < 0)
    {
    return 0;
    }
  interactions_root_id = H5Gopen(apr_root_id, "Interactions", H5P_DEFAULT);
  if(interactions_root_id >=0)
    {
    attr1 = Create_Interaction_Compound();
    interactions = new interaction[2];

    dataset1 = H5Dopen(interactions_root_id, "IsotrInfWind", H5P_DEFAULT);
    if(dataset1 >=0 )
      {
      status = H5Dread(dataset1, attr1, H5S_ALL, H5S_ALL, H5P_DEFAULT, &interactions[0]);
      }
    else
      {
      std::cerr << "error reading Interactions::IsotrInfWind\n";
      }
    dataset2 = H5Dopen(interactions_root_id, "Simple Accretors", H5P_DEFAULT);
    if(dataset2 >=0 )
      {
      status = H5Dread(dataset2, attr1, H5S_ALL, H5S_ALL, H5P_DEFAULT, &interactions[1]);
      }
    else
      {
      std::cerr << "error reading Interactions::Simple Accretors\n";
      }
    H5Dclose(dataset1);
    H5Dclose(dataset2);
    H5Gclose(interactions_root_id);
    }

  models_root_id = H5Gopen(apr_root_id, "Stars", H5P_DEFAULT);

  hid_t dataspace;
  int ModelsArraySize;
  hsize_t dims_out[1]; // we assume rank = 1
  dataset1 = H5Dopen(models_root_id, "Stars: Models&Interaction", H5P_DEFAULT);
  if(dataset1 >=0)
    {
    dataspace = H5Dget_space(dataset1);
    status = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    ModelsArraySize = dims_out[0];
    //cout << "Found " << dims_out[0] << " models of stars\n";
    H5Sclose(dataspace);
    attr2 = Create_StarModel_Compound();

    models = new model[ModelsArraySize];
    status = H5Dread(dataset1, attr2, H5S_ALL, H5S_ALL, H5P_DEFAULT, models);  
    H5Tclose(attr2);
    H5Dclose(dataset1);
    }

  newstar *newstars=NULL;
  star *stars=NULL;
  dataset2 = H5Dopen(models_root_id, "Stars: Present State", H5P_DEFAULT);
  if(dataset2 >=0) // old-style stars before november 6, 2008
    {
    std::cerr << "vtkAMRAmazeReaderInternal::BuildStars(Old-style STARS)\n";
    dataspace = H5Dget_space(dataset2);
    status = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    this->NumberOfSphericallySymmetricStars = nb_stars = dims_out[0];

    stars = new star[nb_stars];
    H5Sclose(dataspace);

    attr2 = Create_Star_Compound();
    status = H5Dread(dataset2, attr2, H5S_ALL, H5S_ALL, H5P_DEFAULT, stars);

    H5Tclose(attr2);
    H5Dclose(dataset2);
    //cout << "Found " << nb_stars << " stars\n";

    for (i=0; i < nb_stars; i++)
      {
      //vtkDoubleArray *velo = vtkDoubleArray::New();
      vtkNew<vtkDoubleArray> velo;
      velo->SetNumberOfComponents(3);
      velo->SetNumberOfTuples(THETARES*(PHIRES-2)+2);
      velo->SetName("Velocity");
      //vtkDoubleArray *mass = vtkDoubleArray::New();
      vtkNew<vtkDoubleArray> mass;
      mass->SetNumberOfComponents(1);
      mass->SetNumberOfTuples(THETARES*(PHIRES-2)+2);
      mass->SetName("Mass");
      //vtkDoubleArray *temp = vtkDoubleArray::New();
      vtkNew<vtkDoubleArray> temp;
      temp->SetNumberOfComponents(1);
      temp->SetNumberOfTuples(THETARES*(PHIRES-2)+2);
      temp->SetName((const char*)PVlabels["Temperature"].c_str());
      //vtkDoubleArray *lum = vtkDoubleArray::New();
      vtkNew<vtkDoubleArray> lum;
      lum->SetNumberOfComponents(1);
      lum->SetNumberOfTuples(THETARES*(PHIRES-2)+2);
      lum->SetName("Luminosity");

      vtkSphereSource *ss = vtkSphereSource::New();
      ss->DebugOff();
      ss->SetCenter(stars[i].Position[0]/this->LengthScaleFactor, stars[i].Position[1]/this->LengthScaleFactor, stars[i].Position[2]/this->LengthScaleFactor);
      ss->SetThetaResolution(THETARES);
      ss->SetPhiResolution(PHIRES);
      stars[i].Radius *= 1; //6.96e10;

      vtkCharArray *name = vtkCharArray::New();
      name->SetNumberOfComponents(1);

      if(strstr(stars[i].Interaction, "Wind")) // look for substring "Wind"
        {
        for(j=0; j < ModelsArraySize; j++)
          {
          if(strstr(models[j].IntActType, "Wind"))
            {
            std::cerr << "found the model: " << models[j].IntActType;
            std::cerr << "multiply by " << interactions[j].CompRadius << endl;
            stars[i].Radius *= interactions[j].CompRadius;
            name->SetNumberOfTuples(strlen("Wind"));
            name->SetName("Wind");
            }
          }
        }
      else
        {
        for(j=0; j < ModelsArraySize; j++)
          {
          if(strstr(models[j].IntActType, "Accretor"))
            {
            std::cerr << "found the model: " << models[j].IntActType;
            std::cerr << "multiply by " << interactions[j].CompRadius << endl;
            stars[i].Radius *= interactions[j].CompRadius;
            name->SetNumberOfTuples(strlen("Accretor"));
            name->SetName("Accretor");
            }
          }
        }
      ss->SetRadius(stars[i].Radius);
      std::cerr << "radius = "<< stars[i].Radius << endl;
      ss->Update();
      for(int k=0; k < THETARES*(PHIRES-2)+2; k++)
        {
        velo->SetTypedTuple(k, stars[i].Velocity);
        mass->SetValue(k, stars[i].Mass);
        if(this->LogData)
          temp->SetValue(k, log10(stars[i].Temperature));
        else
          temp->SetValue(k, stars[i].Temperature);
        lum->SetValue(k, stars[i].Luminosity);
        }
      ss->GetOutput()->GetFieldData()->AddArray(name);
      ss->GetOutput()->GetPointData()->AddArray(velo);
      ss->GetOutput()->GetPointData()->AddArray(mass);
      ss->GetOutput()->GetPointData()->AddArray(temp);
      ss->GetOutput()->GetPointData()->AddArray(lum);
      name->Delete();
      //velo->Delete();
      //mass->Delete();
      //temp->Delete();
      //lum->Delete();
      this->Stars.push_back(ss->GetOutput());
      //ss->Delete();
      }
    } // end of old-style stars before november 6, 2008

  else if((StarsDS = H5Dopen(models_root_id, "Stars", H5P_DEFAULT)) >= 0)
    {
    std::cerr << "\n\n";
    dataspace = H5Dget_space(StarsDS);
    status = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    nb_stars = dims_out[0];

    newstars = new newstar[nb_stars];
    H5Sclose(dataspace);

    attr2 = Create_NewStar_Compound();
    status = H5Dread(StarsDS, attr2, H5S_ALL, H5S_ALL, H5P_DEFAULT, newstars);

    for(i=0; i < nb_stars; i++)
      {
      /*
      std::cout << __LINE__ << ": " << std::format("{:s}\n"
                               //"  StarTime {:e}\n"
                               "  Position {:e}, {:e}, {:e}\n"
                               "  Velocity {:e}, {:e}, {:e}\n"
                               ,
                                   newstars[i].InteractionModel,
                                   //newstars[i].StarTime,
                                   newstars[i].Position[0], newstars[i].Position[1], newstars[i].Position[2],
                                   newstars[i].Velocity[0], newstars[i].Velocity[1], newstars[i].Velocity[2]
                                   ) << std::endl;
     */
/*
      cerr << "star " << i << endl
           << "  StarTime "             << newstars[i].StarTime  << endl
           << "  CompRadiusFrac "       << newstars[i].CompRadiusFrac  << endl
           << "  Mass "                 << newstars[i].Mass  << endl
           << "  SpectralType "         << newstars[i].SpectralType  << endl
           << "  Position "             << newstars[i].Position[0]  << ", " << newstars[i].Position[1]  << ", " << newstars[i].Position[2]  << endl
           << "  Velocity "             << newstars[i].Velocity[0]  << ", " << newstars[i].Velocity[1]  << ", " << newstars[i].Velocity[2]<< endl
           << "  StarModel "            << newstars[i].StarModel  <<endl
           << "  StellarEvolution "     << newstars[i].StellarEvolution  << endl
           << "  StarModelFileName "    << newstars[i].StarModelFileName  << endl
           << "  InteractionModel "     << newstars[i].InteractionModel  << endl
           << "  IActionEvolution "     << newstars[i].IActionEvolution  << endl
           << "  IActionModelFileName " << newstars[i].IActionModelFileName  << endl
           << endl;
*/
      if(strstr(newstars[i].StarModel, "SpherSymStar"))
        NumberOfSphericallySymmetricStars++;
      else
        NumberOfAxisSymmetricStars++;
      }
    H5Tclose(attr2);

    if(this->NumberOfSphericallySymmetricStars)
      {
      hid_t SpherSymStarCurrent_id = Create_SpherSymStar_Compound();
      SpherSymStarCurrent *spherStarData = new SpherSymStarCurrent;
      for (i=0; i < nb_stars; i++)
        {
        double Radius;
        std::string starname = std::format("Star {}", i+1);
        hid_t attr1, dataspace, starN_G, starN_DS;
        hsize_t dims_out[1]; // we assume rank = 1

        starN_G = H5Gopen(models_root_id, starname.c_str(), H5P_DEFAULT);
        if(starN_G >=0)
          {
          starN_DS = H5Dopen(starN_G, "SpherSymStar_Current", H5P_DEFAULT);
          if(starN_DS >=0)
            {
            status = H5Dread(starN_DS, SpherSymStarCurrent_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, spherStarData);
            H5Dclose(starN_DS);
/*
           cerr
           << spherStarData->Radius  << "\t"
           << spherStarData->Temperature  << "\t"
           << spherStarData->Luminosity  << "\t"
           << spherStarData->Omega[0]  << ", " << spherStarData->Omega[1]  << ", " << spherStarData->Omega[2] << "\t"
           << spherStarData->BField[0]  << ", " << spherStarData->BField[1]  << ", " << spherStarData->BField[2]<< endl;
*/
            Radius = spherStarData->Radius;
            }
          else
            Radius=1.0;
          H5Gclose(starN_G);
          }
        else
          {
          std::cerr << "Error opening " << starname << endl;
          }

        vtkNew<vtkDoubleArray> velo;
        velo->SetNumberOfComponents(3);
        velo->SetNumberOfTuples(THETARES*(PHIRES-2)+2);
        velo->SetName("Velocity");
        vtkNew<vtkDoubleArray> mass;
        mass->SetNumberOfComponents(1);
        mass->SetNumberOfTuples(THETARES*(PHIRES-2)+2);
        mass->SetName("Mass");
        char *p = strchr(newstars[i].InteractionModel, ' ');
        *p = '\0';
        vtkCharArray *name = vtkCharArray::New();
        name->SetNumberOfComponents(1);

        name->SetNumberOfTuples(strlen(newstars[i].InteractionModel));
        name->SetName(newstars[i].InteractionModel);

        vtkSphereSource *ss = vtkSphereSource::New();
        ss->SetThetaResolution(THETARES);
        ss->SetPhiResolution(PHIRES);
        if(strstr(newstars[i].StarModel, "SpherSymStar"))
          {
          double c[3];
          if(this->LengthScale)
            {
            ss->SetCenter(newstars[i].Position[0]/this->LengthScaleFactor,
                          newstars[i].Position[1]/this->LengthScaleFactor,
                          newstars[i].Position[2]/this->LengthScaleFactor);
            ss->GetCenter(c);
            }
          else
            {
            ss->SetCenter(newstars[i].Position[0],
                          newstars[i].Position[1],
                          newstars[i].Position[2]);
            ss->GetCenter(c);
            }
          ss->SetRadius(Radius * newstars[i].CompRadiusFrac * 6.96e10 /this->LengthScaleFactor );
          std::cout << __LINE__ << ": " << std::format("{:14s} at Position ({:e}, {:e}, {:e}), Radius = {:e}\n",
                                   newstars[i].InteractionModel, c[0], c[1], c[2], Radius);
          }
        ss->Update();
        for(vtkIdType k=0; k < THETARES*(PHIRES-2)+2; k++)
          {
          velo->SetTypedTuple(k, newstars[i].Velocity);
          mass->SetValue(k, newstars[i].Mass);
          }
        ss->GetOutput()->GetFieldData()->AddArray(name);
        ss->GetOutput()->GetPointData()->AddArray(velo);
        ss->GetOutput()->GetPointData()->AddArray(mass);
        name->Delete();
        this->Stars.push_back(ss->GetOutput());
        //ss->Delete();
        }
      delete spherStarData;
      H5Tclose(SpherSymStarCurrent_id);
      }

  // look if there is any Axi-sym star and create a multi-block holder, then quit the loop

    i=0;
    while(newstars && (i < nb_stars))
      {
      if(strstr(newstars[i].StarModel, "AxiSymStar"))
        {
        //AxiSymStars = vtkMultiBlockDataSet::New();
        //Stars->SetBlock(1, AxiSymStars);
        //AxiSymStars->Delete();
        //Stars->GetMetaData((unsigned int)1)->Set(vtkCompositeDataSet::NAME(), "AxiSymStars");
        }
      i++;
      }

    if(this->NumberOfAxisSymmetricStars) // AxiSymStars)
      {
      hid_t AxiSymStarCurrent_id = Create_AxiSymStar_Compound();

      int J=0; // counter of axis-sym stars
      for(int I=0; I < nb_stars; I++)
        {
        if(strstr(newstars[I].StarModel, "AxiSymStar"))
          {
          std::string starname = std::format("Star {}", J+1);

          hid_t attr1, dataspace, starN_G, starN_DS;
          hsize_t dims_out[1]; // we assume rank = 1
          //std::cerr << "opening " << starname <<"\n";
          starN_G = H5Gopen(models_root_id, starname.c_str(), H5P_DEFAULT);
          AxiSymStarCurrent *axiStarData=NULL;
          int AngleResolution, AxisDirection=0;
          if(starN_G >=0)
            {
            attr1 = H5Aopen_name(starN_G, "AxiSymStar: NTheta");
            status = H5Aread(attr1, H5T_NATIVE_INT, &AngleResolution);
            H5Aclose(attr1);
            attr1 = H5Aopen_name(starN_G, "AxiSymStar: axis-direction");
            if(attr1 >= 0)
              H5Aread(attr1, H5T_NATIVE_INT, &AxisDirection);
            H5Aclose(attr1);

            starN_DS = H5Dopen(starN_G, "AxiSymStar_Current", H5P_DEFAULT);
            dataspace = H5Dget_space(starN_DS);
            status = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
            if(AngleResolution != dims_out[0])
              std::cerr << "sanity check: NTheta resolution is mis-read?\n";

            AngleResolution = dims_out[0];
            //std::cout << "Found phi array of size " << AngleResolution << " for " << starname << endl;
            H5Sclose(dataspace);

            axiStarData = new AxiSymStarCurrent[AngleResolution];
            status = H5Dread(starN_DS, AxiSymStarCurrent_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, axiStarData);
            H5Tclose(AxiSymStarCurrent_id);
            H5Dclose(starN_DS);
            H5Gclose(starN_G);
            }
/*
           cerr << "phi\tRadius\tTemperature\tLuminosity\tEpsilon\tOmega\tBField\n\n";
          for(int k=0; k < AngleResolution; k++)
            {
           cerr
           << axiStarData[k].Theta  << "\t"
           << axiStarData[k].Radius  << "\t"
           << axiStarData[k].Temperature  << "\t"
           << axiStarData[k].Luminosity  << "\t"
           << axiStarData[k].Epsilon  << "\t"
           << axiStarData[k].Omega[0]  << ", " << axiStarData[k].Omega[1]  << ", " << axiStarData[k].Omega[2] << "\t"
           << axiStarData[k].BField[0]  << ", " << axiStarData[k].BField[1]  << ", " << axiStarData[k].BField[2]<< endl;
            }
*/
          vtkPolyData *AxiSymStar = this->AxisSymStarSource(&newstars[I], axiStarData,
                                                      AngleResolution);
          if(AxisDirection == 0 || AxisDirection == 3)
            {
            this->Stars.push_back(AxiSymStar);
            }
          else if(AxisDirection == 1)
            {
            vtkTransform *tf = vtkTransform::New();
            tf->Translate(0.0, 0.0, 0.0);
            tf->Scale(1.0, 1.0, 1.0);
            tf->RotateX(90.0);
            vtkTransformFilter *tfpd = vtkTransformFilter::New();
            tfpd->SetTransform(tf);
            tfpd->SetInputData(AxiSymStar);
            tfpd->Update();
            tf->Delete();
            AxiSymStar->Delete();
            this->Stars.push_back(tfpd->GetPolyDataOutput());
            }
          else if(AxisDirection == 2)
            {
            vtkTransform *tf = vtkTransform::New();
            tf->Translate(newstars[I].Position[0], 0.0, newstars[I].Position[0]);
            tf->Scale(1.0, 1.0, 1.0);
            tf->RotateY(90.0);
            vtkTransformFilter *tfpd = vtkTransformFilter::New();
            tfpd->SetTransform(tf);
            tfpd->SetInputData(AxiSymStar);
            tfpd->Update();
            tf->Delete();
            AxiSymStar->Delete();
            this->Stars.push_back(tfpd->GetPolyDataOutput());
            }
          J++;
          }
        } // do this for all stars
      }
    H5Gclose(StarsDS);
    }
 // end of new-style star after November 6, 2008

  if(interactions)
    delete [] interactions;
  if(models)
    delete [] models;
  if(stars)
    delete [] stars;
  H5Gclose(models_root_id);

  H5Gclose(apr_root_id);
  H5Fclose(this->file_id);
  this->file_id = 0;

  H5Eset_auto2(H5E_DEFAULT, func, client_data);
  //cerr << "exitBuildStars() with " << this->NumberOfSphericallySymmetricStars << " spheric-symmetric stars, and " << this->NumberOfAxisSymmetricStars << " axis-symmetric stars\n";

  return nb_stars;
}

VTK_ABI_NAMESPACE_END
