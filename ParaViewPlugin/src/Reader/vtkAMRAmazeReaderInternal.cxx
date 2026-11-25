#include "vtkAMRAmazeReaderInternal.h"

#include "vtkCharArray.h"
#include "vtkCellArray.h"
#include "vtkDoubleArray.h"
#include "vtkNew.h"
#include "vtkObject.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkSphereSource.h"
#include "vtkTransform.h"
#include "vtkTransformFilter.h"
#include "vtkUniformGrid.h"

#include <cmath>
#include <cstring>
#include <format>
#include <iostream>
#include <vector>

VTK_ABI_NAMESPACE_BEGIN

namespace
{

struct interaction
{
  int    StarNumber;
  int    NTIME;
  int    NAbund;
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
  int  NTime;
  int  NTimePos;
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

hid_t CreateStringType(size_t length)
{
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, length);
  H5Tset_strpad(type, H5T_STR_NULLTERM);
  return type;
}

hid_t CreateDoubleArrayType()
{
  const hsize_t dim = 3;
  return H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, &dim);
}

hid_t Create_Interaction_Compound()
{
  hid_t id = H5Tcreate(H5T_COMPOUND, sizeof(interaction));

  H5Tinsert(id, "Star Number", HOFFSET(interaction, StarNumber), H5T_NATIVE_INT);
  H5Tinsert(id, "NTIME", HOFFSET(interaction, NTIME), H5T_NATIVE_INT);
  H5Tinsert(id, "NAbund", HOFFSET(interaction, NAbund), H5T_NATIVE_INT);
  H5Tinsert(id, "CompRadius", HOFFSET(interaction, CompRadius), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "MassLoss", HOFFSET(interaction, MassLoss), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "VInf", HOFFSET(interaction, VInf), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Temp", HOFFSET(interaction, Temp), H5T_NATIVE_DOUBLE);

  return id;
}

hid_t Create_StarModel_Compound()
{
  hid_t labelstring = CreateStringType(20);

  hid_t id = H5Tcreate(H5T_COMPOUND, sizeof(model));

  H5Tinsert(id, "Model", HOFFSET(model, Model), labelstring);
  H5Tinsert(id, "ModelFileName", HOFFSET(model, ModelFileName), labelstring);
  H5Tinsert(id, "IntActType", HOFFSET(model, IntActType), labelstring);
  H5Tinsert(id, "IntActModel", HOFFSET(model, IntActModel), labelstring);
  H5Tinsert(id, "I-ModelFileName", HOFFSET(model, IModelFileName), labelstring);
  H5Tinsert(id, "NTime", HOFFSET(model, NTime), H5T_NATIVE_INT);
  H5Tinsert(id, "NTimePos", HOFFSET(model, NTimePos), H5T_NATIVE_INT);

  H5Tclose(labelstring);

  return id;
}

hid_t Create_NewStar_Compound()
{
  hid_t labelstring = CreateStringType(100);
  hid_t array1_id = CreateDoubleArrayType();

  hid_t id = H5Tcreate(H5T_COMPOUND, sizeof(newstar));
  H5Tinsert(id, "StarTime [y]", HOFFSET(newstar, StarTime), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "CompRadiusFrac", HOFFSET(newstar, CompRadiusFrac), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Mass [msol]", HOFFSET(newstar, Mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Spectral Type", HOFFSET(newstar, SpectralType), labelstring);
  H5Tinsert(id, "Position", HOFFSET(newstar, Position), array1_id);
  H5Tinsert(id, "Velocity", HOFFSET(newstar, Velocity), array1_id);
  H5Tinsert(id, "Star Model", HOFFSET(newstar, StarModel), labelstring);
  H5Tinsert(id, "Stella rEvolution", HOFFSET(newstar, StellarEvolution), labelstring);
  H5Tinsert(id, "Star Model FileName", HOFFSET(newstar, StarModelFileName), labelstring);
  H5Tinsert(id, "Interaction Model", HOFFSET(newstar, InteractionModel), labelstring);
  H5Tinsert(id, "I-ActionEvolution", HOFFSET(newstar, IActionEvolution), labelstring);
  H5Tinsert(id, "I-ActionModelFileName", HOFFSET(newstar, IActionModelFileName), labelstring);
  H5Tclose(array1_id);
  H5Tclose(labelstring);

  return id;
}

hid_t Create_Star_Compound()
{
  hid_t labelstring = CreateStringType(20);
  hid_t array1_id = CreateDoubleArrayType();

  hid_t id = H5Tcreate(H5T_COMPOUND, sizeof(star));

  H5Tinsert(id, "Type", HOFFSET(star, Type), labelstring);
  H5Tinsert(id, "Position", HOFFSET(star, Position), array1_id);
  H5Tinsert(id, "Velocity", HOFFSET(star, Velocity), array1_id);
  H5Tinsert(id, "Radius", HOFFSET(star, Radius), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Mass", HOFFSET(star, Mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Temperature", HOFFSET(star, Temperature), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Luminosity", HOFFSET(star, Luminosity), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Rotation", HOFFSET(star, Rotation), array1_id);
  H5Tinsert(id, "B-Field", HOFFSET(star, BField), array1_id);
  H5Tinsert(id, "Interaction", HOFFSET(star, Interaction), labelstring);
  H5Tclose(array1_id);
  H5Tclose(labelstring);

  return id;
}

hid_t Create_AxiSymStar_Compound()
{
  const hsize_t dim = 3;
  hid_t array1_id = H5Tarray_create2(H5T_NATIVE_FLOAT, 1, &dim);

  hid_t id = H5Tcreate(H5T_COMPOUND, sizeof(AxiSymStarCurrent));

  H5Tinsert(id, "Theta", HOFFSET(AxiSymStarCurrent, Theta), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Radius", HOFFSET(AxiSymStarCurrent, Radius), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Temperature", HOFFSET(AxiSymStarCurrent, Temperature), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Luminosity", HOFFSET(AxiSymStarCurrent, Luminosity), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Epsilon", HOFFSET(AxiSymStarCurrent, Epsilon), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Omega", HOFFSET(AxiSymStarCurrent, Omega), array1_id);
  H5Tinsert(id, "BField", HOFFSET(AxiSymStarCurrent, BField), array1_id);
  H5Tclose(array1_id);

  return id;
}

hid_t Create_SpherSymStar_Compound()
{
  hid_t array1_id = CreateDoubleArrayType();

  hid_t id = H5Tcreate(H5T_COMPOUND, sizeof(SpherSymStarCurrent));

  H5Tinsert(id, "Radius", HOFFSET(SpherSymStarCurrent, Radius), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Temperature", HOFFSET(SpherSymStarCurrent, Temperature), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Luminosity", HOFFSET(SpherSymStarCurrent, Luminosity), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Omega", HOFFSET(SpherSymStarCurrent, Omega), array1_id);
  H5Tinsert(id, "BField", HOFFSET(SpherSymStarCurrent, BField), array1_id);
  H5Tclose(array1_id);

  return id;
}

} // end anonymous namespace

vtkAMRAmazeReaderInternal::vtkAMRAmazeReaderInternal() = default;

vtkAMRAmazeReaderInternal::~vtkAMRAmazeReaderInternal()
{
  this->SetFileName(nullptr);

  for (auto* starPtr : this->Stars)
  {
    starPtr->Delete();
  }
  this->Stars.clear();
}

void vtkAMRAmazeReaderInternal::SetFileName(const char* fileName)
{
  this->FileName = fileName ? fileName : "";

  if (this->file_id > 0)
  {
    H5Fclose(this->file_id);
    this->file_id = 0;
  }
}

/* ---------------------------------------------------------
   returns the number of stars, or a negative value on error
 ---------------------------------------------------------*/
int vtkAMRAmazeReaderInternal::ReadMetaData()
{
  this->Levels.clear();
  this->Labels.clear();

  if (this->FileName.empty())
  {
    return -1;
  }

  this->file_id = H5Fopen(this->FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (this->file_id < 0)
  {
    std::cerr << "file could not be opened. Check filename " << this->FileName << std::endl;
    return -1;
  }

  int nbstars = this->ReadHDF5MetaData();

  this->Levels.resize(this->NumberOfLevels);
  this->Labels.resize(this->NumberOfComponents);
  this->ReadHDF5VariablesMetaData();

  this->Grids.resize(this->NumberOfGrids);
  this->ReadHDF5GridsMetaData();

  if (this->file_id > 0)
  {
    H5Fclose(this->file_id);
    this->file_id = 0;
  }

  // compute the number of grids per level and build the AMR boxes
  int current_level = -1;
  for (int i = 0; i < this->NumberOfGrids; i++)
  {
    this->Grids[i].amrbox =
      vtkAMRBox(this->Grids[i].layout.box_corners, this->Grids[i].layout.box_corners + 3);
    if (this->Grids[i].layout.level != current_level)
    { // first grid of a given level
      current_level = this->Grids[i].layout.level;
      this->Levels[current_level].GridsPerLevel = 1;
    }
    else
    {
      this->Levels[current_level].GridsPerLevel++;
    }
  }

  return nbstars;
}

constexpr int MAX_FLOATS = 24;

// Data passed to the attribute iterator callback below
struct AttributeData
{
  double values[MAX_FLOATS];
  int    count;
};

// Callback function executed for every attribute found in the group
static herr_t attr_iterator(hid_t g_id, const char* attr_name, const H5A_info_t*, void* op_data)
{
  auto* data = static_cast<AttributeData*>(op_data);

  // Stop iterating if we already found MAX_FLOATS floats
  if (data->count >= MAX_FLOATS)
  {
    return 1; // a positive value stops the iteration successfully
  }

  hid_t attr_id = H5Aopen(g_id, attr_name, H5P_DEFAULT);
  if (attr_id < 0)
  {
    std::cerr << "Warning: Failed to open attribute " << attr_name << "\n";
    return 0; // skip errors and continue iterating
  }

  double val = 0.0;
  if (H5Aread(attr_id, H5T_NATIVE_DOUBLE, &val) >= 0)
  {
    data->values[data->count++] = val;
  }

  H5Aclose(attr_id);
  return 0; // continue to the next attribute
}

void vtkAMRAmazeReaderInternal::ReadHDF5GridsMetaData()
{
  hid_t root_id = H5Gopen(this->file_id, "/", H5P_DEFAULT);
  hid_t dataset = H5Dopen(root_id, "Grid Info", H5P_DEFAULT);
  if (dataset < 0)
  {
    std::cerr << "failed to find Grid Info. Returning without action\n";
    H5Gclose(root_id);
    return;
  }

  hid_t adG_grid_id = H5Tcreate(H5T_COMPOUND, sizeof(adG_grid));

  H5Tinsert(adG_grid_id, "grid number", HOFFSET(adG_grid_layout, grid_nr), H5T_NATIVE_INT);
  H5Tinsert(adG_grid_id, "level", HOFFSET(adG_grid_layout, level), H5T_NATIVE_INT);

  hsize_t dim = adG_MAXDIM;
  hid_t array0_id = H5Tarray_create2(H5T_NATIVE_INT, 1, &dim);
  dim = 2 * adG_MAXDIM;
  hid_t array1_id = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, &dim);
  hid_t array2_id = H5Tarray_create2(H5T_NATIVE_INT, 1, &dim);

  H5Tinsert(adG_grid_id, "dimensions", HOFFSET(adG_grid_layout, dimensions), array0_id);
  H5Tinsert(adG_grid_id, "origin", HOFFSET(adG_grid_layout, origin), array1_id);
  H5Tinsert(adG_grid_id, "box corners", HOFFSET(adG_grid_layout, box_corners), array2_id);

  H5Tclose(array0_id);
  H5Tclose(array1_id);
  H5Tclose(array2_id);

  if (H5Dread(dataset, adG_grid_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->Grids[0]) < 0)
  {
    std::cerr << "error reading Grid Info\n";
  }
  H5Tclose(adG_grid_id);
  H5Dclose(dataset);

  for (int levelId = 0; levelId < this->NumberOfLevels; levelId++)
  {
    std::string name = std::format("Level {}", levelId);

    hid_t attr1 = H5Aopen_by_name(root_id, name.c_str(), "Refinement Ratio", H5P_DEFAULT, H5P_DEFAULT);
    H5Aread(attr1, H5T_NATIVE_INT, &this->Levels[levelId].RefRatio);
    H5Aclose(attr1);

    attr1 = H5Aopen_by_name(root_id, name.c_str(), "Cell Edge Length", H5P_DEFAULT, H5P_DEFAULT);
    H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->Levels[levelId].DXs[0]);
    H5Aclose(attr1);

    for (int axis = 1; axis < 3; axis++)
    {
      std::string attrName = std::format("Cell Edge Length {}", axis == 1 ? 'Y' : 'Z');
      attr1 = H5Aexists_by_name(root_id, name.c_str(), attrName.c_str(), H5P_DEFAULT) > 0
        ? H5Aopen_by_name(root_id, name.c_str(), attrName.c_str(), H5P_DEFAULT, H5P_DEFAULT)
        : -1;
      if (attr1 >= 0)
      {
        H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->Levels[levelId].DXs[axis]);
        H5Aclose(attr1);
      }
      else
      {
        this->Levels[levelId].DXs[axis] = this->Levels[levelId].DXs[0];
      }
    }
  }

  // Now read the Stellar Orbit data
  if (H5Lexists(root_id, "APR_StellarSystems", H5P_DEFAULT) > 0)
  {
    hid_t apr_root_id = H5Gopen(root_id, "APR_StellarSystems", H5P_DEFAULT);
    if (apr_root_id >= 0)
    {
      hid_t stars_id = H5Gopen(apr_root_id, "Stars", H5P_DEFAULT);
      if (stars_id >= 0)
      {
        hid_t stellar_orbit_id =
          H5Lexists(stars_id, "Stellar Orbit", H5P_DEFAULT) > 0
            ? H5Gopen(stars_id, "Stellar Orbit", H5P_DEFAULT)
            : -1;
        if (stellar_orbit_id >= 0)
        {
          hid_t fix_binary_orbit_id =
            H5Lexists(stellar_orbit_id, "Fixed Binary Orbit", H5P_DEFAULT) > 0
              ? H5Gopen(stellar_orbit_id, "Fixed Binary Orbit", H5P_DEFAULT)
              : -1;
          if (fix_binary_orbit_id >= 0)
          {
            AttributeData attr_data;
            attr_data.count = 0;

            // Iterate through all attributes in the group.
            // H5_INDEX_NAME: sorts alphabetically.
            hsize_t idx = 0;
            H5Aiterate2(fix_binary_orbit_id, H5_INDEX_NAME, H5_ITER_INC, &idx, attr_iterator,
                        &attr_data);

            if (attr_data.count < MAX_FLOATS)
            {
              std::cerr << "Warning: Expected " << MAX_FLOATS
                        << " attributes, but only found " << attr_data.count << ".\n";
            }
            std::memcpy(&this->FixedBinaryOrbit, attr_data.values, sizeof(Fixed_Binary_Orbit));
            H5Gclose(fix_binary_orbit_id);
          }
          H5Gclose(stellar_orbit_id);
        }
        H5Gclose(stars_id);
      }
      H5Gclose(apr_root_id);
    }
  } // Done reading APR_StellarSystems

  if (H5Lexists(root_id, "Thermodynamics", H5P_DEFAULT) > 0)
  {
    hid_t thermodynamics_id = H5Gopen(root_id, "Thermodynamics", H5P_DEFAULT);
    if (thermodynamics_id >= 0)
    {
      hid_t mem_type = H5Tcreate(H5T_COMPOUND, sizeof(this->thermodynamics_eos));
      H5Tinsert(mem_type, "Adiabatic Exponent", HOFFSET(Thermodynamics_EOS, Adiabatic_Exponent),
                H5T_NATIVE_DOUBLE);
      H5Tinsert(mem_type, "Specific Gasconstant", HOFFSET(Thermodynamics_EOS, Specific_Gasconstant),
                H5T_NATIVE_DOUBLE);
      H5Tinsert(mem_type, "Scaled Specific Gasconstant",
                HOFFSET(Thermodynamics_EOS, Scaled_Specific_Gasconstant), H5T_NATIVE_DOUBLE);
      H5Tinsert(mem_type, "Scale Factor", HOFFSET(Thermodynamics_EOS, Scale_Factor), H5T_NATIVE_DOUBLE);
      H5Tinsert(mem_type, "Mean Particle Mass", HOFFSET(Thermodynamics_EOS, Mean_Particle_Mass),
                H5T_NATIVE_DOUBLE);

      const char* eosName =
        H5Lexists(thermodynamics_id, "Polytropic EOS", H5P_DEFAULT) > 0 ? "Polytropic EOS"
        : H5Lexists(thermodynamics_id, "Ideal Gas EOS", H5P_DEFAULT) > 0 ? "Ideal Gas EOS"
                                                                         : nullptr;
      if (eosName)
      {
        dataset = H5Dopen2(thermodynamics_id, eosName, H5P_DEFAULT);
        if (dataset >= 0)
        {
          if (H5Dread(dataset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->thermodynamics_eos) < 0)
          {
            std::cerr << "Error reading " << eosName << "\n";
          }
          H5Dclose(dataset);
        }
      }
      H5Tclose(mem_type);
      H5Gclose(thermodynamics_id);
    }
  } // Done reading Thermodynamics stuff if present

  H5Gclose(root_id);

  if (this->MappedGrids != MapName::NoMap)
  {
    root_id = H5Gopen(this->file_id, "/Map", H5P_DEFAULT);
    if (root_id >= 0)
    {
      dataset = H5Dopen(root_id, "Map Parameter", H5P_DEFAULT);
      if (dataset < 0)
      {
        std::cerr << "error opening Map_Parameter\n";
      }
      hid_t mapping_id = -1;
      switch (this->MappedGrids)
      {
        case MapName::Sphere_LogR:
        {
          hid_t label1 = CreateStringType(15);

          mapping_id = H5Tcreate(H5T_COMPOUND, sizeof(ParamToPhysicalMapping));

          H5Tinsert(mapping_id, "Space Name", HOFFSET(ParamToPhysicalMapping, label), label1);
          H5Tinsert(mapping_id, "Min Radius", HOFFSET(ParamToPhysicalMapping, Rmin), H5T_NATIVE_DOUBLE);
          H5Tinsert(mapping_id, "Max Radius", HOFFSET(ParamToPhysicalMapping, Rmax), H5T_NATIVE_DOUBLE);
          H5Tinsert(mapping_id, "Min Theta", HOFFSET(ParamToPhysicalMapping, Tmin), H5T_NATIVE_DOUBLE);
          H5Tinsert(mapping_id, "Max Theta", HOFFSET(ParamToPhysicalMapping, Tmax), H5T_NATIVE_DOUBLE);
          H5Tinsert(mapping_id, "Min Phi", HOFFSET(ParamToPhysicalMapping, Pmin), H5T_NATIVE_DOUBLE);
          H5Tinsert(mapping_id, "Max Phi", HOFFSET(ParamToPhysicalMapping, Pmax), H5T_NATIVE_DOUBLE);
          H5Tinsert(mapping_id, "Stretch Factor", HOFFSET(ParamToPhysicalMapping, StretchFactorOBSOLETE),
                    H5T_NATIVE_DOUBLE);
          H5Tclose(label1);

          if (H5Dread(dataset, mapping_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->SphereLogRMappings) < 0)
          {
            std::cerr << "error reading Map_Parameter\n";
          }
          break;
        }

        case MapName::DCR_Cart2Spheres:
        {
          hid_t label1 = CreateStringType(6);
          hid_t label2 = CreateStringType(7);

          mapping_id = H5Tcreate(H5T_COMPOUND, sizeof(DCR_Mapping));

          H5Tinsert(mapping_id, "Min Radius Physical Space", HOFFSET(DCR_Mapping, Rmin), H5T_NATIVE_DOUBLE);
          H5Tinsert(mapping_id, "Max Radius Physical Space", HOFFSET(DCR_Mapping, Rmax), H5T_NATIVE_DOUBLE);
          H5Tinsert(mapping_id, "Map Case", HOFFSET(DCR_Mapping, MapCase), label1);
          H5Tinsert(mapping_id, "Map Lunarity", HOFFSET(DCR_Mapping, MapLunarity), label2);
          H5Tinsert(mapping_id, "Dimensionality", HOFFSET(DCR_Mapping, Dimension), H5T_NATIVE_INT);
          H5Tclose(label1);
          H5Tclose(label2);

          if (H5Dread(dataset, mapping_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->DCR_Mappings) < 0)
          {
            std::cerr << "error reading Map_Parameter\n";
          }
          break;
        }

        default:
          std::cerr << "error finding an implemented Mapping code\n";
          break;
      }
      if (mapping_id >= 0)
      {
        H5Tclose(mapping_id);
      }
      H5Dclose(dataset);
      H5Gclose(root_id);
    }
  }
} // ReadHDF5GridsMetaData()

void vtkAMRAmazeReaderInternal::ReadHDF5VariablesMetaData()
{
  hid_t root_id = H5Gopen(this->file_id, "/", H5P_DEFAULT);
  hid_t dataset = H5Dopen(root_id, "Variable Info", H5P_DEFAULT);

  hid_t labelstring = CreateStringType(adG_LABELLENGTH);
  hid_t unitstring = CreateStringType(adG_UNITLENGTH);

  hid_t adG_component_id = H5Tcreate(H5T_COMPOUND, sizeof(adG_component));

  H5Tinsert(adG_component_id, "vector length", HOFFSET(adG_component, vec_len), H5T_NATIVE_INT);
  H5Tinsert(adG_component_id, "Variable Name", HOFFSET(adG_component, label), labelstring);
  H5Tinsert(adG_component_id, "Variable Unit", HOFFSET(adG_component, unit), unitstring);
  H5Tinsert(adG_component_id, "scale factor", HOFFSET(adG_component, scalefactor), H5T_NATIVE_FLOAT);

  H5Tclose(labelstring);
  H5Tclose(unitstring);

  if (H5Dread(dataset, adG_component_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->Labels[0]) < 0)
  {
    std::cerr << "error reading Variable Info\n";
  }

  H5Tclose(adG_component_id);
  H5Dclose(dataset);
  H5Gclose(root_id);
}

// build the ParaView display names ("label [unit]") for all variables
void vtkAMRAmazeReaderInternal::MakeVariableNames()
{
  for (int c = 0; c < this->NumberOfComponents; c++)
  {
    std::string unit = this->Labels[c].unit;
    std::string label = this->Labels[c].label;

    if (unit.empty())
    { // unit-less variables such as Mach Number
      this->PVlabels[label] = label;
    }
    else
    {
      this->PVlabels[label] = std::format("{} [{}]", label, unit);
    }
  }
}

/* ---------------------------------------------------------
   returns the number of stars
 ---------------------------------------------------------*/
int vtkAMRAmazeReaderInternal::ReadHDF5MetaData()
{
  hid_t root_id = H5Gopen(this->file_id, "/", H5P_DEFAULT);

  hid_t attr1 = H5Aopen_name(root_id, "NumberLevels");
  H5Aread(attr1, H5T_NATIVE_INT, &this->NumberOfLevels);
  H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "NumberOfComponents");
  H5Aread(attr1, H5T_NATIVE_INT, &this->NumberOfComponents);
  H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "NumberOfGrids");
  H5Aread(attr1, H5T_NATIVE_INT, &this->NumberOfGrids);
  H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "Dimensionality");
  H5Aread(attr1, H5T_NATIVE_INT, &this->Dimensionality);
  H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "Time");
  if (H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->AMAZETime) < 0)
  {
    vtkWarningWithObjectMacro(nullptr, "Failed to open Time Attribute");
  }
  H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "Time Scaling Factor");
  if (H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->AMAZETimeScalor) >= 0)
  {
    this->AMAZETime /= this->AMAZETimeScalor;
  }
  else
  {
    vtkWarningWithObjectMacro(nullptr, "Failed to open Time Scaling Factor Attribute");
  }
  H5Aclose(attr1);

  if (H5Aexists(root_id, "Length Scale Factor") > 0)
  {
    attr1 = H5Aopen_name(root_id, "Length Scale Factor");
    if (attr1 >= 0)
    {
      H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->LengthScaleFactor);
      H5Aclose(attr1);
    }
    else
    {
      this->LengthScaleFactor = 1.0;
    }
  }

  switch (this->ScaleChoice)
  {
    case ScaleType::pc:
      this->LengthScaleFactor *= 3.08567782e18;
      std::cout << "!!!\nUsing PARSEC with Length Scale Factor * by 3.08567782e18 = "
                << this->LengthScaleFactor << "!!!\n";
      break;
    case ScaleType::AU:
      this->LengthScaleFactor *= 1.49597870700e13;
      std::cout << "!!!\nUsing AU with Length Scale Factor * by 1.49597870700e13 = "
                << this->LengthScaleFactor << "!!!\n";
      break;
    case ScaleType::RSun:
      this->LengthScaleFactor *= 6.96342e10;
      std::cout << "!!!\nUsing RSun with Length Scale Factor * by 6.96342e10 = "
                << this->LengthScaleFactor << "!!!\n";
      break;
    case ScaleType::NoScale:
      break;
  }

  if (H5Lexists(root_id, "/Map", H5P_DEFAULT))
  {
    hid_t map_id = H5Gopen(root_id, "/Map", H5P_DEFAULT);
    attr1 = H5Aopen_name(map_id, "Map Type");
    hid_t t1 = H5Aget_type(attr1);
    char map_type[30] = "";
    if (H5Aread(attr1, t1, map_type) >= 0)
    {
      if (!strncmp(map_type, "Sphere-LogR", 11))
      {
        this->MappedGrids = MapName::Sphere_LogR;
      }
      else if (!strncmp(map_type, "DCR_Cart2Spheres", 16))
      {
        this->MappedGrids = MapName::DCR_Cart2Spheres;
      }
    }
    H5Tclose(t1);
    H5Aclose(attr1);
    H5Gclose(map_id);
  }

  int nb_stars = 0;

  if (H5Lexists(root_id, "/APR_StellarSystems", H5P_DEFAULT))
  {
    hid_t apr_root_id = H5Gopen(root_id, "/APR_StellarSystems", H5P_DEFAULT);
    hid_t models_root_id = H5Gopen(apr_root_id, "Stars", H5P_DEFAULT);
    hsize_t dims_out[1]; // we assume rank = 1

    hid_t StarsDS = H5Dopen(models_root_id, "Stars", H5P_DEFAULT);
    if (StarsDS >= 0)
    {
      hid_t dataspace = H5Dget_space(StarsDS);
      H5Sget_simple_extent_dims(dataspace, dims_out, nullptr);
      nb_stars = dims_out[0];
      H5Sclose(dataspace);
      H5Dclose(StarsDS);
    }

    H5Gclose(models_root_id);
    H5Gclose(apr_root_id);
  }
  H5Gclose(root_id);
  return nb_stars;
} // end of ReadHDF5MetaData()

// returns the global id for the given block
int vtkAMRAmazeReaderInternal::FindDomainId(int level, int block)
{
  int domain = block;
  for (int l = 0; l < level; l++)
  {
    domain += this->Levels[l].GridsPerLevel;
  }
  if (domain >= this->NumberOfGrids)
  {
    std::cerr << "level or block too high for this dataset\n";
  }
  return domain;
}

void vtkAMRAmazeReaderInternal::GetSpacing(int level, double* spacing)
{
  const double divisor = this->LengthScale ? this->LengthScaleFactor : 1.0;
  spacing[0] = this->Levels[level].DXs[0] / divisor;
  spacing[1] = this->Levels[level].DXs[1] / divisor;
  spacing[2] = this->Dimensionality == 3 ? this->Levels[level].DXs[2] / divisor : 0.0;
}

int vtkAMRAmazeReaderInternal::GetBlockLevel(int domain) const
{
  int gid = domain;
  int l = 0;

  while (gid >= this->Levels[l].GridsPerLevel)
  {
    gid -= this->Levels[l].GridsPerLevel;
    l++;
  }
  return l;
}

// returns the level and block id at that level for the given global id
void vtkAMRAmazeReaderInternal::FindLevelAndBlock(int domain, int& level, int& block) const
{
  int gid = domain;
  int l = 0;

  while (gid >= this->Levels[l].GridsPerLevel)
  {
    gid -= this->Levels[l].GridsPerLevel;
    l++;
  }
  level = l;
  block = gid;
}

vtkUniformGrid* vtkAMRAmazeReaderInternal::GetAMRGrid(int blockIdx)
{
  // need to retrieve my levelId
  int levelId, block_unused;
  this->FindLevelAndBlock(blockIdx, levelId, block_unused);

  const auto& grid = this->Grids[blockIdx];
  // spacing is constant at a given level
  const double divisor = this->LengthScale ? this->LengthScaleFactor : 1.0;

  double spacing[3];
  spacing[0] = this->Levels[levelId].DXs[0] / divisor;
  spacing[1] = this->Levels[levelId].DXs[1] / divisor;
  spacing[2] = grid.layout.dimensions[2] > 1 ? this->Levels[levelId].DXs[2] / divisor : 0.0;

  vtkUniformGrid* ug = vtkUniformGrid::New();
  ug->Initialize();
  ug->SetSpacing(spacing[0], spacing[1], spacing[2]);
  ug->SetOrigin(grid.layout.origin[0] / divisor, grid.layout.origin[1] / divisor,
                grid.layout.origin[2] / divisor);
  ug->SetDimensions(grid.layout.dimensions[0], grid.layout.dimensions[1], grid.layout.dimensions[2]);

  return ug;
} // GetAMRGrid

vtkDoubleArray* vtkAMRAmazeReaderInternal::ReadVisItVar(int domain, const char* varname)
{
  // find which adG_component that is and then go ahead
  int i = 0;
  while (i < this->NumberOfComponents && strcmp(varname, this->Labels[i].label) != 0)
  {
    i++;
  }
  if (i == this->NumberOfComponents) // reached the end without finding it
  {
    return nullptr;
  }

  int level, block;
  this->FindLevelAndBlock(domain, level, block);
  return this->ReadVar(level, block, this->Labels[i]);
}

vtkDoubleArray* vtkAMRAmazeReaderInternal::ReadVar(int levelId, int block, adG_component& variable)
{
  const int domain = this->FindDomainId(levelId, block);
  const auto& grid = this->Grids[domain];

  this->file_id = H5Fopen(this->FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t level_root_id = H5Gopen(this->file_id, std::format("/Level {}", levelId).c_str(), H5P_DEFAULT);
  if (level_root_id < 0)
  {
    std::cerr << __LINE__ << ": ReadVar() bad level_root_id returned\n";
  }

  hid_t grid_root_id =
    H5Gopen(level_root_id, std::format("Grid {}", grid.layout.grid_nr).c_str(), H5P_DEFAULT);
  if (grid_root_id < 0)
  {
    std::cerr << __LINE__ << ": ReadVar(): bad grid_root_id returned\n";
  }

  vtkDoubleArray* scalars = vtkDoubleArray::New();
  // for the 2D case with 2D vectors, we create a 3-tuple anyway, fill in the
  // 3rd component with zeroes, and we must define a hyperslab select to only
  // fill in the first 2 columns.
  const int ncomp = variable.vec_len == 2 ? 3 : variable.vec_len;
  scalars->SetNumberOfComponents(ncomp);
  scalars->SetName(this->PVlabels[variable.label].c_str());

  const int nvals = grid.layout.dimensions[0] * grid.layout.dimensions[1] * grid.layout.dimensions[2];

  hsize_t dims[2] = { static_cast<hsize_t>(nvals), static_cast<hsize_t>(ncomp) };
  hid_t mem_space_id = H5Screate_simple(2, dims, nullptr);

  hsize_t offset[2] = { 0, 0 };
  hsize_t count[2] = { static_cast<hsize_t>(nvals), static_cast<hsize_t>(variable.vec_len) };
  H5Sselect_hyperslab(mem_space_id, H5S_SELECT_SET, offset, nullptr, count, nullptr);

  scalars->SetNumberOfTuples(nvals);
  double* dArray = static_cast<double*>(scalars->GetVoidPointer(0));

  hid_t dataset_id = H5Dopen(grid_root_id, variable.label, H5P_DEFAULT);
  if (dataset_id < 0)
  {
    std::cerr << "error opening HDF5 dataset for var " << variable.label << std::endl;
  }

  if (H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_space_id, H5S_ALL, H5P_DEFAULT, dArray) < 0)
  {
    std::cerr << __LINE__ << " :error reading HDF5 dataset for var " << variable.label << std::endl;
  }

  if (variable.vec_len == 2)
  {
    for (int k = 0; k < nvals * 3; k += 3)
    {
      dArray[k + 2] = 0;
    }
  }

  if (this->DataScale && variable.scalefactor != 1.0)
  {
    for (int k = 0; k < nvals * variable.vec_len; k++)
    {
      dArray[k] /= variable.scalefactor;
    }
  }

  H5Sclose(mem_space_id);
  H5Dclose(dataset_id);
  H5Gclose(grid_root_id);
  H5Gclose(level_root_id);
  if (this->file_id)
  {
    H5Fclose(this->file_id);
    this->file_id = 0;
  }
  return scalars;
}

constexpr int THETARES = 48;
constexpr int PHIRES = 24;

int vtkAMRAmazeReaderInternal::BuildStars()
{
  this->NumberOfSphericallySymmetricStars = 0;
  this->NumberOfAxisSymmetricStars = 0;
  for (auto* starPtr : this->Stars)
  {
    starPtr->Delete();
  }
  this->Stars.clear();
  this->newstars.clear();

  int nb_stars = 0;

  // turn off error reporting just in case there are no stars
  H5E_auto2_t func = nullptr;
  void* client_data = nullptr;
  H5Eget_auto2(H5E_DEFAULT, &func, &client_data);
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);

  this->file_id = H5Fopen(this->FileName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t apr_root_id = H5Gopen(this->file_id, "/APR_StellarSystems", H5P_DEFAULT);
  if (apr_root_id < 0)
  {
    H5Fclose(this->file_id);
    this->file_id = 0;
    H5Eset_auto2(H5E_DEFAULT, func, client_data);
    return 0;
  }

  std::vector<interaction> interactions(2);
  std::vector<model> models;
  std::vector<star> stars;
  int ModelsArraySize = 0;

  hid_t interactions_root_id = H5Gopen(apr_root_id, "Interactions", H5P_DEFAULT);
  if (interactions_root_id >= 0)
  {
    hid_t compound_id = Create_Interaction_Compound();

    hid_t dataset1 = H5Dopen(interactions_root_id, "IsotrInfWind", H5P_DEFAULT);
    if (dataset1 >= 0)
    {
      H5Dread(dataset1, compound_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &interactions[0]);
      H5Dclose(dataset1);
    }
    else
    {
      std::cerr << "error reading Interactions::IsotrInfWind\n";
    }
    hid_t dataset2 = H5Dopen(interactions_root_id, "Simple Accretors", H5P_DEFAULT);
    if (dataset2 >= 0)
    {
      H5Dread(dataset2, compound_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &interactions[1]);
      H5Dclose(dataset2);
    }
    else
    {
      std::cerr << "error reading Interactions::Simple Accretors\n";
    }
    H5Tclose(compound_id);
    H5Gclose(interactions_root_id);
  }

  hid_t models_root_id = H5Gopen(apr_root_id, "Stars", H5P_DEFAULT);

  hsize_t dims_out[1]; // we assume rank = 1
  hid_t dataset1 = H5Dopen(models_root_id, "Stars: Models&Interaction", H5P_DEFAULT);
  if (dataset1 >= 0)
  {
    hid_t dataspace = H5Dget_space(dataset1);
    H5Sget_simple_extent_dims(dataspace, dims_out, nullptr);
    ModelsArraySize = dims_out[0];
    H5Sclose(dataspace);

    hid_t compound_id = Create_StarModel_Compound();
    models.resize(ModelsArraySize);
    H5Dread(dataset1, compound_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, models.data());
    H5Tclose(compound_id);
    H5Dclose(dataset1);
  }

  hid_t dataset2 = H5Dopen(models_root_id, "Stars: Present State", H5P_DEFAULT);
  if (dataset2 >= 0) // old-style stars before november 6, 2008
  {
    std::cerr << "vtkAMRAmazeReaderInternal::BuildStars(Old-style STARS)\n";
    hid_t dataspace = H5Dget_space(dataset2);
    H5Sget_simple_extent_dims(dataspace, dims_out, nullptr);
    this->NumberOfSphericallySymmetricStars = nb_stars = dims_out[0];
    H5Sclose(dataspace);

    stars.resize(nb_stars);
    hid_t compound_id = Create_Star_Compound();
    H5Dread(dataset2, compound_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, stars.data());
    H5Tclose(compound_id);
    H5Dclose(dataset2);

    for (int i = 0; i < nb_stars; i++)
    {
      vtkNew<vtkDoubleArray> velo;
      velo->SetNumberOfComponents(3);
      velo->SetNumberOfTuples(THETARES * (PHIRES - 2) + 2);
      velo->SetName("Velocity");
      vtkNew<vtkDoubleArray> mass;
      mass->SetNumberOfComponents(1);
      mass->SetNumberOfTuples(THETARES * (PHIRES - 2) + 2);
      mass->SetName("Mass");
      vtkNew<vtkDoubleArray> temp;
      temp->SetNumberOfComponents(1);
      temp->SetNumberOfTuples(THETARES * (PHIRES - 2) + 2);
      temp->SetName(this->PVlabels["Temperature"].c_str());
      vtkNew<vtkDoubleArray> lum;
      lum->SetNumberOfComponents(1);
      lum->SetNumberOfTuples(THETARES * (PHIRES - 2) + 2);
      lum->SetName("Luminosity");

      vtkSphereSource* ss = vtkSphereSource::New();
      ss->DebugOff();
      ss->SetCenter(stars[i].Position[0] / this->LengthScaleFactor,
                    stars[i].Position[1] / this->LengthScaleFactor,
                    stars[i].Position[2] / this->LengthScaleFactor);
      ss->SetThetaResolution(THETARES);
      ss->SetPhiResolution(PHIRES);

      vtkCharArray* name = vtkCharArray::New();
      name->SetNumberOfComponents(1);

      if (strstr(stars[i].Interaction, "Wind")) // look for substring "Wind"
      {
        for (int j = 0; j < ModelsArraySize; j++)
        {
          if (strstr(models[j].IntActType, "Wind"))
          {
            stars[i].Radius *= interactions[j].CompRadius;
            name->SetNumberOfTuples(strlen("Wind"));
            name->SetName("Wind");
          }
        }
      }
      else
      {
        for (int j = 0; j < ModelsArraySize; j++)
        {
          if (strstr(models[j].IntActType, "Accretor"))
          {
            stars[i].Radius *= interactions[j].CompRadius;
            name->SetNumberOfTuples(strlen("Accretor"));
            name->SetName("Accretor");
          }
        }
      }
      ss->SetRadius(stars[i].Radius);
      ss->Update();
      for (int k = 0; k < THETARES * (PHIRES - 2) + 2; k++)
      {
        velo->SetTypedTuple(k, stars[i].Velocity);
        mass->SetValue(k, stars[i].Mass);
        temp->SetValue(k, stars[i].Temperature);
        lum->SetValue(k, stars[i].Luminosity);
      }
      ss->GetOutput()->GetFieldData()->AddArray(name);
      ss->GetOutput()->GetPointData()->AddArray(velo);
      ss->GetOutput()->GetPointData()->AddArray(mass);
      ss->GetOutput()->GetPointData()->AddArray(temp);
      ss->GetOutput()->GetPointData()->AddArray(lum);
      name->Delete();

      vtkPolyData* starPolyData = ss->GetOutput();
      starPolyData->Register(nullptr); // keep alive after the source is deleted
      this->Stars.push_back(starPolyData);
      ss->Delete();
    }
  } // end of old-style stars before november 6, 2008

  else if (hid_t StarsDS = H5Dopen(models_root_id, "Stars", H5P_DEFAULT); StarsDS >= 0)
  {
    // new-style stars after November 6, 2008
    hid_t dataspace = H5Dget_space(StarsDS);
    H5Sget_simple_extent_dims(dataspace, dims_out, nullptr);
    nb_stars = dims_out[0];
    H5Sclose(dataspace);

    this->newstars.resize(nb_stars);
    hid_t compound_id = Create_NewStar_Compound();
    H5Dread(StarsDS, compound_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, this->newstars.data());

    for (int i = 0; i < nb_stars; i++)
    {
      if (strstr(this->newstars[i].StarModel, "SpherSymStar"))
      {
        this->NumberOfSphericallySymmetricStars++;
      }
      else
      {
        this->NumberOfAxisSymmetricStars++;
      }
    }
    H5Tclose(compound_id);

    if (this->NumberOfSphericallySymmetricStars)
    {
      hid_t SpherSymStarCurrent_id = Create_SpherSymStar_Compound();
      SpherSymStarCurrent spherStarData;
      for (int i = 0; i < nb_stars; i++)
      {
        double Radius = 1.0;
        std::string starname = std::format("Star {}", i + 1);

        hid_t starN_G = H5Gopen(models_root_id, starname.c_str(), H5P_DEFAULT);
        if (starN_G >= 0)
        {
          hid_t starN_DS = H5Dopen(starN_G, "SpherSymStar_Current", H5P_DEFAULT);
          if (starN_DS >= 0)
          {
            H5Dread(starN_DS, SpherSymStarCurrent_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &spherStarData);
            Radius = spherStarData.Radius;
            H5Dclose(starN_DS);
          }
          H5Gclose(starN_G);
        }
        else
        {
          std::cerr << "Error opening " << starname << std::endl;
        }

        vtkNew<vtkDoubleArray> velo;
        velo->SetNumberOfComponents(3);
        velo->SetNumberOfTuples(THETARES * (PHIRES - 2) + 2);
        velo->SetName("Velocity");
        vtkNew<vtkDoubleArray> mass;
        mass->SetNumberOfComponents(1);
        mass->SetNumberOfTuples(THETARES * (PHIRES - 2) + 2);
        mass->SetName("Mass");

        char* p = strchr(this->newstars[i].InteractionModel, ' ');
        if (p)
        {
          *p = '\0';
        }
        vtkCharArray* name = vtkCharArray::New();
        name->SetNumberOfComponents(1);
        name->SetNumberOfTuples(strlen(this->newstars[i].InteractionModel));
        name->SetName(this->newstars[i].InteractionModel);

        vtkSphereSource* ss = vtkSphereSource::New();
        ss->SetThetaResolution(THETARES);
        ss->SetPhiResolution(PHIRES);
        if (strstr(this->newstars[i].StarModel, "SpherSymStar"))
        {
          double center[3];
          if (this->LengthScale)
          {
            ss->SetCenter(this->newstars[i].Position[0] / this->LengthScaleFactor,
                          this->newstars[i].Position[1] / this->LengthScaleFactor,
                          this->newstars[i].Position[2] / this->LengthScaleFactor);
          }
          else
          {
            ss->SetCenter(this->newstars[i].Position[0], this->newstars[i].Position[1],
                          this->newstars[i].Position[2]);
          }
          ss->GetCenter(center);
          ss->SetRadius(Radius * this->newstars[i].CompRadiusFrac * 6.96e10 / this->LengthScaleFactor);
          std::cout << __LINE__ << ": " << std::format("{:14s} at Position ({:e}, {:e}, {:e}), Radius = {:e}\n",
                                   newstars[i].InteractionModel, center[0], center[1], center[2], Radius);
        }
        ss->Update();
        for (vtkIdType k = 0; k < THETARES * (PHIRES - 2) + 2; k++)
        {
          velo->SetTypedTuple(k, this->newstars[i].Velocity);
          mass->SetValue(k, this->newstars[i].Mass);
        }
        ss->GetOutput()->GetFieldData()->AddArray(name);
        ss->GetOutput()->GetPointData()->AddArray(velo);
        ss->GetOutput()->GetPointData()->AddArray(mass);
        name->Delete();

        vtkPolyData* starPolyData = ss->GetOutput();
        starPolyData->Register(nullptr); // keep alive after the source is deleted
        this->Stars.push_back(starPolyData);
        ss->Delete();
      }
      H5Tclose(SpherSymStarCurrent_id);
    }

    if (this->NumberOfAxisSymmetricStars)
    {
      hid_t AxiSymStarCurrent_id = Create_AxiSymStar_Compound();

      int J = 0; // counter of axis-sym stars
      for (int I = 0; I < nb_stars; I++)
      {
        if (!strstr(this->newstars[I].StarModel, "AxiSymStar"))
        {
          continue;
        }
        std::string starname = std::format("Star {}", J + 1);

        hid_t starN_G = H5Gopen(models_root_id, starname.c_str(), H5P_DEFAULT);
        std::vector<AxiSymStarCurrent> axiStarData;
        int AngleResolution = 0, AxisDirection = 0;
        if (starN_G >= 0)
        {
          hid_t attr = H5Aopen_name(starN_G, "AxiSymStar: NTheta");
          H5Aread(attr, H5T_NATIVE_INT, &AngleResolution);
          H5Aclose(attr);
          attr = H5Aopen_name(starN_G, "AxiSymStar: axis-direction");
          if (attr >= 0)
          {
            H5Aread(attr, H5T_NATIVE_INT, &AxisDirection);
            H5Aclose(attr);
          }

          hid_t starN_DS = H5Dopen(starN_G, "AxiSymStar_Current", H5P_DEFAULT);
          hid_t space = H5Dget_space(starN_DS);
          H5Sget_simple_extent_dims(space, dims_out, nullptr);
          if (AngleResolution != static_cast<int>(dims_out[0]))
          {
            std::cerr << "sanity check: NTheta resolution is mis-read?\n";
          }
          AngleResolution = dims_out[0];
          H5Sclose(space);

          axiStarData.resize(AngleResolution);
          H5Dread(starN_DS, AxiSymStarCurrent_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, axiStarData.data());
          H5Dclose(starN_DS);
          H5Gclose(starN_G);
        }

        vtkPolyData* AxiSymStar =
          this->AxisSymStarSource(&this->newstars[I], axiStarData.data(), AngleResolution);
        if (AxisDirection == 0 || AxisDirection == 3)
        {
          this->Stars.push_back(AxiSymStar);
        }
        else
        {
          vtkNew<vtkTransform> tf;
          if (AxisDirection == 1)
          {
            tf->RotateX(90.0);
          }
          else // AxisDirection == 2
          {
            tf->Translate(this->newstars[I].Position[0], 0.0, this->newstars[I].Position[0]);
            tf->RotateY(90.0);
          }
          vtkNew<vtkTransformFilter> tfpd;
          tfpd->SetTransform(tf.Get());
          tfpd->SetInputData(AxiSymStar);
          tfpd->Update();
          AxiSymStar->Delete();
          vtkPolyData* starPolyData = tfpd->GetPolyDataOutput();
          starPolyData->Register(nullptr); // keep alive after the filter is deleted
          this->Stars.push_back(starPolyData);
        }
        J++;
      } // do this for all stars
      H5Tclose(AxiSymStarCurrent_id);
    }
    H5Dclose(StarsDS);
  } // end of new-style star after November 6, 2008

  H5Gclose(models_root_id);
  H5Gclose(apr_root_id);
  H5Fclose(this->file_id);
  this->file_id = 0;

  H5Eset_auto2(H5E_DEFAULT, func, client_data);

  return nb_stars;
}

vtkPolyData* vtkAMRAmazeReaderInternal::AxisSymStarSource(newstar* astar,
                                                          AxiSymStarCurrent* axiStarData,
                                                          int AngleResolution)
{
  vtkPolyData* AxiSymStar = vtkPolyData::New();

  double Center[3];
  double Radius;
  if (this->LengthScale)
  {
    Radius = astar->CompRadiusFrac / this->LengthScaleFactor;
  }
  else
  {
    Radius = astar->CompRadiusFrac;
  }
  Center[0] = astar->Position[0];
  Center[1] = astar->Position[1];
  Center[2] = astar->Position[2];

  const int thetaResolution = AngleResolution - 2;
  const int phiResolution = AngleResolution + 1;

  const int numPts = thetaResolution * phiResolution + 2;

  vtkPoints* newPoints = vtkPoints::New();
  newPoints->SetDataTypeToDouble();
  newPoints->Initialize();
  newPoints->Reserve(numPts);

  vtkCellArray* newPolys = vtkCellArray::New();
  newPolys->Initialize();

  double x[3], deltaphi, theta, phi, radius;

  // Create north pole
  x[0] = Center[0];
  x[1] = Center[1];
  x[2] = Center[2] + Radius * axiStarData[0].Radius;
  newPoints->InsertNextPoint(x);

  // Create south pole
  x[0] = Center[0];
  x[1] = Center[1];
  x[2] = Center[2] - Radius * axiStarData[AngleResolution - 1].Radius;
  newPoints->InsertNextPoint(x);

  const int numPoles = 2;
  deltaphi = (2.0 * M_PI) / phiResolution;
  // Create intermediate points, looping around the Z axis
  for (int i = 0; i < phiResolution; i++)
  {
    phi = i * deltaphi;
    for (int j = 1; j < AngleResolution - 1; j++)
    {
      theta = M_PI_2 - axiStarData[j].Theta; // theta should range between 0 and M_PI
      radius = Radius * axiStarData[j].Radius;
      x[0] = Center[0] + radius * cos(phi) * sin(theta);
      x[1] = Center[1] + radius * sin(phi) * sin(theta);
      x[2] = Center[2] + radius * cos(theta);
      newPoints->InsertNextPoint(x);
    }
  }

  // Generate mesh connectivity
  const int base = thetaResolution * phiResolution;

  // around north pole
  for (int i = 0; i < phiResolution; i++)
  {
    newPolys->InsertNextCell({ thetaResolution * i + numPoles,
                               (thetaResolution * (i + 1) % base) + numPoles, 0 });
  }

  // around south pole
  const int numOffset = thetaResolution - 1 + numPoles;
  for (int i = 0; i < phiResolution; i++)
  {
    newPolys->InsertNextCell({ thetaResolution * i + numOffset, numPoles - 1,
                               ((thetaResolution * (i + 1)) % base) + numOffset });
  }

  // bands in-between poles
  for (int i = 0; i < phiResolution; i++)
  {
    for (int j = 0; j < (thetaResolution - 1); j++)
    {
      const vtkIdType p0 = thetaResolution * i + j + numPoles;
      const vtkIdType p1 = p0 + 1;
      const vtkIdType p2 = ((thetaResolution * (i + 1) + j) % base) + numPoles + 1;
      newPolys->InsertNextCell({ p0, p1, p2 });
      newPolys->InsertNextCell({ p0, p2, p2 - 1 });
    }
  }

  newPoints->Squeeze();
  AxiSymStar->SetPoints(newPoints);
  newPoints->Delete();

  newPolys->Squeeze();
  AxiSymStar->SetPolys(newPolys);
  newPolys->Delete();

  // now the data
  vtkDoubleArray* Temperature = vtkDoubleArray::New();
  Temperature->SetNumberOfComponents(1);
  Temperature->SetNumberOfTuples((AngleResolution - 2) * (AngleResolution + 1) + 2);
  Temperature->SetName(this->PVlabels["Temperature"].c_str());
  AxiSymStar->GetPointData()->AddArray(Temperature);
  int p0 = 0;
  Temperature->SetValue(p0++, axiStarData[0].Temperature);
  Temperature->SetValue(p0++, axiStarData[AngleResolution - 1].Temperature);
  for (int p1 = 0; p1 < (AngleResolution + 1); p1++)
  {
    for (int p2 = 1; p2 < AngleResolution - 1; p2++)
    {
      Temperature->SetValue(p0++, axiStarData[p2].Temperature);
    }
  }
  Temperature->Delete();

  vtkNew<vtkDoubleArray> mass;
  mass->SetNumberOfComponents(1);
  mass->SetNumberOfTuples(1);
  mass->SetName("Mass");
  vtkNew<vtkDoubleArray> velo;
  velo->SetNumberOfComponents(3);
  velo->SetNumberOfTuples(1);
  velo->SetName("Velocity");

  AxiSymStar->GetFieldData()->AddArray(velo);
  AxiSymStar->GetFieldData()->AddArray(mass);

  mass->SetValue(0, astar->Mass);
  velo->SetTypedTuple(0, astar->Velocity);

  return AxiSymStar;
}

VTK_ABI_NAMESPACE_END
