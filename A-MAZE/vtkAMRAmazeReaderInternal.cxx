#include "vtkAMRAmazeReaderInternal.h"

#include "vtkAMRBox.h"
#include "vtkByteSwap.h"
#include "vtkStringArray.h"
#include "vtkDataArray.h"
#include "vtkCellArray.h"
#include "vtkDataSetAttributes.h"
#include "vtkDataSetWriter.h"
#include "vtkDoubleArray.h"
#include "vtkSphereSource.h"
#include "vtkFloatArray.h"
#include "vtkErrorCode.h"

#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkStructuredGrid.h"
#include "vtkRectilinearGrid.h"
#include "vtkUniformGrid.h"
#include <stdio.h>
#include <math.h>
#include "vtkTimerLog.h"

//#define SINGLE_OUTPUT_PORT 1
//#define PARALLEL_DEBUG 1
#include <stddef.h>

#include <vector>
#include <string>
#include <map>
#include <cassert>
#include <sstream>

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
  array1_id = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim, NULL);

  id = H5Tcreate(H5T_COMPOUND, sizeof(newstar));
  H5Tinsert(id, "StarTime [y]",             HOFFSET(newstar, StarTime), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "CompRadiusFrac",       HOFFSET(newstar, CompRadiusFrac), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Mass [msol]",                 HOFFSET(newstar, Mass), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Spectral Type",         HOFFSET(newstar, SpectralType), labelstring);
  H5Tinsert(id, "Position",             HOFFSET(newstar, Position), array1_id);
  H5Tinsert(id, "Velocity",             HOFFSET(newstar, Velocity), array1_id);
  H5Tinsert(id, "Star Model",            HOFFSET(newstar, StarModel), labelstring);
  H5Tinsert(id, "Stella rEvolution",     HOFFSET(newstar, StellarEvolution), labelstring);
  H5Tinsert(id, "Star Model FileName",    HOFFSET(newstar, StarModelFileName), labelstring);
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
  array1_id = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim, NULL);

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
  hsize_t array1_id = H5Tarray_create(H5T_NATIVE_FLOAT, 1, dim, NULL);

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
  hsize_t array1_id = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim, NULL);

  hid_t   id = H5Tcreate(H5T_COMPOUND, sizeof(SpherSymStarCurrent));

  H5Tinsert(id, "Radius",      HOFFSET(SpherSymStarCurrent, Radius),      H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Temperature", HOFFSET(SpherSymStarCurrent, Temperature), H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Luminosity",  HOFFSET(SpherSymStarCurrent, Luminosity),  H5T_NATIVE_DOUBLE);
  H5Tinsert(id, "Omega",       HOFFSET(SpherSymStarCurrent, Omega),       array1_id);
  H5Tinsert(id, "BField",      HOFFSET(SpherSymStarCurrent, BField),      array1_id);
  H5Tclose(array1_id);

  return id;
}

vtkPolyData * vtkAMRAmazeReaderInternal::AxisSymStarSource(newstar *astar,
                                      struct AxiSymStarCurrent *axiStarData,
                                      int AngleResolution)
{
  vtkPolyData *AxiSymStar = vtkPolyData::New();

  double Radius, Center[3];

  if(this->LengthScale)
    {
    Radius = astar->CompRadiusFrac / this->LengthScaleFactor;
    //Radius = astar->CompRadiusFrac * 6.96e10 / this->LengthScaleFactor;  
    Center[0] = astar->Position[0];// / this->LengthScaleFactor;
    Center[1] = astar->Position[1];// / this->LengthScaleFactor;
    Center[2] = astar->Position[2];// / this->LengthScaleFactor;
    }
  else
    {
    Radius = astar->CompRadiusFrac;
    //Radius = astar->CompRadiusFrac * 6.96e10;
    Center[0] = astar->Position[0];
    Center[1] = astar->Position[1];
    Center[2] = astar->Position[2];
    }

  double Starttheta = 0.0;  // Rolf uses PI/2 for the North Pole, 
  double Endtheta = 180.0;  //      and -PI/2 for South Pole
  int LatLongTessellation = 0;
  int numPts, numPolys;
  vtkPoints *newPoints;
  vtkFloatArray *newNormals;
  vtkCellArray *newPolys;
  double x[3], n[3], deltaphi, theta, phi, radius, norm;
  int base, numPoles=0;
  vtkIdType pts[4];

  int thetaResolution = AngleResolution - 2;
  int phiResolution = AngleResolution + 1;

  numPts = thetaResolution * phiResolution + 2;
  // creating triangles
  numPolys = thetaResolution * 2 * phiResolution;
  //cerr << "  numPts "             << numPts  << ", numPolys " << numPolys  << endl;
  newPoints = vtkPoints::New();
  newPoints->SetDataTypeToDouble();
  newPoints->Allocate(numPts);
  newNormals = vtkFloatArray::New();
  newNormals->SetNumberOfComponents(3);
  newNormals->Allocate(3*numPts);
  newNormals->SetName("Normals");

  newPolys = vtkCellArray::New();
  newPolys->Allocate(newPolys->EstimateSize(numPolys, 3));

  //cerr << "  Center "             << Center[0]  << ", " << Center[1]  << ", " << Center[2]  << endl;

  // Create north pole
  x[0] = Center[0];
  x[1] = Center[1];
  x[2] = Center[2] + Radius * axiStarData[0].Radius;
  newPoints->InsertNextPoint(x);
  //cerr << "North Pole at: "  << x[0]  << ", " << x[1]  << ", " << x[2]  << endl;
  x[0] = x[1] = 0.0; x[2] = 1.0;
  newNormals->InsertNextTuple(x);
  numPoles++;

    // Create south pole
  x[0] = Center[0];
  x[1] = Center[1];
  x[2] = Center[2] - Radius * axiStarData[AngleResolution-1].Radius;
  newPoints->InsertNextPoint(x);
  //cerr << "South Pole at: "  << x[0]  << ", " << x[1]  << ", " << x[2]  << endl;
  x[0] = x[1] = 0.0; x[2] = -1.0;
  newNormals->InsertNextTuple(x);
  numPoles++;

  deltaphi = (2.0 * M_PI) / phiResolution;
  //cerr << "loop "  << phiResolution << " around the Z axis"<< endl;
  // Create intermediate points
  for (int i=0; i < phiResolution; i++)
    {
    phi = i * deltaphi;
    //cerr << i << " phi = " << phi << endl;
    for (int j=1; j< AngleResolution-1; j++)
      {
      theta = M_PI_2 - axiStarData[j].Theta; // theta should range between 0 and M_PI
      //if(i==0) cerr << theta << endl;
      radius = Radius * axiStarData[j].Radius;
      n[0] = radius * cos((double)phi) * sin((double)theta);
      n[1] = radius * sin((double)phi)* sin((double)theta);
      n[2] = radius * cos((double)theta);
      x[0] = n[0] + Center[0];
      x[1] = n[1] + Center[1];
      x[2] = n[2] + Center[2];
      newPoints->InsertNextPoint(x);

      if ( (norm = vtkMath::Norm(n)) == 0.0 )
        {
        norm = 1.0;
        }
      n[0] /= norm; n[1] /= norm; n[2] /= norm;
      newNormals->InsertNextTuple(n);
      }
    }

   // Generate mesh connect  H5Eset_auto(func, client_data);ivity
   base = thetaResolution * phiResolution;

   // around north pole
  for (int i=0; i < phiResolution; i++)
   {
   pts[0] = thetaResolution*i + numPoles;
   pts[1] = (thetaResolution*(i+1) % base) + numPoles;
   pts[2] = 0;
   newPolys->InsertNextCell(3, pts);
   }

// around south pole
   int numOffset = thetaResolution - 1 + numPoles;

  for (int i=0; i < phiResolution; i++)
    {
    pts[0] = thetaResolution*i + numOffset;
    pts[2] = ((thetaResolution*(i+1)) % base) + numOffset;
    pts[1] = numPoles - 1;
    newPolys->InsertNextCell(3, pts);
    }

  // bands in-between poles
  for (int i=0; i < phiResolution; i++)
    {
    for (int j=0; j < (thetaResolution-1); j++)
      {
      pts[0] = thetaResolution*i + j + numPoles;
      pts[1] = pts[0] + 1;
      pts[2] = ((thetaResolution*(i+1)+j) % base) + numPoles + 1;
      if ( !LatLongTessellation )
         {
         newPolys->InsertNextCell(3, pts);
         pts[1] = pts[2];
         pts[2] = pts[1] - 1;
         newPolys->InsertNextCell(3, pts);
         }
      else
        {
        pts[3] = pts[2] - 1;
        newPolys->InsertNextCell(4, pts);
        }
      }
    }

  newPoints->Squeeze();
  AxiSymStar->SetPoints(newPoints);
  newPoints->Delete();

  newNormals->Squeeze();
  //AxiSymStar->GetPointData()->SetNormals(newNormals);
  newNormals->Delete();

  newPolys->Squeeze();
  AxiSymStar->SetPolys(newPolys);
  newPolys->Delete();

// now the data
  vtkDoubleArray *Temperature = vtkDoubleArray::New();
  Temperature->SetNumberOfComponents(1);
  Temperature->SetNumberOfTuples((AngleResolution-2) * (AngleResolution + 1) + 2);
  Temperature->SetName("Temperature");
  Temperature->SetName((const char*)PVlabels["Temperature"].c_str());
  AxiSymStar->GetPointData()->AddArray(Temperature);
  int p0=0;
  Temperature->SetValue(p0++, axiStarData[0].Temperature);
  Temperature->SetValue(p0++, axiStarData[AngleResolution-1].Temperature);
  for(int p1=0; p1 < (AngleResolution + 1); p1++)
    {
    for(int p2=1; p2 < AngleResolution-1; p2++)
      {
      Temperature->SetValue(p0++, axiStarData[p2].Temperature);
      }
    }
  if(this->LogData)
    {
     for(int p1=0; p1 < Temperature->GetNumberOfTuples(); p1++)
       {
       Temperature->SetValue(p1, log10(Temperature->GetValue(p1)));
       }
    }
  Temperature->Delete();

  vtkDoubleArray *mass = vtkDoubleArray::New();
  mass->SetNumberOfComponents(1);
  mass->SetNumberOfTuples(1);
  mass->SetName("Mass");
  vtkDoubleArray *velo = vtkDoubleArray::New();
  velo->SetNumberOfComponents(3);
  velo->SetNumberOfTuples(1);
  velo->SetName("Velocity");

  AxiSymStar->GetFieldData()->AddArray(velo);
  AxiSymStar->GetFieldData()->AddArray(mass);

  mass->SetValue(0, astar->Mass);
  velo->SetTypedTuple(0, astar->Velocity);
  velo->Delete();
  mass->Delete();

  return AxiSymStar;
}

//vtkCxxRevisionMacro(vtkAMRAmazeReaderInternal, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkAMRAmazeReaderInternal);

vtkAMRAmazeReaderInternal::vtkAMRAmazeReaderInternal()
{
  this->FileName = NULL;
  this->file_id = 0;
  this->Dimensionality = 0;
  this->NumberOfLevels = 0;
  this->NumberOfGrids = 0;
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
  this->LengthScale = 1; // GUI will ALWAYS overwrite that value
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
  for(i=0; i < this->Stars.size(); i++)
    {
    //(this->Stars[i])->Delete();
    }
  this->Stars.clear();
  
  this->Labels.clear();
  
  this->Grids.clear();
  
  this->VarNamesToLog.clear();
  CloseHDF5File("~vtkAMRAmazeReaderInternal");
}

void vtkAMRAmazeReaderInternal::SetFileName(char * fileName )
{
  this->FileName = fileName;
  //cerr << __LINE__ << " Settting this->Filename = " << fileName << "\n";
}

//----------------------------------------------------------------------------
void vtkAMRAmazeReaderInternal::ReadMetaData()
{
  // Check to see if we have read it
  if ( this->NumberOfLevels > 0 )
    {
    return;
    }
  int levelId, GridId, i, node_veclen;
  double time, time_scalor;

  this->Levels.clear();
  this->Labels.clear();
  if (!this->FileName)
    {
    //this->SetErrorCode(vtkErrorCode::NoFileNameError);
    //vtkErrorMacro(<< "Must specify adG file");
    cerr << __LINE__ << " this->FileName is NULL " << ")\n";
    return;
    }
  OpenHDF5File("ReadMetaData");
  this->ReadHDF5MetaData();
  //cerr << "done with ReadHDF5MetaData() " << endl;
  this->SetTime(this->GetTime() / this->TimeScalor);

  //info->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &this->Time, 1);

  this->Levels.resize(this->NumberOfLevels);

  this->Labels.resize(this->NumberOfComponents);

  this->ReadHDF5VariablesMetaData();
  //cerr << "done with ReadHDF5VariablesMetaData() " << endl;
  this->Grids.resize(this->NumberOfGrids);

  this->ReadHDF5GridsMetaData(false);
  //cerr << "done with ReadHDF5GridsMetaData() " << endl;

  this->CheckVarSize(0, 0, this->Labels[0]);

  CloseHDF5File("ReadMetaData");

///////////////////////////////////////////////////////////////////
  int current_level = -1;
  std::vector<adG_grid> &grid = this->Grids;

  for (i = 0; i < this->NumberOfGrids; i++)
    {
    this->Grids[i].amrbox = vtkAMRBox(grid[i].box_corners, grid[i].box_corners+3);
    if(this->Grids[i].level != current_level)
      { // first grid of a given level
      current_level = this->Grids[i].level;
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

  cerr<< "\tDimensionality: " << this->Dimensionality << "\n\tNumberOfComponents: " << this->NumberOfComponents << "\n\tNumberOfLevels: " << this->NumberOfLevels << "\n\tNumberOfGrids: " << this->NumberOfGrids << endl;
 grid = this->Grids;
 int size = grid[0].dimensions[0];

  for (i = 0; i < this->NumberOfLevels; i++)
    {
    cerr<< "\t\tLevel " << i << " has " << this->Levels[i].GridsPerLevel << " grids" << endl;
    }
  return;
}

void vtkAMRAmazeReaderInternal::WriteVTIFile(const char *fileName, adG_grid &grid, double *spacing)
{
  assert( "fileName should not be NULL" && (fileName != NULL) );
  FILE *fp = fopen(fileName, "w");
  fprintf(fp, "<VTKFile type=\"ImageData\" version=\"2.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "  <ImageData WholeExtent=\"0 %d 0 %d 0 %d\" Origin=\"%g %g %g\" Spacing=\"%g %g %g\">\n", grid.dimensions[0]-1, grid.dimensions[1]-1, grid.dimensions[2]-1, grid.origin[0], grid.origin[1], grid.origin[2], spacing[0], spacing[1], spacing[2]);
  fprintf(fp, "    <Piece Extent=\"0 %d 0 %d 0 %d\">\n", grid.dimensions[0]-1, grid.dimensions[1]-1, grid.dimensions[2]-1);
  fprintf(fp, "      <PointData>\n");
  fprintf(fp, "      </PointData>\n");
  fprintf(fp, "      <CellData>\n");
  fprintf(fp, "      </CellData>\n");
  fprintf(fp, "    </Piece>\n");
  fprintf(fp, "  </ImageData>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

void vtkAMRAmazeReaderInternal::WriteXMLFile(const char *fileName)
{
  assert( "fileName should not be NULL" && (fileName != NULL) );
  FILE *fp = fopen(fileName, "w");
  fprintf(fp, "<VTKFile type=\"vtkHierarchicalBoxDataSet\" version=\"1.1\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp, "  <vtkHierarchicalBoxDataSet origin=\"0 0 0\" grid_description=\"XYZ\">\n");
  int I=0;
  for (int i = 0; i < 3;i++) //this->NumberOfLevels; i++)
    {
    double spacing[3];
    this->GetSpacing(i, spacing);
    //cerr<< "\t\tLevel " << i << " has " << this->Levels[i].GridsPerLevel << " grids" << endl;
    fprintf(fp,   "    <Block level=\"%d\" spacing=\"%g %g %g\">\n", i, spacing[0], spacing[1], spacing[2]);
    for (int d = 0; d < this->Levels[i].GridsPerLevel; d++, I++)
      {
      std::ostringstream vtifile;
      vtifile << "/tmp/amr/amr_" << I << ".vti" << ends;
      WriteVTIFile(vtifile.str().c_str(), this->Grids[I], spacing);
      vtkAMRBox &box = this->Grids[i].amrbox;
      int *corners = this->Grids[I].box_corners;
      int lo[3] = {box.GetLoCorner()[0], box.GetLoCorner()[1], box.GetLoCorner()[2]};
      int hi[3] = {box.GetHiCorner()[0], box.GetHiCorner()[1], box.GetHiCorner()[2]};
      fprintf(fp,   "      <DataSet index=\"%d\" amr_box=\"%d %d %d %d %d %d\" file=\"%s\">\n", d, corners[0], corners[3], corners[1], corners[4], corners[2], corners[5], vtifile.str().c_str());
      fprintf(fp,   "      </DataSet>\n");
      }
    fprintf(fp,   "    </Block>\n");
    }
  fprintf(fp, "  </vtkHierarchicalBoxDataSet>\n");
  fprintf(fp, "</VTKFile>\n");
  fclose(fp);
}

int vtkAMRAmazeReaderInternal::OpenHDF5File(const char *funcName)
{
  if(!this->file_id) // not already opened
    {
    this->file_id = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
    if(this->file_id < 0)
      {
      cerr << __LINE__ << "file could not be opened. Check filename = " << this->FileName << ")\n";
      this->file_id = 0;
      return -1;
      }
    //cerr << funcName << "::H5Fopen(" << this->FileName << ")\n";
    }
  return 0;
}

void vtkAMRAmazeReaderInternal::CloseHDF5File(const char *funcName)
{
  if(this->file_id)
    {
    H5Fclose(this->file_id);
    this->file_id = 0;
    //cerr << funcName << "::H5Fclose()\n";
  }
}

/*
Added the Logical to Physical mappers def.
May 23, 2011
*/
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
      OpenHDF5File("ReadHDF5GridsMetaData");
      }
    //cout << "reading the shifted grid info\n";
    root_id = H5Gopen(this->file_id, "/");
    dataset = H5Dopen(root_id, "Shifted Grid Info");
    if(dataset < 0)
      {
      cerr << "failed to find shifted grid info. Returning without action\n";
      H5Gclose(root_id);
      return;
      }
    CloseHDF5File(" ReadHDF5GridsMetaData");
    }
  else
    {
    root_id = H5Gopen(this->file_id, "/");
    dataset = H5Dopen(root_id, "Grid Info");
    }
  adG_grid_id = H5Tcreate (H5T_COMPOUND, sizeof(adG_grid));

  H5Tinsert(adG_grid_id, "grid number", offsetof(adG_grid, grid_nr), H5T_NATIVE_INT);
  H5Tinsert(adG_grid_id, "level", HOFFSET(adG_grid, level), H5T_NATIVE_INT);

  dim[0] = adG_MAXDIM;
  array0_id = H5Tarray_create(H5T_NATIVE_INT, 1, dim, NULL);

  dim[0] = 2 * adG_MAXDIM;
  array1_id = H5Tarray_create(H5T_NATIVE_DOUBLE, 1, dim, NULL);

  dim[0] = 2 * adG_MAXDIM;
  array2_id = H5Tarray_create(H5T_NATIVE_INT, 1, dim, NULL);

  H5Tinsert(adG_grid_id, "dimensions", HOFFSET(adG_grid, dimensions),
            array0_id);
  H5Tinsert(adG_grid_id, "origin", HOFFSET(adG_grid, origin), array1_id);
  H5Tinsert(adG_grid_id, "box corners", HOFFSET(adG_grid, box_corners), array2_id);

  H5Tclose(array0_id);
  H5Tclose(array1_id);
  H5Tclose(array2_id);

  status = H5Dread(dataset, adG_grid_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, &this->Grids[0]);
  H5Tclose(adG_grid_id);
  H5Dclose(dataset);

  for (int levelId=0; levelId < this->NumberOfLevels; levelId++)
    {
    char name[16];
    sprintf(name, "Level %d", levelId);

    //level_root_id = H5Gopen(root_id, name); //levelName.str());
    //if(level_root_id < 0)
      //cerr << "bad level_root_id returned\n";

    attr1 = H5Aopen_by_name(root_id, name, "Refinement Ratio", H5P_DEFAULT, H5P_DEFAULT);
    status = H5Aread(attr1, H5T_NATIVE_INT, &this->Levels[levelId].RefRatio);
    status = H5Aclose(attr1);

    attr1 = H5Aopen_by_name(root_id, name, "Cell Edge Length", H5P_DEFAULT, H5P_DEFAULT);
    status = H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->Levels[levelId].DXs[0]);
    status = H5Aclose(attr1);

    if(H5Aexists_by_name(root_id, name, "Cell Edge Length Y", H5P_DEFAULT) > 0)
      {
      attr1 = H5Aopen_by_name(root_id, name, "Cell Edge Length Y", H5P_DEFAULT, H5P_DEFAULT);
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

    if(H5Aexists_by_name(root_id, name, "Cell Edge Length Z", H5P_DEFAULT) > 0)
      {
      attr1 = H5Aopen_by_name(root_id, name, "Cell Edge Length Z", H5P_DEFAULT, H5P_DEFAULT);
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

  if(this->MappedGrids && (root_id = H5Gopen(this->file_id, "/Map") )>= 0)
    {
    dataset = H5Dopen(root_id, "Map Parameter");
    if(dataset < 0)
      {
      cerr << "error opening Map_Parameter\n";
      }
    switch(this->MappedGrids) {
      case Sphere_LogR:
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
        cerr << "error reading Map_Parameter\n";
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
      
      case DCR_Cart2Spheres:
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
        cerr << "error reading Map_Parameter\n";
        }
      cerr << "\nparam =\t"  << "\n\t";
      cerr << DCR_Mappings.Rmin<< "\n\t";
      cerr << DCR_Mappings.Rmax << "\n\t";
      cerr << DCR_Mappings.MapCase << "\n\t";
      cerr << DCR_Mappings.MapLunarity << "\n\t";
      cerr << DCR_Mappings.Dimension << "\n\n";
      break;
      
      default:
        cerr << "error finding an implemented Mapping code\n";
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

  root_id = H5Gopen(this->file_id, "/");
  dataset = H5Dopen(root_id, "Variable Info");

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
  //cerr << __LINE__ << " MakeVariableNames()\n";
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

void vtkAMRAmazeReaderInternal::ReadHDF5MetaData()
{
  hid_t   root_id;
  hid_t   attr1;
  herr_t  status;
// turn off error reporting
  H5E_auto_t func;
  void *client_data;
  H5Eget_auto(&func, &client_data);
  H5Eset_auto(NULL, NULL);

  root_id = H5Gopen(this->file_id, "/");

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
  status = H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->Time);
  status = H5Aclose(attr1);
    //cerr << "AMR::RequestInformation for time " << this->Time << endl;

  attr1 = H5Aopen_name(root_id, "Time Scaling Factor");
  status = H5Aread(attr1, H5T_NATIVE_DOUBLE, &TimeScalor);
  status = H5Aclose(attr1);

  attr1 = H5Aopen_name(root_id, "Length Scale Factor");
  if(attr1 >0 )
    {
    status = H5Aread(attr1, H5T_NATIVE_DOUBLE, &this->LengthScaleFactor);
    status = H5Aclose(attr1);
    //cerr << "vtkAMRAmazeReaderInternal::RequestInformation for LengthScaleFactor " << this->LengthScaleFactor << endl;
    }
  else
    this->LengthScaleFactor = 1.0;

      switch(this->ScaleChoice)
        {
        case 0: // pc
          this->LengthScaleFactor  *= 3.08567782e18;
          cerr << __LINE__ << "   Using PARSEC with Length Scale Factor * by 3.08567782e18 = " << this->LengthScaleFactor << "\n";
        break;
        case 1: // AU
          this->LengthScaleFactor  *= 1.49597870700e13;
          cerr << __LINE__ << "   Using AU with Length Scale Factor * by 1.49597870700e13 = " << this->LengthScaleFactor << "\n";
        break;
        case 2: // RSun
         this->LengthScaleFactor  *= 6.96342e10;
         cerr << __LINE__ << "   Using RSun with Length Scale Factor * by 6.96342e10 = " << this->LengthScaleFactor << "\n";
        break;
        case 3:
          //cerr << __LINE__ << "   Length Scale Factor  is untouched\n";
        break;
        }

  if(H5Lexists(root_id, "/Map", H5P_DEFAULT))
    {
    hid_t map_id = H5Gopen(root_id, "/Map");
    attr1 = H5Aopen_name(map_id, "Map Type");
    hid_t t1 = H5Aget_type(attr1);
    char map_type[30];
    status = H5Aread(attr1, t1, map_type);
    H5Tclose(t1); cerr << map_type<<endl;
    if(!strncmp(map_type, "Sphere-LogR", 11))
      {
      this->MappedGrids = Sphere_LogR;
      cerr << "Using Mapped Grids:" << this->MappedGrids << endl;
      }
    else if(!strncmp(map_type, "DCR_Cart2Spheres", 16))
      {
      this->MappedGrids = DCR_Cart2Spheres;
      cerr << "Using Mapped Grids:" << this->MappedGrids << endl;
      }
    status = H5Aclose(attr1);
    status = H5Gclose(map_id);
    }
  
  H5Gclose(root_id);
  H5Eset_auto(func, client_data);
}

void vtkAMRAmazeReaderInternal::SetScaleChoice(int _arg)
{
    if (this->ScaleChoice != _arg)
      {
      this->ScaleChoice = ScaleOption(_arg);
/*
      switch(this->ScaleChoice)
        {
        case 0: // pc
          this->LengthScaleFactor  *= 3.08567782e18;
          cerr << __LINE__ << "   Using PARSEC with Length Scale Factor * by 3.08567782e18 = " << this->LengthScaleFactor << "\n";
        break;
        case 1: // AU
          this->LengthScaleFactor  *= 1.49597870700e13;
          cerr << __LINE__ << "   Using AU with Length Scale Factor * by 1.49597870700e13 = " << this->LengthScaleFactor << "\n";
        break;
        case 2: // RSun
         this->LengthScaleFactor  *= 6.96342e10;
         cerr << __LINE__ << "   Using RSun with Length Scale Factor * by 6.96342e10 = " << this->LengthScaleFactor << "\n";
        break;
        case 3:
          cerr << __LINE__ << "   Length Scale Factor  is untouched\n";
        break;
        }
*/
      this->Modified();
      }
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
    return NULL;
}


vtkDoubleArray* vtkAMRAmazeReaderInternal::ReadVar(int levelId, int block, adG_component &variable)
{
  hid_t level_root_id, grid_root_id, dataset_id, status, mem_space_id;
  char lname[16];
  int domain = this->FindDomainId(levelId, block);
  /*cerr << "domain = " << domain 
       << ", level = " << levelId 
       << ", block = " << block 
       << ", varname = " << variable.label 
       << endl;*/
  adG_grid grid = this->Grids[domain];

  sprintf(lname, "/Level %d", levelId);
  OpenHDF5File("ReadVar");
  level_root_id = H5Gopen(this->file_id, lname);
  if(level_root_id < 0)
    cerr << __LINE__ << ": ReadVar() bad level_root_id returned\n";

  sprintf(lname, "Grid %d", grid.grid_nr);  //lname is reused here
  grid_root_id = H5Gopen(level_root_id, lname);
  if(grid_root_id < 0)
    cerr << __LINE__<< ": ReadVar(): bad grid_root_id returned\n";

// for 2D case with 2D vectors, we create a 3-tuple anyway, fill in the 3rd component with zeroes
// and we must define a hyperslab select to only fill in the first 2 columns.
  vtkDoubleArray*scalars = vtkDoubleArray::New();
  scalars->SetNumberOfComponents(variable.vec_len == 2? 3 : variable.vec_len);
  scalars->SetName((const char*)PVlabels[(const char *)variable.label].c_str());
  //cerr << __LINE__ << ": SetName( " << PVlabels[variable.label] << ")\n";
// default naming. Could be over-written by "Log10()"

  int nvals;
  if(GetCellCentered())
    {
    if(grid.dimensions[2] != 1)
      nvals = (grid.dimensions[0]-1) * (grid.dimensions[1]-1) * (grid.dimensions[2]-1);
    else
      nvals = (grid.dimensions[0]-1) * (grid.dimensions[1]-1);
    }
  else
    nvals = grid.dimensions[0] * grid.dimensions[1] * grid.dimensions[2];
/*
  cerr << lname << ":" << PVlabels[(const char *)variable.label] << "("<< variable.vec_len << "," << nvals << ")\n";
  cerr << "nvals = " << grid.dimensions[0] << "x"<<
                    grid.dimensions[1]<< "x"<<
                    grid.dimensions[2]<< endl;
*/
  hsize_t dims[2], count[2], offset[2];
  dims[0] = nvals;
  dims[1] =  variable.vec_len == 2 ? 3 : variable.vec_len;
  mem_space_id = H5Screate_simple(2,  dims, NULL);

  offset[0] = 0; offset[1] = 0;
  count[0] = nvals; count[1] = variable.vec_len;
  H5Sselect_hyperslab (mem_space_id, H5S_SELECT_SET, offset, NULL, count, NULL);

  scalars->SetNumberOfTuples(nvals);
  //cerr << "reading PointDataArray " << variable.label<<endl;
  void *dataArray = scalars->GetVoidPointer(0);

  dataset_id = H5Dopen(grid_root_id, (const char *) variable.label);
  if(dataset_id < 0)
     {
     cerr << "error opening HDF5 dataset for var " << variable.label << endl;
     }

  if(H5Dread(dataset_id, H5T_NATIVE_DOUBLE, mem_space_id, H5S_ALL, H5P_DEFAULT, dataArray) < 0)
     {
     cerr << "951:error reading HDF5 dataset for var " << variable.label << endl;
     }

  double *dArray = (double *)dataArray;
  if(variable.vec_len == 2) for(int k=0; k < nvals * 3; k+=3) dArray[k+2] = 0;

  if(this->DataScale == true && variable.scalefactor != 1.0)
    {
    vtkDebugMacro( << "should divide by scaling factor "<< variable.scalefactor << " for " << *(variable.label));
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
  //CloseHDF5File(" ReadVar");
  return scalars;
}

int vtkAMRAmazeReaderInternal::CheckVarSize(int levelId, int block, adG_component &variable)
{
  hid_t level_root_id, grid_root_id, dataset_id, status, mem_space_id;
  char lname[16];
  int domain = this->FindDomainId(levelId, block);
  /*cerr << "domain = " << domain 
       << ", level = " << levelId 
       << ", block = " << block 
       << ", varname = " << variable.label 
       << endl;*/
  adG_grid grid = this->Grids[domain];

  sprintf(lname, "/Level %d", levelId);
  level_root_id = H5Gopen(this->file_id, lname);
  if(level_root_id < 0)
    cerr << "543 bad level_root_id returned\n";
  else
    {
    sprintf(lname, "Grid %d", grid.grid_nr);  //lname is reused here
    if(H5Lexists(level_root_id, lname, H5P_DEFAULT))
      {
      grid_root_id = H5Gopen(level_root_id, lname);
      if(grid_root_id < 0)
        cerr << "ReadVar(): bad grid_root_id returned\n";

      int nvals = grid.dimensions[0] * grid.dimensions[1] * grid.dimensions[2];
  /*
  cerr << lname << ":" << PVlabels[(const char *)variable.label] << "("<< variable.vec_len << "," << nvals << ")\n";
  cerr << "nvals = " << grid.dimensions[0] << "x"<<
                    grid.dimensions[1]<< "x"<<
                    grid.dimensions[2]<< endl;
*/
      dataset_id = H5Dopen(grid_root_id, (const char *) variable.label);
      if(dataset_id < 0)
        {
        cerr << "error opening HDF5 dataset for var " << variable.label << endl;
        }
      hid_t space_id;
      hsize_t dims[3], maxdims;
      space_id = H5Dget_space(dataset_id);
      H5Sget_simple_extent_dims(space_id, dims, NULL );

      H5Sclose(space_id);
      H5Dclose(dataset_id);
      H5Gclose(grid_root_id);
      H5Gclose(level_root_id);

      if(nvals  == dims[0])
        {
        this->CellCenteredOff();
        return 0;
        }
      else
        {
        this->CellCenteredOn();
        return 1;
        }
      }
    }
}

// The ParaView reader already has the notion of level and block-id within that level
// so we use that for the method's signature and we calculate the domain-id
// domain-id is between 0 and N-1

vtkUniformGrid* vtkAMRAmazeReaderInternal::ReadUniformGrid(int levelId, int block)
{
  int domain = this->FindDomainId(levelId, block);
  adG_grid grid = this->Grids[domain];
  // uniform is uniform, so no need to get all 3 DXs. One is enough

  double dx[3]; // spacing is constant at a given level
  dx[0] = this->Levels[levelId].DXs[0];
  dx[1] = this->Levels[levelId].DXs[1];
  dx[2] = this->Levels[levelId].DXs[2];

  vtkStringArray *sarr = vtkStringArray::New();
  sarr->SetName("GridName");
  sarr->SetNumberOfComponents(1);
  sarr->SetNumberOfTuples(1);
  char name[16];
  sprintf(name, "Grid %d", grid.grid_nr);
  sarr->SetValue(0, name);

  vtkUniformGrid* ug = vtkUniformGrid::New();
  ug->Initialize();
  ug->GetFieldData()->AddArray(sarr);
  sarr->Delete();
  vtkDoubleArray *data = vtkDoubleArray::New();
  data->SetName("Time");
  data->InsertValue(0, this->Time);
  ug->GetFieldData()->AddArray(data);
  data->Delete();
/*
  if(this->LengthScale)
    cerr << "LengthScaleFactor = " << this->LengthScaleFactor<<"\n";
*/
  if(grid.dimensions[2] > 1)  // a 3D grid
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
  //ug->SetWholeExtent(0, grid.dimensions[0]-1,
                     //0, grid.dimensions[1]-1,
                     //0, grid.dimensions[2]-1);
  if(this->LengthScale)
    ug->SetOrigin(grid.origin[0]/this->LengthScaleFactor,
                  grid.origin[1]/this->LengthScaleFactor,
                  grid.origin[2]/this->LengthScaleFactor);
  else
    ug->SetOrigin(grid.origin[0], grid.origin[1], grid.origin[2]);
/*
cerr << "UG (L=" << levelId << ", b=" << block << ", ): at Origin " << grid.origin[0] << ", "<<
                    grid.origin[1]<< ", "<<
                    grid.origin[2]<<endl;
*/
  ug->SetDimensions(grid.dimensions[0],
                    grid.dimensions[1],
                    grid.dimensions[2]);

/*
cerr << "UG: of size " << grid.dimensions[0] << "x"<<
                    grid.dimensions[1]<< "x"<<
                    grid.dimensions[2]<<endl;
*/
  return ug;
} // ReadUniformGrid

void
map_rtp2xyz(double R, double T, double P,
            double *x, double *y, double *z)
{
  *x = R * cos(T);
  *y = R * sin(T); // intermediate y calculation
  *z = *y * sin(P);
  *y *= cos(P);
//cerr << " " << R << " " << T << " " << P << " => " << *x << " " << *y << " " << *z << endl;
}

vtkStructuredGrid* vtkAMRAmazeReaderInternal::ReadStructuredGrid(int domain)
{
  int I, nvals, levelId, blockId;
  double x, y, z;
  adG_grid grid = this->Grids[domain];
  this->FindLevelAndBlock(domain, levelId, blockId);
  double dx[3]; // spacing is constant at a given level
  dx[0] = this->Levels[levelId].DXs[0];
  dx[1] = this->Levels[levelId].DXs[1];
  dx[2] = this->Levels[levelId].DXs[2];

/*
  cerr << "Phi ranges from "    << grid.origin[2] << " to " << grid.origin[5] << " with increment " << dx[2] << endl;
  cerr << "Theta ranges from "  << grid.origin[1] << " to " << grid.origin[4] << " with increment " << dx[1] << endl;
  cerr << "Radius ranges from " << grid.origin[0] << " to " << grid.origin[3] << " with increment " << dx[0] << endl;

*/
  vtkStringArray *sarr = vtkStringArray::New();
  sarr->SetName("GridName");
  sarr->SetNumberOfComponents(1);
  sarr->SetNumberOfTuples(1);
  char name[16];
  sprintf(name, "Grid %d", grid.grid_nr);
  sarr->SetValue(0, name);

  vtkStructuredGrid* sg = vtkStructuredGrid::New();
  sg->GetFieldData()->AddArray(sarr);
  sarr->Delete();
  vtkDoubleArray *data = vtkDoubleArray::New();
  data->SetName("Time");
  data->InsertValue(0, this->Time);
  sg->GetFieldData()->AddArray(data);
  data->Delete();

  sg->SetDimensions(grid.dimensions[0],
                    grid.dimensions[1],
                    grid.dimensions[2]);
  nvals = grid.dimensions[0] * grid.dimensions[1] * grid.dimensions[2];

  vtkDoubleArray *coords = vtkDoubleArray::New();
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(nvals);

  int Iphi, Itheta, Iradius;

  //dx[0] = dx[1]; //. Rolf said the example is wrong 
  double NLevel_R = (1.0 - 0.0) / dx[0];
  //cerr << "NLevel_R " << NLevel_R << endl;
  double Delta_Level_T = (this->SphereLogRMappings[1].Tmax - this->SphereLogRMappings[1].Tmin)/((1.0 - 0.0) / dx[1]);
  //cerr << "Delta_Level_T " << Delta_Level_T << endl;
  double Delta_Level_P = (this->SphereLogRMappings[1].Pmax - this->SphereLogRMappings[1].Pmin)/((1.0 - 0.0) / dx[2]);
  //cerr << "Delta_Level_P " << Delta_Level_P << endl;

  double alpha = pow((this->SphereLogRMappings[1].Rmax/this->SphereLogRMappings[1].Rmin), 1.0/NLevel_R) - 1;
  //cerr << "alpha " << alpha << endl;
  for(Iphi=0; Iphi < grid.dimensions[2]; Iphi++)
    {
    double arg_phi = this->SphereLogRMappings[1].Pmin + (grid.box_corners[2] + Iphi)*Delta_Level_P;
    for(Itheta=0; Itheta < grid.dimensions[1]; Itheta++)
      {
      double arg_theta = this->SphereLogRMappings[1].Tmin + (grid.box_corners[1] + Itheta)*Delta_Level_T;
// so the fastest index is the radial index. 
// cell 0 is at the lower right corner of the 2D map, and then the theta sweep goes from 0 to PI in counter-clockwise fashion
      I = Iphi*(grid.dimensions[1]*grid.dimensions[0]) + Itheta*grid.dimensions[0];
      for(Iradius=0; Iradius < grid.dimensions[0]; Iradius++)
        {
         //if(Itheta==0) cerr << this->SphereLogRMappings[1].Rmin * pow(1.0 + alpha, grid.box_corners[0]+Iradius) << "\n";
        double R = this->SphereLogRMappings[1].Rmin * pow(1.0 + alpha, grid.box_corners[0]+Iradius);
        R /= this->LengthScaleFactor;
        map_rtp2xyz(R, arg_theta, arg_phi, &x, &y, &z);
        //coords->SetTuple3(I+Iradius, param[1].Rmin * pow(1.0 + alpha, grid.box_corners[0]+Iradius), arg_theta, arg_phi);
        coords->SetTuple3(I+Iradius, x, y, z);
        }//if(Itheta==0) cerr << endl;
      }
    }

  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);
  coords->Delete();
  sg->SetPoints(points);
  points->Delete();

  return sg;
} // ReadStructuredGrid

void
map_c2p_fig31(double xc, double yc, double zc, double R1,
        double *xp, double *yp, double *zp)
{
  double d, r;
  if (fabs(xc) > fabs(yc))
     d = fabs(xc);
  else
     d = fabs(yc);
  r = sqrt(xc*xc + yc*yc);
  if (r < 1.e-10)
     r = 1.e-10;
  *xp = R1 * d * xc/r;
  *yp = R1 * d * yc/r;
  *zp = 0.0;
//cerr << " " << R << " " << T << " " << P << " => " << *x << " " << *y << " " << *z << endl;
}

void
map_c2p_fig32a(double xc, double yc, double zc, double R1,
        double *xp, double *yp, double *zp)
{
//inputs are xc, yc, zc, R1
//outputs are xp, yp, zp
  double R, center, D, d, absxc, absyc;
  absxc = fabs(xc);
  absyc = fabs(yc);

  if (absxc > absyc)
     d = absxc;
  else
     d = absyc;

  if (d < 1.e-10)
     d = 1.e-10;
  D = R1 * d / M_SQRT2;
  R = R1 * d;
  R = R*R; // we do this here because we only use R^2 until the end and exit.
  center = D -sqrt(R - D*D);

  D /= d; // used twice below
  *xp = D * absxc;
  *yp = D * absyc;
  
  
  if(absyc >= absxc)
    *yp = center + sqrt(R - *xp * *xp);
  if(absxc >= absyc)
    *xp = center + sqrt(R - *yp * *yp);
  if(xc < 0)
    *xp = -1. * *xp;
  if(yc < 0)
    *yp = -1. * *yp;
  *zp = 0.0;
}

void
map_c2p_fig32b(double xc, double yc, double zc, double R1,
        double *xp, double *yp, double *zp)
{
//inputs are xc, yc, zc, R1
//outputs are xp, yp, zp
  double R, center, D, d, absxc, absyc;
  absxc = fabs(xc);
  absyc = fabs(yc);

  if (absxc > absyc)
     d = absxc;
  else
     d = absyc;

  if (d < 1.e-10)
     d = 1.e-10;
  D = R1 * d / M_SQRT2;
  R = R1;
  R = R*R; // we do this here because we only use R^2 until the end and exit.
  center = D -sqrt(R - D*D);

  D /= d; // used twice below
  *xp = D * absxc;
  *yp = D * absyc;
  
  
  if(absyc >= absxc)
    *yp = center + sqrt(R - *xp * *xp);
  if(absxc >= absyc)
    *xp = center + sqrt(R - *yp * *yp);
  if(xc < 0)
    *xp = -1. * *xp;
  if(yc < 0)
    *yp = -1. * *yp;
  *zp = 0.0;
}

void
map_c2p_fig32c(double xc, double yc, double zc, double R1,
        double *xp, double *yp, double *zp)
{
//inputs are xc, yc, zc, R1
//outputs are xp, yp, zp
  double R, center, D, d, absxc, absyc;
  absxc = fabs(xc);
  absyc = fabs(yc);

  if (absxc > absyc)
     d = absxc;
  else
     d = absyc;

  if (d < 1.e-10)
     d = 1.e-10;
  D = R1 * d * (2.0 - d)/ M_SQRT2;
  R = R1;
  R = R*R; // we do this here because we only use R^2 until the end and exit.
  center = D -sqrt(R - D*D);

  D /= d; // used twice below
  *xp = D * absxc;
  *yp = D * absyc;
  
  
  if(absyc >= absxc)
    *yp = center + sqrt(R - *xp * *xp);
  if(absxc >= absyc)
    *xp = center + sqrt(R - *yp * *yp);
  if(xc < 0)
    *xp = -1. * *xp;
  if(yc < 0)
    *yp = -1. * *yp;
  *zp = 0.0;
}


void
map_c2p_fig32d(double xc, double yc, double zc, double R1,
        double *xp, double *yp, double *zp)
{
//inputs are xc, yc, zc, R1
//outputs are xp, yp, zp
  double d, r, w;
  if (fabs(xc) > fabs(yc))
     d = fabs(xc);
  else
     d = fabs(yc);
  r = sqrt(xc*xc + yc*yc);
  if (r < 1.e-10)
     r = 1.e-10;
  *xp = R1 * d * xc/r;
  *yp = R1 * d * yc/r;
  w = d*d;
  *xp = w*(*xp) + (1.-w) * R1 * xc / M_SQRT2;
  *yp = w*(*yp) + (1.-w) * R1 * yc / M_SQRT2;
  *zp = 0.0;
//cerr << " " << R << " " << T << " " << P << " => " << *x << " " << *y << " " << *z << endl;
}

vtkStructuredGrid* vtkAMRAmazeReaderInternal::ReadStructuredGrid2(int domain)
{
  int I, nvals, levelId, blockId;
  double x, y, z;
  adG_grid grid = this->Grids[domain];
  this->FindLevelAndBlock(domain, levelId, blockId);
  double dx[3]; // spacing is constant at a given level
  dx[0] = this->Levels[levelId].DXs[0];
  dx[1] = this->Levels[levelId].DXs[1];
  dx[2] = this->Levels[levelId].DXs[2];

  //cerr << "Phi ranges from "    << grid.origin[2] << " to " << grid.origin[5] << " with increment " << dx[2] << endl;
  //cerr << "Theta ranges from "  << grid.origin[1] << " to " << grid.origin[4] << " with increment " << dx[1] << endl;
  //cerr << "Radius ranges from " << grid.origin[0] << " to " << grid.origin[3] << " with increment " << dx[0] << endl;

  vtkStringArray *sarr = vtkStringArray::New();
  sarr->SetName("GridName");
  sarr->SetNumberOfComponents(1);
  sarr->SetNumberOfTuples(1);
  char name[16];
  sprintf(name, "Grid %d", grid.grid_nr);
  sarr->SetValue(0, name);

  vtkStructuredGrid* sg = vtkStructuredGrid::New();
  sg->GetFieldData()->AddArray(sarr);
  sarr->Delete();
  vtkDoubleArray *data = vtkDoubleArray::New();
  data->SetName("Time");
  data->InsertValue(0, this->Time);
  sg->GetFieldData()->AddArray(data);
  data->Delete();

  sg->SetDimensions(grid.dimensions[0],
                    grid.dimensions[1],
                    grid.dimensions[2]);
  nvals = grid.dimensions[0] * grid.dimensions[1] * grid.dimensions[2];

  vtkDoubleArray *coords = vtkDoubleArray::New();
  coords->SetNumberOfComponents(3);
  coords->SetNumberOfTuples(nvals);

  int Ix, Iy, Iz;
  double Ox, Oy, Oz, xc, yc, zc;

  if(this->LengthScale)
    {
    Ox = grid.origin[0]/this->LengthScaleFactor;
    Oy = grid.origin[1]/this->LengthScaleFactor;
    Oz = grid.origin[2]/this->LengthScaleFactor;
    }
  else
    {
    Ox = grid.origin[0];
    Oy = grid.origin[1];
    Oz = grid.origin[2];
    }
  double R = this->DCR_Mappings.Rmax;
  //R /= this->LengthScaleFactor;
  if(!strncmp(this->DCR_Mappings.MapCase, "CASE_A", 6))
    {
    cerr << "CASE_A"    << endl;
    for(Iz=0; Iz < grid.dimensions[2]; Iz++)
    {
    for(Iy=0; Iy < grid.dimensions[1]; Iy++)
      {
      I = Iz*(grid.dimensions[1]*grid.dimensions[0]) + Iy*grid.dimensions[0];
      for(Ix=0; Ix < grid.dimensions[0]; Ix++)
        {
        if(!strncmp(this->DCR_Mappings.MapLunarity, "FULL", 4))
          {
          xc = Ox + (2*Ix-(grid.dimensions[0]-1))*dx[0];
          yc = Oy + (2*Iy-(grid.dimensions[1]-1))*dx[1];
          zc = 0;
          }
        else if(!strncmp(this->DCR_Mappings.MapLunarity, "HALF", 4))
          {
          xc = Ox + (2*Ix-(grid.dimensions[0]-1))*dx[0];
          yc = Oy + Iy*dx[1];
          zc = 0;
          }
        else if(!strncmp(this->DCR_Mappings.MapLunarity, "QUARTER", 7))
          {
          xc = Ox + Ix*dx[0];
          yc = Oy + Iy*dx[1];
          zc = 0;
          }
        map_c2p_fig32a(xc, yc, zc, R, &x, &y, &z);
        coords->SetTuple3(I+Ix, x, y, z);
        }
      }
    }
    }
  if(!strncmp(this->DCR_Mappings.MapCase, "CASE_B", 6))
    {
    cerr << "CASE_B"    << endl;
    for(Iz=0; Iz < grid.dimensions[2]; Iz++)
    {
    for(Iy=0; Iy < grid.dimensions[1]; Iy++)
      {
      I = Iz*(grid.dimensions[1]*grid.dimensions[0]) + Iy*grid.dimensions[0];
      for(Ix=0; Ix < grid.dimensions[0]; Ix++)
        {
        if(!strncmp(this->DCR_Mappings.MapLunarity, "FULL", 4))
          {
          xc = Ox + (2*Ix-(grid.dimensions[0]-1))*dx[0];
          yc = Oy + (2*Iy-(grid.dimensions[1]-1))*dx[1];
          zc = 0;
          }
        else if(!strncmp(this->DCR_Mappings.MapLunarity, "HALF", 4))
          {
          xc = Ox + (2*Ix-(grid.dimensions[0]-1))*dx[0];
          yc = Oy + Iy*dx[1];
          zc = 0;
          }
        else if(!strncmp(this->DCR_Mappings.MapLunarity, "QUARTER", 7))
          {
          xc = Ox + Ix*dx[0];
          yc = Oy + Iy*dx[1];
          zc = 0;
          }
        map_c2p_fig32b(xc, yc, zc, R, &x, &y, &z);
        coords->SetTuple3(I+Ix, x, y, z);
        }
      }
    }
    }
  if(!strncmp(this->DCR_Mappings.MapCase, "CASE_C", 6))
    {
    cerr << "CASE_C"    << endl;
    for(Iz=0; Iz < grid.dimensions[2]; Iz++)
    {
    for(Iy=0; Iy < grid.dimensions[1]; Iy++)
      {
      I = Iz*(grid.dimensions[1]*grid.dimensions[0]) + Iy*grid.dimensions[0];
      for(Ix=0; Ix < grid.dimensions[0]; Ix++)
        {
        if(!strncmp(this->DCR_Mappings.MapLunarity, "FULL", 4))
          {
          xc = Ox + (2*Ix-(grid.dimensions[0]-1))*dx[0];
          yc = Oy + (2*Iy-(grid.dimensions[1]-1))*dx[1];
          zc = 0;
          }
        else if(!strncmp(this->DCR_Mappings.MapLunarity, "HALF", 4))
          {
          xc = Ox + (2*Ix-(grid.dimensions[0]-1))*dx[0];
          yc = Oy + Iy*dx[1];
          zc = 0;
          }
        else if(!strncmp(this->DCR_Mappings.MapLunarity, "QUARTER", 7))
          {
          xc = Ox + Ix*dx[0];
          yc = Oy + Iy*dx[1];
          zc = 0;
          }
        map_c2p_fig32c(xc, yc, zc, R, &x, &y, &z);
        coords->SetTuple3(I+Ix, x, y, z);
        }
      }
    }
    }
  else if(!strncmp(this->DCR_Mappings.MapCase, "CASE_D", 6))
    {
    cerr << "CASE_D"    << endl;
    for(Iz=0; Iz < grid.dimensions[2]; Iz++)
    {
    for(Iy=0; Iy < grid.dimensions[1]; Iy++)
      {
      I = Iz*(grid.dimensions[1]*grid.dimensions[0]) + Iy*grid.dimensions[0];
      for(Ix=0; Ix < grid.dimensions[0]; Ix++)
        {
        if(!strncmp(this->DCR_Mappings.MapLunarity, "FULL", 4))
          {
          xc = Ox + (2*Ix-(grid.dimensions[0]-1))*dx[0];
          yc = Oy + (2*Iy-(grid.dimensions[1]-1))*dx[1];
          zc = 0;
          }
        else if(!strncmp(this->DCR_Mappings.MapLunarity, "HALF", 4))
          {
          xc = Ox + (2*Ix-(grid.dimensions[0]-1))*dx[0];
          yc = Oy + Iy*dx[1];
          zc = 0;
          }
        else if(!strncmp(this->DCR_Mappings.MapLunarity, "QUARTER", 7))
          {
          xc = Ox + Ix*dx[0];
          yc = Oy + Iy*dx[1];
          zc = 0;
          }
        map_c2p_fig32d(xc, yc, zc, R, &x, &y, &z);
        //coords->SetTuple3(I+Iradius, param[1].Rmin * pow(1.0 + alpha, grid.box_corners[0]+Iradius), arg_theta, arg_phi);
        coords->SetTuple3(I+Ix, x, y, z);
        }//if(Itheta==0) cerr << endl;
      }
    }
    }

  vtkPoints *points = vtkPoints::New();
  points->SetData(coords);
  coords->Delete();
  sg->SetPoints(points);
  points->Delete();

  return sg;
} // ReadStructuredGrid2

vtkRectilinearGrid* vtkAMRAmazeReaderInternal::ReadRectilinearGrid(int domain)
{
  int i, levelId, blockId;

  adG_grid grid = this->Grids[domain];
  this->FindLevelAndBlock(domain, levelId, blockId);

  double dx[3]; // spacing is constant at a given level
  dx[0] = this->Levels[levelId].DXs[0];
  dx[1] = this->Levels[levelId].DXs[1];
  dx[2] = this->Levels[levelId].DXs[2];
  vtkRectilinearGrid* rg = vtkRectilinearGrid::New();
/*
  vtkCharArray *nameArray = vtkCharArray::New();
  nameArray->SetName("Name");
  char *name = nameArray->WritePointer(0, 20);
  sprintf(name, "Grid %d", grid.grid_nr);

  rg->GetFieldData()->AddArray(nameArray);
  nameArray->Delete();
  vtkDoubleArray *data = vtkDoubleArray::New();
  data->SetName("Time");
  data->InsertValue(0, this->Time);
  rg->GetFieldData()->AddArray(data);
  data->Delete();
*/
  vtkDoubleArray *xcoords = vtkDoubleArray::New();
  vtkDoubleArray *ycoords = vtkDoubleArray::New();
  vtkDoubleArray *zcoords = vtkDoubleArray::New();

  xcoords->SetNumberOfComponents(1);
  ycoords->SetNumberOfComponents(1);
  zcoords->SetNumberOfComponents(1);

  xcoords->SetNumberOfTuples(grid.dimensions[0]);
  ycoords->SetNumberOfTuples(grid.dimensions[1]);

  if(this->LengthScale)
    {
  //cerr << "LengthScaleFactor = " << this->LengthScaleFactor<<"\n";
    dx[0] = dx[0]/this->LengthScaleFactor;
    dx[1] = dx[1]/this->LengthScaleFactor;
    dx[2] = dx[2]/this->LengthScaleFactor;
    }
  //cerr << "RectGrid(dx0,dx1,dx2) "<< dx[0] << ", "<< dx[1] << "," << dx[2]<<"\n";
  if(grid.dimensions[2] > 1)  // a 3D grid
    {
    zcoords->SetNumberOfTuples(grid.dimensions[2]);
    }
  else
    {
    zcoords->SetNumberOfTuples(1);
    }

  //rg->SetWholeExtent(0, grid.dimensions[0]-1,
                     //0, grid.dimensions[1]-1,
                     //0, grid.dimensions[2]-1);
  double origin;
  if(this->LengthScale)
    origin = grid.origin[0]/this->LengthScaleFactor;
  else
    origin = grid.origin[0];

  for(i=0; i < grid.dimensions[0]; i++)
    {
    xcoords->SetValue(i, (origin + i*dx[0]) );
    }
  //cerr << "Grid(x0,x1, y0,y1,z0,z1) "<< domain << ": "<<origin << "," << (origin + (i-1)*dx[0])<<", ";
  rg->SetXCoordinates(xcoords);
  xcoords->Delete();
/////////////////////////////////////////
  if(this->LengthScale)
    origin = grid.origin[1]/this->LengthScaleFactor;
  else
    origin = grid.origin[1];

  for(i=0; i < grid.dimensions[1]; i++)
    {
    ycoords->SetValue(i, origin + i*dx[1] );
    }
  //cerr << origin << "," << origin + (i-1)*dx[1]<<", ";
  rg->SetYCoordinates(ycoords);
  ycoords->Delete();
/////////////////////////////////////////
  if(this->LengthScale)
    origin = grid.origin[2]/this->LengthScaleFactor;
  else
    origin = grid.origin[2];

  for(i=0; i < grid.dimensions[2]; i++)
    {
    zcoords->SetValue(i, origin + i*dx[2] );
    }
  //cerr << origin << "," << origin + (i-1)*dx[2]<<"\n";
  rg->SetZCoordinates(zcoords);
  zcoords->Delete();
/////////////////////////////////////////

  rg->SetDimensions(grid.dimensions[0],
                    grid.dimensions[1],
                    grid.dimensions[2]);

/*
cerr << "RG: of size " << grid.dimensions[0] << "x"<<
                    grid.dimensions[1]<< "x"<<
                    grid.dimensions[2]<<endl;
*/

  return rg;
} // ReadRectilinearGrid()

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
/*
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
*/

// returns the global id for the given block
int vtkAMRAmazeReaderInternal::FindDomainId(int level, int block)
{
  int domain = 0;
  for(int l=0; l < level; l++)
    domain += this->Levels[l].GridsPerLevel;
  domain += block;
  if (domain >= this->NumberOfGrids)
    cerr << "level or block too high for this dataset\n";
  return domain;
}

void vtkAMRAmazeReaderInternal::GetSpacing(int level, double *spacing)
{
  int domain = this->FindDomainId(level, 0);
  adG_grid g = this->Grids[domain];

    if(GetDimensionality() == 3)  // a 3D grid
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

void vtkAMRAmazeReaderInternal::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  //os << "FileName: " << (this->FileName? this->FileName:"(none)") << "\n";
  os << "FileName: " << this->FileName  << "\n";
  os << indent << "Dimensionality: " << this->Dimensionality << endl;
  os << indent << "NumberOfLevels: " << this->GetNumberOfLevels() << endl;

  int i, nb_of_grids=0;
  for(i=0; i < this->NumberOfLevels; i++)
    {
    nb_of_grids += this->Levels[i].GridsPerLevel;
    }
  if(nb_of_grids != this->NumberOfGrids)
    cerr << "mismatch between number of grids and levels\n";

  for(i=0; i < this->NumberOfGrids; i++)
    {
    int block, level, domain;
    this->FindLevelAndBlock(i, level, block);
    //cout << setw(4) << i << " : level " << setw(2) << level << ", block " << setw(3) << block << endl;
    domain = this->FindDomainId(level, block);
    //cout << setw(4) << domain << " : level " << setw(2) << level << ", block " << setw(3) << block << endl;
    }
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
  H5E_auto_t func;
  void *client_data;
  H5Eget_auto(&func, &client_data);
  H5Eset_auto(NULL, NULL);
  OpenHDF5File("BuildStars");
  apr_root_id = H5Gopen(this->file_id, "/APR_StellarSystems");

  if(apr_root_id < 0)
    {
    cerr << __LINE__ << "/APR_StellarSystems not found in HDF5 file\n";
    return 0;
    }
  interactions_root_id = H5Gopen(apr_root_id, "Interactions");
  if(interactions_root_id >=0)
    {
    attr1 = Create_Interaction_Compound();
    interactions = new interaction[2];

    dataset1 = H5Dopen(interactions_root_id, "IsotrInfWind");
    if(dataset1 >=0 )
      {
      status = H5Dread(dataset1, attr1, H5S_ALL, H5S_ALL, H5P_DEFAULT, &interactions[0]);
      }
    else
      {
      cerr << "error reading Interactions::IsotrInfWind\n";
      }
    dataset2 = H5Dopen(interactions_root_id, "Simple Accretors");
    if(dataset2 >=0 )
      {
      status = H5Dread(dataset2, attr1, H5S_ALL, H5S_ALL, H5P_DEFAULT, &interactions[1]);
      }
    else
      {
      cerr << "error reading Interactions::Simple Accretors\n";
      }
    H5Dclose(dataset1);
    H5Dclose(dataset2);
    H5Gclose(interactions_root_id);
    }

  models_root_id = H5Gopen(apr_root_id, "Stars");
  if(models_root_id < 0)
    {
    cerr << __LINE__ << "/APR_StellarSystems/Stars not found in HDF5 file\n";
    return 0;
    }
  hid_t dataspace;
  int ModelsArraySize;
  hsize_t dims_out[1]; // we assume rank = 1
  dataset1 = H5Dopen(models_root_id, "Stars: Models&Interaction");
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
  dataset2 = H5Dopen(models_root_id, "Stars: Present State");
  if(dataset2 >=0) // old-style stars before november 6, 2008
    {
    cerr << "vtkAMRAmazeReaderInternal::BuildStars(Old-style STARS)\n";
    dataspace = H5Dget_space(dataset2);
    status = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    this->NumberOfSphericallySymmetricStars = nb_stars = dims_out[0];

    stars = new star[nb_stars];
    H5Sclose(dataspace);

    attr2 = Create_Star_Compound();
    status = H5Dread(dataset2, attr2, H5S_ALL, H5S_ALL, H5P_DEFAULT, stars);

    H5Tclose(attr2);
    H5Dclose(dataset2);
    cerr << __LINE__ << " :Found " << nb_stars << " stars\n";

    for (i=0; i < nb_stars; i++)
      {
      vtkDoubleArray *velo = vtkDoubleArray::New();
      velo->SetNumberOfComponents(3);
      velo->SetNumberOfTuples(THETARES*(PHIRES-2)+2);
      velo->SetName("Velocity");
      vtkDoubleArray *mass = vtkDoubleArray::New();
      mass->SetNumberOfComponents(1);
      mass->SetNumberOfTuples(THETARES*(PHIRES-2)+2);
      mass->SetName("Mass");
      vtkDoubleArray *temp = vtkDoubleArray::New();
      temp->SetNumberOfComponents(1);
      temp->SetNumberOfTuples(THETARES*(PHIRES-2)+2);
      temp->SetName((const char*)PVlabels["Temperature"].c_str());
      vtkDoubleArray *lum = vtkDoubleArray::New();
      lum->SetNumberOfComponents(1);
      lum->SetNumberOfTuples(THETARES*(PHIRES-2)+2);
      lum->SetName("Luminosity");

      vtkSphereSource *ss = vtkSphereSource::New();
      ss->DebugOff();
      ss->SetCenter(stars[i].Position[0]/this->LengthScaleFactor, stars[i].Position[1]/this->LengthScaleFactor, stars[i].Position[2]/this->LengthScaleFactor);
      ss->SetThetaResolution(THETARES);
      ss->SetPhiResolution(PHIRES);
      stars[i].Radius *= 1; //6.96e10;

      vtkStringArray *name = vtkStringArray::New();
      name->SetName("Wind or Accretor");
      name->SetNumberOfComponents(1);
      name->SetNumberOfTuples(1);

      if(strstr(stars[i].Interaction, "Wind")) // look for substring "Wind"
        {
        for(j=0; j < ModelsArraySize; j++)
          {
          if(strstr(models[j].IntActType, "Wind"))
            {
            cerr << "found the model: " << models[j].IntActType;
            cerr << "multiply by " << interactions[j].CompRadius << endl;
            stars[i].Radius *= interactions[j].CompRadius;
            name->SetValue(0, "Wind");
            }
          }
        }
      else
        {
        for(j=0; j < ModelsArraySize; j++)
          {
          if(strstr(models[j].IntActType, "Accretor"))
            {
            cerr << "found the model: " << models[j].IntActType;
            cerr << "multiply by " << interactions[j].CompRadius << endl;
            stars[i].Radius *= interactions[j].CompRadius;
            name->SetValue(0, "Accretor");
            }
          }
        }
      ss->SetRadius(stars[i].Radius);
      cerr << "radius = "<< stars[i].Radius << endl;
      ss->Update();
      for(vtkIdType k=0; k < THETARES*(PHIRES-2)+2; k++)
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
      name->Delete(); velo->Delete();mass->Delete();temp->Delete();lum->Delete();
      this->Stars.push_back(ss->GetOutput());
      //ss->Delete();
      }
/*
for(int i=0; i < this->Stars.size(); i++)
  {
  vtkDataSetWriter *writer = vtkDataSetWriter::New();
  writer->SetInput(this->Stars[i]);
  char f[65];
  sprintf(f,"/tmp/tmp/foo.%02d.vtk", i);
  cerr << "writing " << this->Stars[i]->GetNumberOfPoints() << " to disk\n";
  writer->SetFileName((const char*)f);
  writer->Write();
  writer->Delete();
  }*/
    } // end of old-style stars before november 6, 2008

  else if((StarsDS = H5Dopen(models_root_id, "Stars")) >= 0)
    {
    cerr << "\n\n";
    dataspace = H5Dget_space(StarsDS);
    status = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
    nb_stars = dims_out[0];

    newstars = new newstar[nb_stars];
    H5Sclose(dataspace);

    attr2 = Create_NewStar_Compound();
    status = H5Dread(StarsDS, attr2, H5S_ALL, H5S_ALL, H5P_DEFAULT, newstars);

    for(i=0; i < nb_stars; i++)
      {
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

      if(strstr(newstars[i].StarModel, "SpherSymStar"))
        NumberOfSphericallySymmetricStars++;
      else
        NumberOfAxisSymmetricStars++;
      }

    cerr << __LINE__ << " :found " << this->NumberOfSphericallySymmetricStars << " spheric-symmetric stars, and " << this->NumberOfAxisSymmetricStars << " axis-symmetric stars\n";
    H5Tclose(attr2);

    if(this->NumberOfSphericallySymmetricStars)
      {
      hid_t SpherSymStarCurrent_id = Create_SpherSymStar_Compound();
      SpherSymStarCurrent *spherStarData = new SpherSymStarCurrent;
      for (i=0; i < nb_stars; i++)
        {
        double Radius;
        char starname[16];
        sprintf(starname, "Star %d", i+1);
        hid_t attr1, dataspace, starN_G, starN_DS;
        hsize_t dims_out[1]; // we assume rank = 1

        starN_G = H5Gopen(models_root_id, starname);
        if(starN_G >=0)
          {
          starN_DS = H5Dopen(starN_G, "SpherSymStar_Current");
          if(starN_DS >=0)
            {
            status = H5Dread(starN_DS, SpherSymStarCurrent_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, spherStarData);
            H5Dclose(starN_DS);

           cerr
           << spherStarData->Radius  << "\t"
           << spherStarData->Temperature  << "\t"
           << spherStarData->Luminosity  << "\t"
           << spherStarData->Omega[0]  << ", " << spherStarData->Omega[1]  << ", " << spherStarData->Omega[2] << "\t"
           << spherStarData->BField[0]  << ", " << spherStarData->BField[1]  << ", " << spherStarData->BField[2]<< endl;

            Radius = spherStarData->Radius;
            }
          else
            Radius=1.0;
          H5Gclose(starN_G);
          }
        else
          {
          cerr << "Error opening " << starname << endl;
          }

        vtkDoubleArray *velo = vtkDoubleArray::New();
        velo->SetNumberOfComponents(3);
        velo->SetNumberOfTuples(1);
        velo->SetName("Velocity");

        vtkDoubleArray *mass = vtkDoubleArray::New();
        mass->SetNumberOfComponents(1);
        mass->SetNumberOfTuples(1);
        mass->SetName("Mass");

        vtkDoubleArray *radius = vtkDoubleArray::New();
        radius->SetNumberOfComponents(1);
        radius->SetNumberOfTuples(1);
        radius->SetName("Radius");

        velo->SetTypedTuple(0, newstars[i].Velocity);
        mass->SetValue(0, newstars[i].Mass);


        char *p = strchr(newstars[i].InteractionModel, ' ');
        *p = '\0';
        vtkStringArray *name = vtkStringArray::New();
        name->SetNumberOfComponents(1);
        name->SetNumberOfTuples(1);
        name->SetName("InteractionModel");
        name->SetValue(0, newstars[i].InteractionModel);

        //vtkSphereSource *ss = vtkSphereSource::New();
        vtkPolyData *ss = vtkPolyData::New();
        vtkPoints   *starpts = vtkPoints::New();
        ss->SetPoints(starpts);
        starpts->Delete();
        starpts->SetNumberOfPoints(1);

        if(strstr(newstars[i].StarModel, "SpherSymStar"))
          {
          if(this->LengthScale)
            {
            starpts->SetPoint(0, newstars[i].Position[0]/this->LengthScaleFactor,
                          newstars[i].Position[1]/this->LengthScaleFactor,
                          newstars[i].Position[2]/this->LengthScaleFactor);
            }
          else
            {
            starpts->SetPoint(0, newstars[i].Position[0],
                          newstars[i].Position[1],
                          newstars[i].Position[2]);
            }

          radius->SetValue(0, Radius * newstars[i].CompRadiusFrac * 6.96e10 /this->LengthScaleFactor);
          }
        vtkIdType pts[1]={0};
        ss->Allocate(1, 1);
        ss->InsertNextCell(VTK_VERTEX, 1, pts);
        ss->GetFieldData()->AddArray(name);
        ss->GetPointData()->AddArray(velo);
        ss->GetPointData()->AddArray(mass);
        ss->GetPointData()->AddArray(radius);
        name->Delete();
        mass->Delete();
        velo->Delete();
        radius->Delete();
        this->Stars.push_back(ss);
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
          char starname[16];
          sprintf(starname, "Star %d", J+1);

          hid_t attr1, dataspace, starN_G, starN_DS;
          hsize_t dims_out[1]; // we assume rank = 1
          //cerr << "opening " << starname<<"\n";
          starN_G = H5Gopen(models_root_id, starname);
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

            starN_DS = H5Dopen(starN_G, "AxiSymStar_Current");
            dataspace = H5Dget_space(starN_DS);
            status = H5Sget_simple_extent_dims(dataspace, dims_out, NULL);
            if(AngleResolution != dims_out[0])
              cerr << "sanity check: NTheta resolution is mis-read?\n";

            AngleResolution = dims_out[0];
            //cout << "Found phi array of size " << AngleResolution << " for " << starname << endl;
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
            vtkTransformPolyDataFilter *tfpd = vtkTransformPolyDataFilter::New();
            tfpd->SetTransform(tf);
            tfpd->SetInputData(AxiSymStar);
            tfpd->Update();
            tf->Delete();
            AxiSymStar->Delete();
            this->Stars.push_back(tfpd->GetOutput());
            }
          else if(AxisDirection == 2)
            {
            vtkTransform *tf = vtkTransform::New();
            tf->Translate(newstars[I].Position[0], 0.0, newstars[I].Position[0]);
            tf->Scale(1.0, 1.0, 1.0);
            tf->RotateY(90.0);
            vtkTransformPolyDataFilter *tfpd = vtkTransformPolyDataFilter::New();
            tfpd->SetTransform(tf);
            tfpd->SetInputData(AxiSymStar);
            tfpd->Update();
            tf->Delete();
            AxiSymStar->Delete();
            this->Stars.push_back(tfpd->GetOutput());
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
  //CloseHDF5File("BuildStars");
  H5Eset_auto(func, client_data);
  //cerr << __LINE__ << " :exit BuildStars() with " << this->NumberOfSphericallySymmetricStars << " spheric-symmetric stars, and " << this->NumberOfAxisSymmetricStars << " axis-symmetric stars\n";

  return nb_stars;
} // ::BuildStars

/*
for(int i=0; i < this->Stars.size(); i++)
  {
  vtkDataSetWriter *writer = vtkDataSetWriter::New();
  writer->SetInputData(this->Stars[i]);
  char f[65];
  sprintf(f,"/tmp/tmp/foo.%02d.vtk", i);
  cerr << "writing " << this->Stars[i]->GetNumberOfPoints() << " to disk\n";
  writer->SetFileName((const char*)f);
  writer->Write();
  writer->Delete();
} 
*/
