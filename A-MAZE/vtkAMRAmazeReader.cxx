#include "vtkAMRAmazeReader.h"
#include "vtkObjectFactory.h"
#include "vtkOverlappingAMR.h"
#include "vtkUniformGrid.h"
#include "vtkDataArraySelection.h"
#include "vtkDataArray.h"
#include "vtkPolyData.h"

#include "vtkByteSwap.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkCellArray.h"

#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkErrorCode.h"
#include "vtkInformation.h"
#include "vtkAMRInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"

#include "vtkPointData.h"
#include "vtkPoints.h"

#include "vtkMultiBlockDataSet.h"
#include "vtkStringArray.h"

#include "vtkTimerLog.h"
#include "vtkAMRUtilities.h"
#define PARALLEL_DEBUG 1

#include <vector>
#include <string>
#include <map>
#include <sstream>
//vtkCxxRevisionMacro(vtkAMRAmazeReader, "$Revision: 1.2 $");
vtkStandardNewMacro(vtkAMRAmazeReader);

//
// http://www.paraview.org/ParaView/index.php/Multi-Resolution_Rendering_with_Overlapping_AMR
//
vtkAMRAmazeReader::vtkAMRAmazeReader()
{
  this->DataScaleOn();
  this->ShiftedGridOff();
  this->DebugOff();
  this->MaxLevelWrite = -1;

  this->LevelRead[0] = -1;
  this->LevelRead[1] = -1;

  this->LevelRange[0] = -1;
  this->LevelRange[1] = -1;
  this->SetNumberOfInputPorts(0);

  this->SetNumberOfOutputPorts(2);
// this is for port number 1 which we do in all cases.
  vtkMultiBlockDataSet *pd = vtkMultiBlockDataSet::New();
  pd->ReleaseData();
  this->GetExecutive()->SetOutputData(1, pd);
  pd->Delete();
// this is for port number 1 which we do in all cases.
  this->myreader = vtkAMRAmazeReaderInternal::New();
  this->LogDataOn();
  this->Initialize();
}

vtkAMRAmazeReader::~vtkAMRAmazeReader()
{
  if(this->FileName != NULL)
    {
    delete [] this->FileName;
    this->FileName = NULL;
    }

  this->myreader->Delete();
  this->myreader = NULL;
}

void vtkAMRAmazeReader::SetLogData(int _arg)
{
  if (this->GetLogData() != _arg)
    {
    this->LogData = _arg;
    this->myreader->SetLogData(_arg);
    this->Modified();
    }
}

void vtkAMRAmazeReader::SetScaleChoice(int _arg)
{
  if (this->GetScaleChoice() != _arg)
    {
    this->ScaleChoice = _arg;
    this->myreader->SetScaleChoice(_arg);
    //cerr << __LINE__ << ": SetScaleChoice(" << _arg << ")\n";
    this->Modified();
    }
}

void vtkAMRAmazeReader::SetLengthScale(int _arg)
{
  if (this->GetLengthScale() != _arg)
    {
    this->LengthScale = _arg;
    this->myreader->SetLengthScale(_arg);
    this->Modified();
    }
}

int vtkAMRAmazeReader::CanReadFile(const char* fname )
{
  int ret = 0;
  if (! fname )
    return ret;
  //cerr << __LINE__ << ": CanReadFile(" << fname << ")\n";
  hid_t f_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t root_id = H5Gopen(f_id, "/");
  if(H5Lexists(root_id, "/Grid Info", H5P_DEFAULT))
    ret = 1;

  H5Gclose(root_id);
  H5Fclose(f_id);
  return ret;
}

void vtkAMRAmazeReader::SetFileName( const char* fileName )
{
  assert( "pre: Internal AMAZE Reader is NULL" && (this->myreader != NULL) );

  if( fileName && strcmp(fileName,"") &&
     ( ( this->FileName == NULL ) || strcmp( fileName, this->FileName ) ) )
    {
    if( this->FileName )
      {
      delete [] this->FileName;
      this->FileName = NULL;
      this->myreader->SetFileName( NULL );
      }

    this->FileName = new char[ strlen(fileName)+1 ];
    strcpy( this->FileName, fileName );
    this->FileName[ strlen( fileName ) ] = '\0';

    //this->IsReady = true;
    this->myreader->SetFileName( this->FileName );
    this->LoadedMetaData = false;

    this->SetUpDataArraySelections();
    this->InitializeArraySelections();
    }

  this->Modified();
}

int vtkAMRAmazeReader::RequestInformation(
  vtkInformation* request, 
  vtkInformationVector** inputVector, 
  vtkInformationVector* outputVector)
{
  //cerr << __LINE__ << " vtkAMRAmazeReader::RequestInformation\n";
  if (!this->Superclass::RequestInformation(request, inputVector, outputVector))
    {
    return 0;
    }
  this->myreader->MakeVariableNames();
  return 1;
}

int vtkAMRAmazeReader::RequestData(
      vtkInformation* request,
      vtkInformationVector** inputVector,
      vtkInformationVector* outputVector)
{
  //cerr << __LINE__ << " vtkAMRAmazeReader::RequestData\n";
  if (!this->Superclass::RequestData(request, inputVector, outputVector))
    {
    return 0;
    }
  vtkInformation* info = outputVector->GetInformationObject(0);
  int piece = info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  if(piece == 0)// should only proc 0 read the stars?
    {
    //cerr << "vtkAMRAmazeReader::RequestData() Load stars\n\n";
    info = outputVector->GetInformationObject(1);
    vtkMultiBlockDataSet* output2 = vtkMultiBlockDataSet::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
    if (!output2)
      {
      cerr << "did not get a correct output2 object\n";
      }
    else
      {
      //cerr << __LINE__ << " vtkAMRAmazeReader::LoadStars\n";
      int nb_stars = vtkAMRAmazeReader::LoadStars(output2);
      }
    }
  return 1;
}

int vtkAMRAmazeReader::LoadStars(vtkMultiBlockDataSet* SpherSymStars)
{
  hid_t    dataset1, dataset2, StarsDS;
  hid_t    attr1, attr2;
  herr_t   status;

  int  i, j, nb_stars, nb_SpherSymStars=0, nb_AxiSymStars=0;
  herr_t error;
  H5E_auto_t func;
  void *client_data;
  H5Eget_auto(&func, &client_data);
  H5Eset_auto(NULL, NULL);

  this->myreader->BuildStars();

  for(int i=0; i < this->myreader->NumberOfSphericallySymmetricStars+this->myreader->NumberOfAxisSymmetricStars; i++)
    {
    SpherSymStars->SetBlock(i, this->myreader->Stars[i]);
    vtkStringArray *mna = vtkStringArray::SafeDownCast(this->myreader->Stars[i]->GetFieldData()->GetAbstractArray("InteractionModel"));

    SpherSymStars->GetMetaData((unsigned int)i)->Set(vtkCompositeDataSet::NAME(), mna->GetValue(0));
    }

  H5Eset_auto(func, client_data);
  return nb_stars;
}

//----------------------------------------------------------------------------
void vtkAMRAmazeReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << "FileName: " << (this->FileName? this->FileName:"(none)") << "\n";
  os << indent << "Dimensionality: " << this->myreader->GetDimensionality() << endl;
  os << indent << "NumberOfLevels: " << this->GetNumberOfLevels() << endl;
  os << indent << "PointData Available: " << this->GetNumberOfComponents() << endl;
  for(int c=0; c < this->GetNumberOfComponents(); c++)
    {
    if(this->PointDataArraySelection->ArrayIsEnabled(this->GetPointArrayName(c)))
      {
      os << indent << indent << "Read \"" << this->GetPointArrayName(c) << "\"\n";
      }
    else
      {
      os << indent << indent << "Skipped \"" << this->GetPointArrayName(c) << "\"\n";
      }
   }
  os << indent << "Output Ports: " << this->GetNumberOfOutputPorts() << endl;
}


int vtkAMRAmazeReader::FillOutputPortInformation(int port, vtkInformation* info)
{
  //cerr << __LINE__ << " vtkAMRAmazeReader::FillOutputPortInformation port = "  << port  << endl;
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkHierarchicalBoxDataSet");
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  return 1;
}

// Get the second output which contains the stars
//vtkPolyData* vtkAMRAmazeReader::GetStarsOutput()
vtkMultiBlockDataSet* vtkAMRAmazeReader::GetStarsOutput()
{
  if (this->GetNumberOfOutputPorts() < 3)
    {
    return NULL;
    }

  //return vtkPolyData::SafeDownCast(this->GetExecutive()->GetOutputData(1));
  return vtkMultiBlockDataSet::SafeDownCast(this->GetExecutive()->GetOutputData(2));
}


int vtkAMRAmazeReader::GetBlockLevel(const int domain)
{
  int gid = domain;
  int l=0;

  while(gid >= this->myreader->Levels[l].GridsPerLevel)
    {
    gid -= this->myreader->Levels[l].GridsPerLevel;
    l++;
    }
  return l;
}

void vtkAMRAmazeReader::SetUpDataArraySelections()
{
  assert( "pre: Internal AMAZE Reader is NULL" && (this->myreader != NULL) );
  this->myreader->ReadMetaData();

  for(int i=0; i < this->GetNumberOfComponents(); i++)
    {
    //cerr << __LINE__ << " PointDataArraySelection->AddArray(" << this->myreader->Labels[i].label << ")\n";
    this->PointDataArraySelection->AddArray((const char *)this->myreader->Labels[i].label);
    int count = strcspn((const char *)this->myreader->Labels[i].unit, " "); // count how many characters in unit different than " "
    this->myreader->Labels[i].unit[count] = '\0';
    }
}

int vtkAMRAmazeReader::FillMetaData( )
{
  this->myreader->ReadMetaData();
  int levelId, GridId;
  double origin[3];
  int *blocksPerLevel = new int[this->GetNumberOfLevels()];
  for (levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
    {
    blocksPerLevel[levelId] = this->GridsPerLevels(levelId);
    }

  //cerr << __LINE__ << " FillMetaData( levels=" << this->GetNumberOfLevels() << ")\n";
  this->Metadata->Initialize(this->GetNumberOfLevels(), blocksPerLevel);
  if(this->myreader->GetDimensionality() == 2)
    this->Metadata->SetGridDescription(VTK_XY_PLANE);
  else
    this->Metadata->SetGridDescription(VTK_XYZ_GRID);
  delete [] blocksPerLevel;

  this->Metadata->SetOrigin(this->myreader->Grids[0].origin);

  std::vector< int > b2level( this->myreader->NumberOfLevels+1, 0 );
  std::vector<adG_grid> &grid = this->myreader->Grids;

  double spacing[3];
  int globalBoxId = 0;
  for (levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
    {
    this->myreader->GetSpacing(levelId, spacing);
    //cerr << __LINE__ << " :Metadata->SetSpacing(" << levelId << ", " << spacing[0]<< ", " << spacing[1]<< ", " << spacing[2] << ")\n";
    this->Metadata->SetSpacing(levelId, spacing);
    for (GridId=0; GridId < this->GridsPerLevels(levelId); GridId++)
      {
      this->Metadata->SetAMRBox(levelId, GridId, this->myreader->Grids[globalBoxId].amrbox);
      this->Metadata->SetAMRBlockSourceIndex(levelId, GridId, globalBoxId);
      globalBoxId++;
      }
    }

  this->Metadata->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),this->GetTime());

  return( 1 );
}

void vtkAMRAmazeReader::ReadMetaData()
{
  this->myreader->ReadMetaData();
}

vtkUniformGrid* vtkAMRAmazeReader::GetAMRGrid( const int blockIdx )
{
  assert( "pre: Internal AMAZE Reader is NULL" && (this->myreader != NULL) );

  int levelId;// = GetBlockLevel(blockIdx);
  int boxId;// = this->myreader->FindDomainId(levelId, blockIdx);
  this->myreader->FindLevelAndBlock(blockIdx, levelId, boxId);
/*
  cerr << __LINE__ << " GetAMRGrid: blockIdx = " << blockIdx
                   << " levelId = " << levelId
                   << " boxId = " << boxId << "\n";*/
  vtkUniformGrid* ug = this->myreader->ReadUniformGrid(levelId, boxId);

  return( ug );
}

void vtkAMRAmazeReader::GetAMRGridPointData(const int blockIdx, vtkUniformGrid *block, const char *field)
{
  //cerr << __LINE__ << " GetAMRGridPointData: blockIdx = " << blockIdx
                   //<< " field = " << field << "\n";
  vtkDataArray *data = this->myreader->ReadVisItVar(blockIdx, field);
  //cerr << __LINE__ << " GetAMRGridPointData: got " << data->GetName() << "\n";
  //data->SetName(field);
  //cerr << __LINE__ << " data->SetName " << field << "\n";
  block->GetPointData()->AddArray(data);
  data->Delete();
}


// ANYTHING BELOW HERE is OLD
#ifdef OLD
//----------------------------------------------------------------------------
int vtkAMRAmazeReader::RequestInformation(
  vtkInformation* request, 
  vtkInformationVector** inputVector, 
  vtkInformationVector* outputVector)
{
  if (!this->Superclass::RequestInformation(request, inputVector, outputVector))
    {
    return 0;
    }

  //cerr << "Begin RequestInformation\n";
//cerr << "AMR:Information\n";
  vtkInformation* info = outputVector->GetInformationObject(0);
  //info->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

  //vtkInformation* info1 = outputVector->GetInformationObject(1);
  //info1->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

  vtkOverlappingAMR *output = static_cast<vtkOverlappingAMR *>(
    info->Get(vtkDataObject::DATA_OBJECT()));

  //FILE *fp=NULL;
  int levelId, GridId, i, node_veclen;
  double time, time_scalor;

  if ( this->FileName == NULL || this->FileName[0] == '\0'  )
    {
    //this->SetErrorCode(vtkErrorCode::NoFileNameError);
    //vtkErrorMacro(<< "Must specify adG file");
    return 1;
    }
  if(this->GetLengthScale())
    this->myreader->LengthScaleOn();
  else
    this->myreader->LengthScaleOff();

  this->myreader->ScaleChoice = (ScaleOption)this->ScaleChoice;

  this->myreader->SetFileName(this->FileName);

  this->ReadMetaData();
  this->LevelRange[0] = 0;
  this->LevelRange[1] = this->myreader->NumberOfLevels-1;

  this->SetTime(this->GetTime() / this->myreader->TimeScalor);
  double localTime = this->GetTime();
  info->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &localTime, 1);

  for(i=0; i < this->GetNumberOfComponents(); i++)
    {
    this->PointDataArraySelection->AddArray((const char *)this->myreader->Labels[i].label);
    int count = strcspn((const char *)this->myreader->Labels[i].unit, " "); // count how many characters in unit different than " "
    this->myreader->Labels[i].unit[count] = '\0';
    if(count == 0)
      {
      //cerr << "found PointDataArray " << this->myreader->Labels[i].label << "\n";
      }
    else
      {
      //cerr << "found PointDataArray " << this->myreader->Labels[i].label << " [" << this->myreader->Labels[i].unit << "]\n";
      }
    }
  int *blocksPerLevel = new int[this->GetNumberOfLevels()];
  for (levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
    {
    blocksPerLevel[levelId] = this->GridsPerLevels(levelId);
    }
  output->Initialize(this->GetNumberOfLevels(), blocksPerLevel);
  output->SetOrigin(this->myreader->Grids[0].origin);
  if(this->myreader->GetDimensionality() == 2)
    output->SetGridDescription(VTK_XY_PLANE);
  else
    output->SetGridDescription(VTK_XYZ_GRID);
  delete [] blocksPerLevel;

  int current_level = -1;
  std::vector<adG_grid> &grid = this->myreader->Grids;

  // pre-6.0 output->SetNumberOfLevels(this->GetNumberOfLevels());

  double spacing[3];
  int globalBoxId = 0;
  for (levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
    {
    this->myreader->GetSpacing(levelId, spacing);
    output->SetSpacing(levelId, spacing);
cerr << "output->SetSpacing(" << levelId << ", " << spacing[0]<< ", " << spacing[1]<< ", " << spacing[2] << ")\n";
    // old output->SetNumberOfDataSets(levelId, this->GridsPerLevels[levelId]);
    if ( levelId >= this->myreader->MinLevelRead && levelId <= this->myreader->MaxLevelRead )
      {
      // old pre-6.0 output->SetNumberOfDataSets(levelId, this->GridsPerLevels(levelId));

      for (GridId=0; GridId < this->GridsPerLevels(levelId); GridId++)
        {
        output->SetAMRBox(levelId,
                           GridId,
                           this->myreader->Grids[globalBoxId].amrbox);
        //output->SetDataSet(levelId, GridId, 0);
        //cerr << "gid " << globalBoxId << " " << this->myreader->Grids[globalBoxId].amrbox.GetDimensionality() << endl;
        globalBoxId++;
        }
      }
    else
      {
      // what to do? output->SetNumberOfDataSets(levelId, 0);
      globalBoxId += this->GridsPerLevels(levelId);
      //cerr << "skip " << this->GridsPerLevels[levelId]  << " blocks\n";
      }
    }

  vtkDebugMacro(<< "Dimensionality: " << this->myreader->GetDimensionality() << "\nNumberOfComponents: " << this->GetNumberOfComponents() << "\nNumberOfLevels: " << this->GetNumberOfLevels() << "\nNumberOfGrids: " << this->GetNumberOfBlocks());

  output->GenerateParentChildInformation();

  if(output)
    outputVector->GetInformationObject(0)->Set(vtkCompositeDataPipeline::COMPOSITE_DATA_META_DATA(), output);
  else
    outputVector->GetInformationObject(0)->Remove(vtkCompositeDataPipeline::COMPOSITE_DATA_META_DATA());

  //output->GetAMRInfo()->Print(cerr);
  return 1;
} // end of ExecuteInformation


int vtkAMRAmazeReader::RequestData(
  vtkInformation*, vtkInformationVector**, vtkInformationVector* outputVector)
{
  int  n, nb_stars, record;

  herr_t   status;
  hid_t    root_id, level_root_id, attr1;
  this->MinLevelRead = this->LevelRead[0];
  this->MaxLevelRead = this->LevelRead[1];
  //cerr << "read only beetween " << this->MinLevelRead << " and " << this->MaxLevelRead << endl;

  this->UpdateProgress(0.0);
  if(this->LogData == true)
    this->myreader->LogDataOn();
  else
    this->myreader->LogDataOff();

  this->myreader->MakeVariableNames();

  if(this->ShiftedGrid == true)
    {
     this->myreader->ReadHDF5GridsMetaData(true);
     cout << "switching ON the shifted grid info\n";
    }

  if(this->LengthScale == true)
    this->myreader->LengthScaleOn();
  else
    this->myreader->LengthScaleOff();

  vtkInformation* info = outputVector->GetInformationObject(0);
  bool has_block_requests =
    info->Has(vtkCompositeDataPipeline::LOAD_REQUESTED_BLOCKS()) != 0;
  vtkOverlappingAMR* output = vtkOverlappingAMR::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
  cerr << "has_block_requests = " << has_block_requests << endl;
  if(has_block_requests)
  {
    //this->Internal->UpdateIndices = std::set<int>();
    int length = info->Length(vtkCompositeDataPipeline::UPDATE_COMPOSITE_INDICES());
    if (length > 0)
      {
      cerr << "requested " << length << " block" <<endl;
      int* idx = info->Get(vtkCompositeDataPipeline::UPDATE_COMPOSITE_INDICES());
      for(int k=0; k < length; k++)
        cerr << k << ", " << endl;
      cerr << endl;
      //this->Internal->UpdateIndices = std::set<int>(idx, idx+length);
      }
  }
  //double steps[1] = 1.234;
  //output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEPS(), steps, 1);
				  
  //vtkTimerLog *timer0 = vtkTimerLog::New();
  //timer0->StartTimer();
  int piece = info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int numberOfPieces = info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  double localTime = this->GetTime();
  output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), localTime);

  unsigned int numLevels = this->GetNumberOfLevels();
  int i, levelId, g = 0;
#ifdef PARALLEL_DEBUG
  std::ostringstream fname;
  fname << "/tmp/out." << piece << ".txt" << ends;
  ofstream errs;
  errs.open(fname.str().c_str(),ios::app);
  //delete [] fname.str();
  errs << "piece " << piece << " out of " << numberOfPieces << endl;
#endif

  int *blocksPerLevel = new int[this->GetNumberOfLevels()];
  for (levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
    {
    blocksPerLevel[levelId] = this->GridsPerLevels(levelId);
    }
  output->Initialize(this->GetNumberOfLevels(), blocksPerLevel);
  output->SetOrigin(this->myreader->Grids[0].origin);
  if(this->myreader->GetDimensionality() == 2)
    output->SetGridDescription(VTK_XY_PLANE);
  else
    output->SetGridDescription(VTK_XYZ_GRID);
  delete [] blocksPerLevel;

  // Make sure to set the number of levels and number of datasets for
  // each level. This guarantees that this is set properly on all
  // processes, even ones that will not load any data. If this is not
  // done, things that depend on the number of blocks being set
  // properly (such as GenerateVisibilityArrays() will fail)
  if(this->MaxLevelRead == 0)
    { // most likely we did not use the Paraview GUI which forces a valid value
    this->MaxLevelRead = numLevels;
    }

  //output->SetNumberOfLevels(this->GetNumberOfLevels());

  int TotNumBoxes = 0;
  for (int levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
    {
    TotNumBoxes += this->GridsPerLevels(levelId);
#ifdef PARALLEL_DEBUG
    errs << "Level " << levelId << ": " << this->GridsPerLevels(levelId) << " blocks => Tot = " <<   TotNumBoxes << endl;
#endif
    }

// establish what levels are actually read
  int lastLevel = this->MaxLevelRead;
  if (lastLevel > this->GetNumberOfLevels())
    {
    lastLevel = this->GetNumberOfLevels();
    }
  int firstLevel = this->MinLevelRead;
  if (firstLevel < 0)
    {
    firstLevel = 0;
    }
  if (firstLevel > lastLevel)
    {
    firstLevel = lastLevel;
    }

// we now recount, including only the levels to be read
  int endBoxId=0;
  TotNumBoxes = 0;
  for (int levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
    {
     if (levelId <= lastLevel)
       endBoxId += this->GridsPerLevels(levelId);
    if((levelId >= firstLevel) && (levelId <= lastLevel))
      {
      TotNumBoxes += this->GridsPerLevels(levelId);
#ifdef PARALLEL_DEBUG
      errs <<  "NewLevel " << levelId << ": " << this->GridsPerLevels(levelId) << " blocks => Tot = " <<  TotNumBoxes << endl;
#endif
      }
    }
#ifdef PARALLEL_DEBUG
      errs << "firstLevel = " << firstLevel << ", lastLevel = " <<  lastLevel << endl;
#endif

  int c;
  std::vector<adG_grid> &grid = this->myreader->Grids;
  int globalBoxId = 0;
  double spacing[3];

// setup following for all grids, even if not all procs participate to 
  for (int levelId=0; levelId< this->GetNumberOfLevels(); levelId++)
    {
    if (levelId >= firstLevel && levelId <= lastLevel )
      {
      this->myreader->GetSpacing(levelId, spacing);
      output->SetSpacing(levelId, spacing);
      for (int boxId=0; boxId < this->GridsPerLevels(levelId); boxId++, globalBoxId++)
        {
        output->SetAMRBox(levelId, boxId, this->myreader->Grids[globalBoxId].amrbox);
        }
      }
    }
  globalBoxId = 0;
  for (int levelId=0; levelId< this->GetNumberOfLevels(); levelId++)
    {
    if (levelId >= firstLevel && levelId <= lastLevel )
      {
#ifdef PARALLEL_DEBUG
      errs << "Level = " << levelId << ", numBoxes = " <<  this->GridsPerLevels(levelId) << endl;
#endif

      for (int boxId=0; boxId < this->GridsPerLevels(levelId); boxId++, globalBoxId++)
        {
        if( ((endBoxId - globalBoxId) % numberOfPieces) != piece)
          {
	  output->SetDataSet(levelId, boxId, NULL);
#ifdef PARALLEL_DEBUG
      errs << "  skiped BoxId = " << boxId << " (global id = " << globalBoxId << ")\n";
#endif
          continue;
          }
        vtkUniformGrid* ug = this->myreader->ReadUniformGrid(levelId, boxId);
        output->SetDataSet(levelId, boxId, ug);
        //cerr << "gid2 " << globalBoxId << " " << this->myreader->Grids[globalBoxId].amrbox.GetDimensionality() << endl;
        ug->Delete();
#ifdef PARALLEL_DEBUG
        errs << "  read BoxId = " << boxId << " (global id = " << globalBoxId << ")\n";
#endif
        for(c=0; c < this->GetNumberOfComponents(); c++)
          {
          if(this->PointDataArraySelection->ArrayIsEnabled(this->GetPointArrayName(c)))
            { // used to test this->PointDataArraySelection->GetArraySetting(c)
            //vtkDoubleArray* scalars = this->myreader->ReadVar(levelId, boxId, this->myreader->Labels[c]);
            vtkDoubleArray* scalars = this->myreader->ReadVisItVar(globalBoxId, (const char*)this->myreader->Labels[c].label);

            ug->GetPointData()->AddArray(scalars);
            if (!ug->GetPointData()->GetScalars())
              {
              ug->GetPointData()->SetScalars(scalars);
              }
            scalars->Delete();
            } // variable was selected for reading
          }  // end of reading of all components

    this->UpdateProgress (0.85*double(globalBoxId)/double(endBoxId));
      } // all grids at that particular level
#define ALL 1

      } // if this level is included between MinLevel and MaxLevel
    else
      {
      // what to do ? output->SetNumberOfDataSets(levelId, 0);

      globalBoxId += this->GridsPerLevels(levelId);
      }
    } //end of processing all levels
  //cerr << __LINE__ << "Entering GenerateParentChildInformation()\n";
  //output->GenerateParentChildInformation();
  cerr << __LINE__ << "Entering BlankCells()\n";
  vtkAMRUtilities::BlankCells(output);

  //timer0->StopTimer();
  //cerr << "time elapsed to read data arrays: " << timer0->GetElapsedTime() << endl;
  //timer0->Delete();
  
  //vtkTimerLog *timer = vtkTimerLog::New();
  //timer->StartTimer();

  //timer->StopTimer();
  //cerr << "time elapsed to generate visibility arrays: " << timer->GetElapsedTime() << endl;
  //timer->Delete();

#ifdef PARALLEL_DEBUG
    errs.close();
#endif

  if(piece == 0)// should only proc 0 read the stars?
    {
    //cerr << "vtkAMRAmazeReader::RequestData() Load stars\n\n";

      vtkInformation* info = outputVector->GetInformationObject(1);
      vtkMultiBlockDataSet* output2 = vtkMultiBlockDataSet::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
      if (!output2)
        {
        cerr << "did not get a correct output2 object\n";
        }
      else
        {
        nb_stars = vtkAMRAmazeReader::LoadStars(root_id, output2);
        }
    }

  this->UpdateProgress(1.0);
  return 1;
} // RequestData

const char* vtkAMRAmazeReader::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}

int vtkAMRAmazeReader::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

void vtkAMRAmazeReader::SetPointArrayStatus(const char* name, int status)
{
  if(status)
    {
    this->PointDataArraySelection->EnableArray(name);
    }
  else
    {
    this->PointDataArraySelection->DisableArray(name);
    }
}

int vtkAMRAmazeReader::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}

void vtkAMRAmazeReader::Enable(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}

void vtkAMRAmazeReader::Disable(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}

void vtkAMRAmazeReader::EnableAll()
{
  this->PointDataArraySelection->EnableAllArrays();
}

void vtkAMRAmazeReader::DisableAll()
{
  this->PointDataArraySelection->DisableAllArrays();
}
#endif  // OLD stuff
