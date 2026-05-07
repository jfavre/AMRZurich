#include "vtkAMRAmazeReader.h"

#include "vtkByteSwap.h"
#include "vtkCharArray.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkDataArray.h"
#include "vtkCellArray.h"
#include "vtkDataArraySelection.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkErrorCode.h"
#include "vtkOverlappingAMR.h"
#include "vtkInformation.h"
#include "vtkAMRInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkUniformGrid.h"

#include "vtkTimerLog.h"
#include "vtkAMRUtilities.h"
#define PARALLEL_DEBUG 1

#include <vector>
#include <string>
#include <map>
#include <sstream>
#include <iostream>

using namespace std;

vtkStandardNewMacro(vtkAMRAmazeReader);

//
// http://www.paraview.org/ParaView/index.php/Multi-Resolution_Rendering_with_Overlapping_AMR
//
vtkAMRAmazeReader::vtkAMRAmazeReader()
{
  this->FileName = nullptr;
  this->nbstars = 0;
  this->LogDataOn();
  this->DataScaleOn();
  this->ShiftedGridOff();
  this->PointDataArraySelection = vtkDataArraySelection::New();
  this->DebugOff();
  this->MaxLevelWrite = -1;

  this->LevelRead[0] = -1;
  this->LevelRead[1] = -1;

  this->LevelRange[0] = -1;
  this->LevelRange[1] = -1;
  this->SetNumberOfInputPorts(0);
  this->LengthScaleFactor = 1e13; // bogus
  this->SetNumberOfOutputPorts(2);
// this is for port number 1 which we do in all cases.
  vtkMultiBlockDataSet *pd = vtkMultiBlockDataSet::New();
  pd->ReleaseData();
  this->GetExecutive()->SetOutputData(1, pd);
  pd->Delete();
// this is for port number 1 which we do in all cases.
  this->Internal = vtkAMRAmazeReaderInternal::New();
}

vtkAMRAmazeReader::~vtkAMRAmazeReader()
{
  if(this->FileName != NULL)
    {
    delete [] this->FileName;
    this->FileName = nullptr;
    }

  this->PointDataArraySelection->Delete();
  this->Internal->Delete();
}

vtkTypeBool vtkAMRAmazeReader::CanReadFile(const char* fname )
{
  if (! fname )
    return false;
    /*
  cout << "vtkAMRAmazeReader::CanReadFile\n";
  hid_t f_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t root_id = H5Gopen(f_id, "/", H5P_DEFAULT);
  if(H5Lexists(root_id, "/Grid Info", H5P_DEFAULT))
    {
      H5Gclose(root_id);
      H5Fclose(f_id);
      return true;
    }
  else
    {
    H5Gclose(root_id);
    H5Fclose(f_id);
    return false;
    }
    */return true;
}


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

  vtkInformation* info = outputVector->GetInformationObject(0);
  //info->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

  //vtkInformation* info1 = outputVector->GetInformationObject(1);
  //info1->Set(vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(), -1);

  vtkOverlappingAMR *output = static_cast<vtkOverlappingAMR *>(
    info->Get(vtkDataObject::DATA_OBJECT()));
  if (output == nullptr)
  {
    output = vtkOverlappingAMR::New();
    //cout << __LINE__ << " : Got a New vtkOverlappingAMR* = " << output << std::endl;
  }

  int levelId;
  double time, time_scalor;

  if ( this->FileName == NULL || this->FileName[0] == '\0'  )
  {
    //this->SetErrorCode(vtkErrorCode::NoFileNameError);
    //vtkErrorMacro(<< "Must specify adG file");
    return 1;
  }
  if(this->GetLengthScale())
    this->Internal->LengthScaleOn();
  else
    this->Internal->LengthScaleOff();

  this->Internal->ScaleChoice = (ScaleOption)this->ScaleChoice;
  //cerr << __LINE__ << "vtkAMRAmazeReader::RequestInformation() Fname = " << this->FileName << "\n";
  this->Internal->SetFileName(this->FileName);

  this->nbstars = this->Internal->ReadMetaData();

  this->LevelRange[0] = 0;
  this->LevelRange[1] = this->Internal->NumberOfLevels-1;

  double localTime = this->GetTime();
  double timeRange[2] = {localTime, localTime};
  info->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  info->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());

  //cerr << __LINE__ << "Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS() = " << localTime << "\n";
  info->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &localTime, 1);
  output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), localTime);
  info->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
  
  this->SetUpDataArraySelections();
  
  std::vector<unsigned int> blocksPerLevel;
  blocksPerLevel.reserve(this->GetNumberOfLevels());
  for (levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
  {
    blocksPerLevel.push_back(this->GridsPerLevels(levelId));
  }
  output->Initialize(blocksPerLevel);
  output->SetOrigin(this->Internal->Grids[0].layout.origin);
  if(this->Internal->GetDimensionality() == 2)
    output->SetGridDescription(vtkStructuredData::VTK_STRUCTURED_XY_PLANE);
  else
    output->SetGridDescription(vtkStructuredData::VTK_STRUCTURED_XYZ_GRID);

  int current_level = -1;
  std::vector<adG_grid> &grid = this->Internal->Grids;

  double spacing[3];
  int globalBoxId = 0;
  for (levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
  {
    this->Internal->GetSpacing(levelId, spacing);
    output->SetSpacing(levelId, spacing);
    //cout << __LINE__ << ": output->SetSpacing(" << levelId << ", " << spacing[0]<< ", " << spacing[1]<< ", " << spacing[2] << ")\n";
    // old output->SetNumberOfDataSets(levelId, this->GridsPerLevels[levelId]);
    if ( levelId >= this->Internal->MinLevelRead && levelId <= this->Internal->MaxLevelRead )
    {
      for (auto GridId=0; GridId < this->GridsPerLevels(levelId); GridId++)
      {
        output->SetAMRBox(levelId,
                           GridId,
                           this->Internal->Grids[globalBoxId].amrbox);
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

  vtkDebugMacro(<< "\nDimensionality: " << this->Internal->GetDimensionality()
                << "\nNumberOfComponents: " << this->GetNumberOfComponents()
                << "\nNumberOfLevels: " << this->GetNumberOfLevels()
                << "\nNumberOfBlocks: " << this->GetNumberOfBlocks()
                );

  output->GenerateParentChildInformation();

  if(output)
    outputVector->GetInformationObject(0)->Set(vtkCompositeDataPipeline::COMPOSITE_DATA_META_DATA(), output);
  else
    outputVector->GetInformationObject(0)->Remove(vtkCompositeDataPipeline::COMPOSITE_DATA_META_DATA());

  return 1;
} // end of ExecuteInformation

//------------------------------------------------------------------------------
void vtkAMRAmazeReader::GetAMRGridPointData(int blockIdx, vtkUniformGrid* block, const char* field)
{
  vtkDoubleArray* scalars = this->Internal->ReadVisItVar(blockIdx, field);

  block->GetPointData()->AddArray(scalars);
  if (!block->GetPointData()->GetScalars())
  {
    block->GetPointData()->SetScalars(scalars);
  }
  scalars->Delete();
}

//------------------------------------------------------------------------------
void vtkAMRAmazeReader::SetUpDataArraySelections()
{
  for(auto i=0; i < this->GetNumberOfComponents(); i++)
  {
    this->PointDataArraySelection->AddArray((const char *)this->Internal->Labels[i].label);
    int count = strcspn((const char *)this->Internal->Labels[i].unit, " "); // count how many characters in unit different than " "
    this->Internal->Labels[i].unit[count] = '\0';
  }
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::GetBlockLevel(int blockIdx)
{
  return this->Internal->GetBlockLevel(blockIdx);
}
//------------------------------------------------------------------------------
vtkUniformGrid* vtkAMRAmazeReader::GetAMRGrid(int blockIdx)
{
  return this->Internal->GetAMRGrid(blockIdx);
}

int vtkAMRAmazeReader::RequestData(
  vtkInformation*, vtkInformationVector**, vtkInformationVector* outputVector)
{
  int local_nb_stars;

  this->MinLevelRead = this->LevelRead[0];
  this->MaxLevelRead = this->LevelRead[1];
  //cerr << "read only beetween " << this->MinLevelRead << " and " << this->MaxLevelRead << endl;

  this->UpdateProgress(0.0);
  if(this->LogData)
    this->Internal->LogDataOn();
  else
    this->Internal->LogDataOff();

  this->Internal->MakeVariableNames();

  if(this->ShiftedGrid)
  {
    this->Internal->ReadHDF5GridsMetaData(true);
    cout << "switching ON the shifted grid info\n";
  }

  if(this->LengthScale)
    this->Internal->LengthScaleOn();
  else
    this->Internal->LengthScaleOff();

  vtkInformation* info = outputVector->GetInformationObject(0);
  
  double requestedTime = info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
  cout << __LINE__ << ": requestedTime = " << requestedTime << endl;
  int length = info->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

  bool has_block_requests =
    info->Has(vtkCompositeDataPipeline::LOAD_REQUESTED_BLOCKS()) != 0;
  vtkOverlappingAMR* output = vtkOverlappingAMR::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
  //cout << __LINE__ << ": has_block_requests = " << has_block_requests << endl;
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

  int piece = info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  int numberOfPieces = info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), this->Internal->AMAZETime);

  unsigned int numLevels = this->GetNumberOfLevels();

#ifdef PARALLEL_DEBUG
  std::ostringstream fname;
  fname << "/capstor/scratch/cscs/jfavre/Amaze/out." << piece << ".txt";

  std::ofstream errs;
  errs.open(fname.str().c_str(), std::ios::app);

// If opening the primary file fails, fall back to /dev/shm
  if (!errs.is_open())
  {
    std::ostringstream altname;
    altname << "/dev/shm/out." << piece << ".txt";
    errs.clear(); // clear fail state before retrying
    errs.open(altname.str().c_str(), std::ios::app);
    if (!errs.is_open())
    {
      std::cerr << "Error: could not open either output file:\n"
                << "  " << fname.str() << "\n"
                << "  " << altname.str() << std::endl;
    }
  }

  errs << "piece " << piece << " out of " << numberOfPieces << endl;
  //errs << "time = " << localTime  << endl;
#endif

  std::vector<unsigned int> blocksPerLevel;
  blocksPerLevel.reserve(this->GetNumberOfLevels());
  for (auto levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
  {
    blocksPerLevel.push_back(this->GridsPerLevels(levelId));
  }
  output->Initialize(blocksPerLevel);
  output->SetOrigin(this->Internal->Grids[0].layout.origin);
  
  if(this->Internal->GetDimensionality() == 2)
    output->SetGridDescription(vtkStructuredData::VTK_STRUCTURED_XY_PLANE);
  else
    output->SetGridDescription(vtkStructuredData::VTK_STRUCTURED_XYZ_GRID);

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

  std::vector<adG_grid> &grid = this->Internal->Grids;
  int globalBoxId = 0;
  double spacing[3];

// setup following for all grids, even if not all procs participate to 
  for (int levelId=0; levelId< this->GetNumberOfLevels(); levelId++)
  {
    if (levelId >= firstLevel && levelId <= lastLevel )
    {
      this->Internal->GetSpacing(levelId, spacing);
      output->SetSpacing(levelId, spacing);
      for (int boxId=0; boxId < this->GridsPerLevels(levelId); boxId++, globalBoxId++)
      {
        output->SetAMRBox(levelId, boxId, this->Internal->Grids[globalBoxId].amrbox);
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
        vtkUniformGrid* ug = this->GetAMRGrid(this->Internal->FindDomainId(levelId, boxId));
        output->SetDataSet(levelId, boxId, ug);
        //cerr << "gid2 " << globalBoxId << " " << this->Internal->Grids[globalBoxId].grid_nr << endl;
        ug->Delete();
#ifdef PARALLEL_DEBUG
        errs << "  read BoxId = " << boxId << " (global id = " << globalBoxId << ")\n";
#endif
        for(auto c=0; c < this->GetNumberOfComponents(); c++)
        {
          if(this->PointDataArraySelection->ArrayIsEnabled(this->GetPointArrayName(c)))
          {
            this->GetAMRGridPointData(globalBoxId, ug, (const char*)this->Internal->Labels[c].label);
          } // variable was selected for reading
        }  // end of reading of all components

        this->UpdateProgress (0.85*double(globalBoxId)/double(endBoxId));
      } // all grids at that particular level

    } // if this level is included between MinLevel and MaxLevel
    else
    {
      // what to do ? output->SetNumberOfDataSets(levelId, 0);
      globalBoxId += this->GridsPerLevels(levelId);
    }
  } //end of processing all levels
  //cerr << __LINE__ << "Entering GenerateParentChildInformation()\n";
  //output->GenerateParentChildInformation();
  vtkAMRUtilities::BlankCells(output);

#ifdef PARALLEL_DEBUG
  errs.close();
#endif
  // the second output holds the Stars Polydata
  info = outputVector->GetInformationObject(1);
  vtkMultiBlockDataSet* output2 = vtkMultiBlockDataSet::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
  if(piece == 0) // only MPI task 0 reads the stars
  {
    local_nb_stars = vtkAMRAmazeReader::LoadStars(output2);
    if(local_nb_stars != this->nbstars)
    {
      cerr << "something went wrong counting stars\n";
      cerr << __LINE__ << "nb_stars = " << local_nb_stars << endl;
    }
  }
  else
  {
    for(auto i=0; i < this->nbstars; i++)
      output2->SetBlock(i, nullptr);
  }

  this->UpdateProgress(1.0);
  cout << "=================================================================\n";
  return 1;
} // RequestData

int vtkAMRAmazeReader::LoadStars(vtkMultiBlockDataSet* SpherSymStars)
{
  this->Internal->BuildStars();
  int nb_stars = this->Internal->NumberOfSphericallySymmetricStars + this->Internal->NumberOfAxisSymmetricStars;
  for(auto i=0; i < nb_stars; i++)
  {
    SpherSymStars->SetBlock(i, this->Internal->Stars[i]);
    SpherSymStars->GetMetaData((unsigned int)i)->Set(vtkCompositeDataSet::NAME(), this->Internal->Stars[i]->GetFieldData()->GetArray(0)->GetName());
  }

  return nb_stars;
}

//----------------------------------------------------------------------------
void vtkAMRAmazeReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);

  os << "FileName: " << (this->FileName? this->FileName:"(none)") << "\n";
  os << indent << "Dimensionality: " << this->Internal->GetDimensionality() << endl;
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

int vtkAMRAmazeReader::FillOutputPortInformation(int port, vtkInformation* info)
{
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkOverlappingAMR");
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  return 1;
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

// Get the second output which contains the stars
//vtkPolyData* vtkAMRAmazeReader::GetStarsOutput()
vtkMultiBlockDataSet* vtkAMRAmazeReader::GetStarsOutput()
{
  if (this->GetNumberOfOutputPorts() < 3)
  {
    return nullptr;
  }
  return vtkMultiBlockDataSet::SafeDownCast(this->GetExecutive()->GetOutputData(2));
}
