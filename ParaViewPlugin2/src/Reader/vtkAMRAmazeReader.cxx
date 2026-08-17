// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
#include "vtkAMRAmazeReader.h"
#include "vtkAMRBox.h"
#include "vtkByteSwap.h"
#include "vtkDataArraySelection.h"
#include "vtkObjectFactory.h"
#include "vtkOverlappingAMR.h"
#include "vtkUniformGrid.h"

#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkLongLongArray.h"
#include "vtkNew.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkShortArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedShortArray.h"

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <map>
#include <vector>
#include <iostream>

#include "vtkAMRAmazeReaderInternal.h"

vtkStandardNewMacro(vtkAMRAmazeReader);

//------------------------------------------------------------------------------
vtkAMRAmazeReader::vtkAMRAmazeReader()
{
  this->nbstars = 0;

  this->DataScaleOn();

  //if(this->PointDataArraySelection == nullptr)
    //this->PointDataArraySelection = vtkDataArraySelection::New();
  this->DebugOff();

  this->SetNumberOfInputPorts(0);

  this->SetNumberOfOutputPorts(2);
// this is for port number 1 which we do in all cases.
  vtkMultiBlockDataSet *pd = vtkMultiBlockDataSet::New();
  pd->ReleaseData();
  this->GetExecutive()->SetOutputData(1, pd);
  pd->Delete();

  this->IsReady = false;
// this is for port number 1 which we do in all cases.
  this->Internal = std::make_unique<vtkAMRAmazeReaderInternal>();
  this->Initialize();
}

//------------------------------------------------------------------------------
vtkAMRAmazeReader::~vtkAMRAmazeReader()
{
  if(this->FileName != NULL)
    {
    delete [] this->FileName;
    this->FileName = nullptr;
    }
}

int vtkAMRAmazeReader::FillOutputPortInformation(int port, vtkInformation* info)
{
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkOverlappingAMR");
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  return 1;
}

int vtkAMRAmazeReader::RequestInformation(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector)
{
  vtkInformation* info = outputVector->GetInformationObject(0);
  
  int status = vtkAMRBaseReader::RequestInformation(request, inputVector, outputVector);
  
  info->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  info->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());
  double localTime = this->Internal->AMAZETime;
  double timeRange[2] = {localTime, localTime};
  info->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &localTime, 1);
  info->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);
  
  return status;
}

int vtkAMRAmazeReader::RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector)
{
  int status = vtkAMRBaseReader::RequestData(request, inputVector, outputVector);
    
  double localTime = this->Internal->AMAZETime;
  this->Metadata->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), localTime);
    
  // the first output holds the AMR structure
  vtkInformation* info0 = outputVector->GetInformationObject(0);
  vtkOverlappingAMR* amr = vtkOverlappingAMR::SafeDownCast(info0->Get(vtkDataObject::DATA_OBJECT()));
  bool res = amr->CheckValidity();

  // the second output holds the Stars Polydata
  vtkInformation* info = outputVector->GetInformationObject(1);
  vtkMultiBlockDataSet* output2 = vtkMultiBlockDataSet::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
  if(1) // only MPI task 0 reads the stars
  {
    int local_nb_stars = vtkAMRAmazeReader::LoadStars(output2);
    if(local_nb_stars != this->nbstars)
    {
      std::cerr << "something went wrong counting stars\n";
      std::cerr << __LINE__ << "nb_stars = " << local_nb_stars << endl;
    }
  }
  else
  {
    for(auto i=0; i < this->nbstars; i++)
      output2->SetBlock(i, nullptr);
  }
  return status;
};

vtkTypeBool vtkAMRAmazeReader::CanReadFile(const char* fname ) const
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


//------------------------------------------------------------------------------
void vtkAMRAmazeReader::PrintSelf(std::ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//------------------------------------------------------------------------------
void vtkAMRAmazeReader::SetFileName(const char* fileName)
{
  assert("pre: Internal Amaze Reader is nullptr" && (this->Internal != nullptr));

  if (fileName && strcmp(fileName, "") != 0 &&
    ((this->FileName == nullptr) || strcmp(fileName, this->FileName) != 0))
  {
    if (this->FileName)
    {
      delete[] this->FileName;
      this->FileName = nullptr;
      this->Internal->SetFileName(nullptr);
    }

    this->FileName = new char[strlen(fileName) + 1];
    strcpy(this->FileName, fileName);
    this->FileName[strlen(fileName)] = '\0';

    this->IsReady = true;
    this->Internal->SetFileName(this->FileName);
    this->LoadedMetaData = false;

    this->SetUpDataArraySelections();
    this->InitializeArraySelections();
  }

  this->Modified();
}

//------------------------------------------------------------------------------
void vtkAMRAmazeReader::ReadMetaData()
{
  assert("pre: Internal Amaze Reader is nullptr" && (this->Internal != nullptr));
  this->Internal->ReadMetaData();
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::GetBlockLevel(int blockIdx)
{
  return this->Internal->GetBlockLevel(blockIdx);
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::GetNumberOfBlocks()
{
  assert("pre: Internal Amaze Reader is nullptr" && (this->Internal != nullptr));
  if (!this->IsReady)
  {
    return 0;
  }
  
  return Internal->NumberOfGrids;
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::GetNumberOfLevels()
{
  assert("pre: Internal Amaze Reader is nullptr" && (this->Internal != nullptr));
  if (!this->IsReady)
  {
    return 0;
  }
  
  this->Internal->ReadMetaData();
  //std::cerr << __LINE__ << " :vtkAMRAmazeReader::GetNumberOfLevels() = " <<  this->Internal->NumberOfLevels << "\n";
  return (this->Internal->NumberOfLevels);
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::FillMetaData()
{
  this->nbstars = this->Internal->ReadMetaData();

  this->SetUpDataArraySelections();

  std::vector<unsigned int> blocksPerLevel;
  blocksPerLevel.reserve(this->GetNumberOfLevels());
  for (int levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
  {
    //std::cout << "Level " <<  levelId << " has " << this->Internal->GridsPerLevels(levelId) << " blocks\n";
    blocksPerLevel.push_back(this->Internal->GridsPerLevels(levelId));
  }
  this->Metadata->Initialize(blocksPerLevel);
  this->Metadata->SetOrigin(this->Internal->Grids[0].layout.origin);
  if(this->Internal->Dimensionality == 2)
    this->Metadata->SetGridDescription(vtkStructuredData::VTK_STRUCTURED_XY_PLANE);
  else
    this->Metadata->SetGridDescription(vtkStructuredData::VTK_STRUCTURED_XYZ_GRID);

  double spacing[3];
  int globalBoxId = 0, refratio=-1;
  for (int levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
  {
    this->Internal->GetSpacing(levelId, spacing);
    this->Metadata->SetSpacing(levelId, spacing);
    if (levelId == this->GetNumberOfLevels() - 1)  // we are at finest level of refinement
    {
      refratio = 1;
    }
    else
    {
      refratio = this->Internal->Levels[levelId+1].RefRatio; // note the offset by +1
    }
    this->Metadata->SetRefinementRatio(levelId, refratio);

    //std::cerr << "Level " <<  levelId << ": RefRatio = " << refratio << ", spacing = [ " << spacing[0] << ", " << spacing[1] <<  ", " << spacing[2] << "]\n";

    for (auto GridId=0; GridId < this->Internal->GridsPerLevels(levelId); GridId++)
    {
      this->Metadata->SetAMRBox(levelId, GridId, this->Internal->Grids[globalBoxId].amrbox);
      this->Metadata->SetAMRBlockSourceIndex(levelId, GridId, globalBoxId);
      globalBoxId++;
    }
  }

  vtkDebugMacro(<< "\nDimensionality: "     << this->Internal->Dimensionality
                << "\nNumberOfComponents: " << this->Internal->NumberOfComponents
                << "\nNumberOfLevels: "     << this->GetNumberOfLevels()
                << "\nNumberOfBlocks: "     << this->GetNumberOfBlocks()
                );

  this->Metadata->GenerateParentChildInformation();
  return 1;
}

//------------------------------------------------------------------------------
vtkUniformGrid* vtkAMRAmazeReader::GetAMRGrid(int blockIdx)
{
  //std::cerr << __LINE__ << " :vtkAMRAmazeReader::GetAMRGrid(" <<  blockIdx << ")\n";
  return this->Internal->GetAMRGrid(blockIdx);
}

//------------------------------------------------------------------------------
void vtkAMRAmazeReader::GetAMRGridPointData(int blockIdx, vtkUniformGrid* block, const char* field)
{
  //std::cerr << __LINE__ << " :GetAMRGridPointData(" <<  blockIdx << ", " << field << ")\n";
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
  for(auto i=0; i < this->Internal->NumberOfComponents; i++)
  {
    // std::cerr << __LINE__ << "PointDataArraySelection->AddArray(" <<  this->Internal->Labels[i].label << ")\n";
    // all arrays are added as disabled.
    this->PointDataArraySelection->AddArray((const char *)this->Internal->Labels[i].label, false);
    int count = strcspn((const char *)this->Internal->Labels[i].unit, " "); // count how many characters in unit different than " "
    this->Internal->Labels[i].unit[count] = '\0';
  }
  this->Internal->MakeVariableNames();
}

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
