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
#include "vtkDataSet.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
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
  this->FileName = nullptr;
  this->nbstars = 0;
  this->LogDataOn();
  this->DataScaleOn();
  this->ShiftedGridOff();
  if(this->PointDataArraySelection == nullptr)
    this->PointDataArraySelection = vtkDataArraySelection::New();
  this->DebugOff();
  this->MaxLevelWrite = -1;

  this->LevelRead[0] = -1;
  this->LevelRead[1] = -1;

  this->LevelRange[0] = -1;
  this->LevelRange[1] = -1;
  this->SetNumberOfInputPorts(0);
  this->LengthScaleFactor = 1e13; // bogus
#ifdef MULTI_PORTS
  this->SetNumberOfOutputPorts(2);
// this is for port number 1 which we do in all cases.
  vtkMultiBlockDataSet *pd = vtkMultiBlockDataSet::New();
  pd->ReleaseData();
  this->GetExecutive()->SetOutputData(1, pd);
  pd->Delete();
#endif
  this->IsReady = false;
// this is for port number 1 which we do in all cases.
  this->Internal = new vtkAMRAmazeReaderInternal();
}

//------------------------------------------------------------------------------
vtkAMRAmazeReader::~vtkAMRAmazeReader()
{
  if(this->FileName != NULL)
    {
    delete [] this->FileName;
    this->FileName = nullptr;
    }

  //this->PointDataArraySelection->Delete();
  delete this->Internal;
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

    //this->SetUpDataArraySelections();
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
  std::cerr << __LINE__ << " :vtkAMRAmazeReader::FillMetaData()\n";
  this->nbstars = this->Internal->ReadMetaData();

  double localTime = this->Internal->AMAZETime;
  this->Metadata->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), localTime);

  this->SetUpDataArraySelections();

  std::vector<unsigned int> blocksPerLevel;
  blocksPerLevel.reserve(this->GetNumberOfLevels());
  for (int levelId=0; levelId < this->GetNumberOfLevels(); levelId++)
  {
    std::cerr << "Level " <<  levelId << " has " << this->Internal->GridsPerLevels(levelId) << " blocks\n";
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

    std::cerr << "Level " <<  levelId << ": RefRatio = " << refratio << ", spacing = [ " << spacing[0] << ", " << spacing[1] <<  ", " << spacing[2] << "]\n";

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
    std::cerr << __LINE__ << " :vtkAMRAmazeReader::GetAMRGrid(" <<  blockIdx << ")\n";
    return this->Internal->GetAMRGrid(blockIdx);
}

//------------------------------------------------------------------------------
void vtkAMRAmazeReader::GetAMRGridPointData(int blockIdx, vtkUniformGrid* block, const char* field)
{
  std::cerr << __LINE__ << " :GetAMRGridPointData(" <<  blockIdx << ", " << field << ")\n";
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
    std::cerr << __LINE__ << "PointDataArraySelection->AddArray(" <<  this->Internal->Labels[i].label << ")\n";
    // all arrays are added as disabled.
    this->PointDataArraySelection->AddArray((const char *)this->Internal->Labels[i].label);
    int count = strcspn((const char *)this->Internal->Labels[i].unit, " "); // count how many characters in unit different than " "
    this->Internal->Labels[i].unit[count] = '\0';
  }
  this->Internal->MakeVariableNames();
}

