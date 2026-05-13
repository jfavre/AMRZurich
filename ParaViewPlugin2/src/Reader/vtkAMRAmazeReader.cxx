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

#include "vtk_hdf5.h"

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

//------------------------------------------------------------------------------
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
  return Internal->NumberOfGrids;
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::GetNumberOfLevels()
{
return 0;
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::FillMetaData()
{
return 0;
}

//------------------------------------------------------------------------------
vtkUniformGrid* vtkAMRAmazeReader::GetAMRGrid(int blockIdx)
{
  return this->Internal->GetAMRGrid(blockIdx);
}

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
  for(auto i=0; i < this->Internal->GetNumberOfComponents(); i++)
  {
    this->PointDataArraySelection->AddArray((const char *)this->Internal->Labels[i].label);
    int count = strcspn((const char *)this->Internal->Labels[i].unit, " "); // count how many characters in unit different than " "
    this->Internal->Labels[i].unit[count] = '\0';
  }
}

