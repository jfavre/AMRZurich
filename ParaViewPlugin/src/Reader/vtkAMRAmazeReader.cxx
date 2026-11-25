// SPDX-FileCopyrightText: Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
// SPDX-License-Identifier: BSD-3-Clause
#include "vtkAMRAmazeReader.h"

#include "vtkAMRAmazeReaderInternal.h"

#include "vtkDataArraySelection.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkNew.h"
#include "vtkObjectFactory.h"
#include "vtkOverlappingAMR.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUniformGrid.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <vector>

vtkStandardNewMacro(vtkAMRAmazeReader);

//------------------------------------------------------------------------------
vtkAMRAmazeReader::vtkAMRAmazeReader()
{
  this->DebugOff();
  this->SetNumberOfInputPorts(0);

  // Output port 0 holds the AMR structure, output port 1 holds the stars
  this->SetNumberOfOutputPorts(2);
  vtkMultiBlockDataSet* starsOutput = vtkMultiBlockDataSet::New();
  starsOutput->ReleaseData();
  this->GetExecutive()->SetOutputData(1, starsOutput);
  starsOutput->Delete();

  this->Internal = std::make_unique<vtkAMRAmazeReaderInternal>();
  this->Initialize();
}

//------------------------------------------------------------------------------
vtkAMRAmazeReader::~vtkAMRAmazeReader() = default;

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::FillOutputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
  {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkOverlappingAMR");
  }
  else if (port == 1)
  {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  }
  return 1;
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::RequestInformation(vtkInformation* request,
                                          vtkInformationVector** inputVector,
                                          vtkInformationVector* outputVector)
{
  if (!this->FileName)
  {
    // nothing to do until a file has been selected (ParaView triggers
    // information updates as soon as the reader proxy is created)
    return 1;
  }

  vtkInformation* info = outputVector->GetInformationObject(0);

  int status = vtkAMRBaseReader::RequestInformation(request, inputVector, outputVector);

  const double localTime = this->Internal->AMAZETime;
  const double timeRange[2] = { localTime, localTime };
  info->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
  info->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());
  info->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(), &localTime, 1);
  info->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange, 2);

  return status;
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::RequestData(vtkInformation* request,
                                   vtkInformationVector** inputVector,
                                   vtkInformationVector* outputVector)
{
  int status = vtkAMRBaseReader::RequestData(request, inputVector, outputVector);

  // the first output holds the AMR structure
  vtkInformation* info0 = outputVector->GetInformationObject(0);
  vtkOverlappingAMR* amr = vtkOverlappingAMR::SafeDownCast(info0->Get(vtkDataObject::DATA_OBJECT()));

  const double omega = this->Internal->FixedBinaryOrbit.Omega; // if read from file, never 0.0
  if (omega > 0.0)
  {
    vtkNew<vtkDoubleArray> com;
    com->SetName("CenterOfMass");
    com->SetNumberOfComponents(3);
    com->InsertNextValue(this->Internal->FixedBinaryOrbit.CenterofMass_x);
    com->InsertNextValue(this->Internal->FixedBinaryOrbit.CenterofMass_y);
    com->InsertNextValue(this->Internal->FixedBinaryOrbit.CenterofMass_z);
    amr->GetFieldData()->AddArray(com);

    vtkNew<vtkDoubleArray> omegaArray;
    omegaArray->SetName("Omega");
    omegaArray->SetNumberOfComponents(1);
    omegaArray->InsertNextValue(omega);
    amr->GetFieldData()->AddArray(omegaArray);
  }

  // if read from file, never 0.0
  const double meanParticleMass = this->Internal->thermodynamics_eos.Mean_Particle_Mass;
  if (meanParticleMass > 0.0)
  {
    vtkNew<vtkDoubleArray> eos;
    eos->SetName("Thermodynamics_EOS");
    eos->SetNumberOfComponents(5);
    eos->InsertNextValue(this->Internal->thermodynamics_eos.Adiabatic_Exponent);
    eos->InsertNextValue(this->Internal->thermodynamics_eos.Specific_Gasconstant);
    eos->InsertNextValue(this->Internal->thermodynamics_eos.Scaled_Specific_Gasconstant);
    eos->InsertNextValue(this->Internal->thermodynamics_eos.Scale_Factor);
    eos->InsertNextValue(meanParticleMass);
    amr->GetFieldData()->AddArray(eos);
  }

  // the second output holds the stars poly data
  vtkInformation* info1 = outputVector->GetInformationObject(1);
  vtkMultiBlockDataSet* starsOutput =
    vtkMultiBlockDataSet::SafeDownCast(info1->Get(vtkDataObject::DATA_OBJECT()));

  const int local_nb_stars = this->LoadStars(starsOutput);
  if (local_nb_stars != this->NumberOfStars)
  {
    std::cerr << "something went wrong counting stars: expected " << this->NumberOfStars
              << ", read " << local_nb_stars << "\n";
  }

  // star properties published as FieldData to be used in programmable filters
  const int nbStarInfos =
    std::min(this->NumberOfStars, static_cast<int>(this->Internal->newstars.size()));
  for (int s = 0; s < nbStarInfos; s++)
  {
    const newstar& star = this->Internal->newstars[s];
    vtkNew<vtkDoubleArray> starInfo;
    starInfo->SetName(star.InteractionModel);
    starInfo->SetNumberOfComponents(1 + 3 + 3);
    starInfo->InsertNextValue(star.Mass);
    starInfo->InsertNextValue(star.Position[0]);
    starInfo->InsertNextValue(star.Position[1]);
    starInfo->InsertNextValue(star.Position[2]);
    starInfo->InsertNextValue(star.Velocity[0]);
    starInfo->InsertNextValue(star.Velocity[1]);
    starInfo->InsertNextValue(star.Velocity[2]);
    amr->GetFieldData()->AddArray(starInfo);
  }

  return status;
}

//------------------------------------------------------------------------------
vtkTypeBool vtkAMRAmazeReader::CanReadFile(const char* fname) const
{
  return fname && H5Fis_hdf5(fname) > 0;
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

  if (fileName && *fileName != '\0' &&
      (this->FileName == nullptr || strcmp(fileName, this->FileName) != 0))
  {
    delete[] this->FileName;
    this->FileName = nullptr;
    this->Internal->SetFileName(nullptr);

    const size_t length = strlen(fileName);
    this->FileName = new char[length + 1];
    strcpy(this->FileName, fileName);

    this->IsReady = true;
    this->Internal->SetFileName(this->FileName);

    this->SetUpDataArraySelections();
    this->InitializeArraySelections();
  }

  this->Modified();
}

//------------------------------------------------------------------------------
void vtkAMRAmazeReader::ReadMetaData()
{
  assert("pre: Internal Amaze Reader is nullptr" && (this->Internal != nullptr));
  if (this->Internal->ReadMetaData() >= 0)
  {
    this->LoadedMetaData = true;
  }
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
  return this->IsReady ? this->Internal->NumberOfGrids : 0;
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::GetNumberOfLevels()
{
  assert("pre: Internal Amaze Reader is nullptr" && (this->Internal != nullptr));
  if (!this->IsReady)
  {
    return 0;
  }

  if (!this->LoadedMetaData)
  {
    this->Internal->ReadMetaData();
    this->LoadedMetaData = true;
  }

  return this->Internal->NumberOfLevels;
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::FillMetaData()
{
  if (this->NumberOfStars != -1) // the metadata has already been read
  {
    return 1;
  }

  this->NumberOfStars = this->Internal->ReadMetaData();
  this->SetUpDataArraySelections();

  std::vector<unsigned int> blocksPerLevel;
  blocksPerLevel.reserve(this->GetNumberOfLevels());
  for (int levelId = 0; levelId < this->GetNumberOfLevels(); levelId++)
  {
    blocksPerLevel.push_back(this->Internal->GridsPerLevels(levelId));
  }

  this->Metadata->Initialize(blocksPerLevel);
  this->Metadata->SetOrigin(this->Internal->Grids[0].layout.origin);
  if (this->Internal->Dimensionality == 2)
  {
    this->Metadata->SetGridDescription(vtkStructuredData::VTK_STRUCTURED_XY_PLANE);
  }
  else
  {
    this->Metadata->SetGridDescription(vtkStructuredData::VTK_STRUCTURED_XYZ_GRID);
  }

  double spacing[3];
  int globalBoxId = 0;
  for (int levelId = 0; levelId < this->GetNumberOfLevels(); levelId++)
  {
    this->Internal->GetSpacing(levelId, spacing);
    this->Metadata->SetSpacing(levelId, spacing);

    // the refinement ratio of a level is the ratio towards the next finer level;
    // the finest level gets 1
    const int refratio = levelId == this->GetNumberOfLevels() - 1
      ? 1
      : this->Internal->Levels[levelId + 1].RefRatio;
    this->Metadata->SetRefinementRatio(levelId, refratio);

    //std::cerr << "Level " <<  levelId << ": RefRatio = " << refratio
    //            << ", spacing = [ " << spacing[0] << ", " << spacing[1] <<  ", " << spacing[2] << "]\n";

    const int nbBlocks = this->Internal->GridsPerLevels(levelId);
    for (int gridId = 0; gridId < nbBlocks; gridId++)
    {
      this->Metadata->SetAMRBox(levelId, gridId, this->Internal->Grids[globalBoxId].amrbox);
      this->Metadata->SetAMRBlockSourceIndex(levelId, gridId, globalBoxId);
      globalBoxId++;
    }
  }

  vtkDebugMacro(<< "\nDimensionality: " << this->Internal->Dimensionality
                << "\nNumberOfComponents: " << this->Internal->NumberOfComponents
                << "\nNumberOfLevels: " << this->GetNumberOfLevels()
                << "\nNumberOfBlocks: " << this->GetNumberOfBlocks());

  this->Metadata->GenerateParentChildInformation();
  this->Metadata->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), this->Internal->AMAZETime);

  return 1;
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
  if (!scalars)
  {
    return;
  }

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
  for (int i = 0; i < this->Internal->NumberOfComponents; i++)
  {
    // all arrays are added as disabled
    this->PointDataArraySelection->AddArray(this->Internal->Labels[i].label, false);

    // keep only the first word of the unit string
    const int count = strcspn(this->Internal->Labels[i].unit, " ");
    this->Internal->Labels[i].unit[count] = '\0';
  }
  this->Internal->MakeVariableNames();
}

//------------------------------------------------------------------------------
int vtkAMRAmazeReader::LoadStars(vtkMultiBlockDataSet* starsDataSet)
{
  this->Internal->BuildStars();

  const int nbStars =
    this->Internal->NumberOfSphericallySymmetricStars + this->Internal->NumberOfAxisSymmetricStars;
  for (int i = 0; i < nbStars; i++)
  {
    starsDataSet->SetBlock(i, this->Internal->Stars[i]);
    starsDataSet->GetMetaData(i)->Set(vtkCompositeDataSet::NAME(),
                                      this->Internal->Stars[i]->GetFieldData()->GetArray(0)->GetName());
  }

  return nbStars;
}
