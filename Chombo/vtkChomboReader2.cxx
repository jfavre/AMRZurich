/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkChomboReader2.cxx,v $
  Language:  C++
  Date:      $Date: 2005/05/19 18:54:22 $
  Version:   $Revision: 1.1 $

  Copyright (c) 1993-2002 Ken Martin, Will Schroeder, Bill Lorensen 
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkChomboReader2.h"

#include "vtkAMRBox.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkCommand.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkHierarchicalBoxDataSet.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkObjectFactory.h"
#include "vtkPolyData.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUniformGrid.h"
#include "vtkDataArraySelection.h"

#include <vtkstd/vector>
#include <vtkstd/string>

#include <hdf5.h>

vtkCxxRevisionMacro(vtkChomboReader2, "$Revision: 1.1 $");
vtkStandardNewMacro(vtkChomboReader2);

struct vtkChomboReader2Internals
{
  vtkChomboReader2Internals(vtkChomboReader2* reader);
  ~vtkChomboReader2Internals();

  int GetIntAttribute( hid_t locID, const char* attrName, int& val );
  int GetRealTypeAttribute( hid_t locID, const char* attrName, double& val );
  int GetStringAttribute( hid_t locID, const char* attrName, vtkstd::string& val );
  int ReadBoxes( hid_t locID , vtkstd::vector<vtkAMRBox>& boxes );
  void CreateBoxDataType();

  hid_t BoxDataType;
  vtkChomboReader2* Reader;
  typedef vtkstd::vector<vtkAMRBox> LevelBoxesType;
  typedef vtkstd::vector<LevelBoxesType> BoxesType;
  BoxesType Boxes;

  typedef vtkstd::vector<hsize_t> LevelOffsetsType;
  typedef vtkstd::vector<LevelOffsetsType> OffsetsType;
  OffsetsType Offsets;

  vtkstd::vector<double> DXs;
  vtkstd::vector<vtkstd::string> ComponentNames;
};

//----------------------------------------------------------------------------
vtkChomboReader2Internals::vtkChomboReader2Internals(vtkChomboReader2* reader)
{
  this->BoxDataType = -1;
  this->Reader = reader;
}

//----------------------------------------------------------------------------
vtkChomboReader2Internals::~vtkChomboReader2Internals()
{
  if ( this->BoxDataType >= 0 )
    {
    H5Tclose( this->BoxDataType );
    }
}

//----------------------------------------------------------------------------
void vtkChomboReader2Internals::CreateBoxDataType()
{
  if ( this->BoxDataType >= 0)
    {
    int size = H5Tget_size(this->BoxDataType);
    int dimensionality = size/2/sizeof(int);
    if (dimensionality == this->Reader->Dimensionality)
      {
      return;
      }
    H5Tclose( this->BoxDataType );
    }

  hid_t boxType = H5Tcreate( H5T_COMPOUND, 
                             2*this->Reader->Dimensionality*sizeof(int) );
  char loCornerComponentNames[3][5] = {"lo_i", "lo_j", "lo_k"};
  char hiCornerComponentNames[3][5] = {"hi_i", "hi_j", "hi_k"};

  int i;
  for( i=0; i<this->Reader->Dimensionality; i++ )
    {
    H5Tinsert( boxType, 
               loCornerComponentNames[i],
               i*sizeof(int),
               H5T_NATIVE_INT );
    }
  for( i=0; i<this->Reader->Dimensionality; i++ )
    {
    H5Tinsert( boxType, 
               hiCornerComponentNames[i],
               (this->Reader->Dimensionality+i) * sizeof(int),
               H5T_NATIVE_INT );
    }
  this->BoxDataType = boxType;
}

//----------------------------------------------------------------------------
int vtkChomboReader2Internals::GetRealTypeAttribute( 
  hid_t locID, const char* attrName, double& val )
{
  hid_t attribute = H5Aopen_name( locID, attrName );
  if ( attribute < 0 )
    {
    vtkGenericWarningMacro( << "Failed to open attribute " << attrName );
    return 0;
    }

  if ( this->Reader->RealType == vtkChomboReader2::FLOAT )
    {
    float tmp;
    if( H5Aread( attribute, H5T_NATIVE_FLOAT, &tmp ) != 0 )
      {
      vtkGenericWarningMacro( << "Failed to read attribute " << attrName );
      H5Aclose( attribute );
      return 0;
      }
    val = tmp;
    }
  else
    {
    double tmp;
    if( H5Aread( attribute, H5T_NATIVE_DOUBLE, &tmp ) != 0 )
      {
      vtkGenericWarningMacro( << "Failed to read attribute " << attrName );
      H5Aclose( attribute );
      return 0;
      }
    val = tmp;
    }
  H5Aclose( attribute );
  return 1;
}

//----------------------------------------------------------------------------
int vtkChomboReader2Internals::GetIntAttribute( 
  hid_t locID, const char* attrName, int& val )
{
  hid_t attribute = H5Aopen_name( locID, attrName );
  if ( attribute < 0 )
    {
    vtkGenericWarningMacro( << "Failed to open attribute " << attrName);
    H5Aclose( attribute );
    return 0;
    }

  if( H5Aread( attribute, H5T_NATIVE_INT, &val ) != 0 )
    {
    vtkGenericWarningMacro( << "Failed to read attribute " << attrName );
    H5Aclose( attribute );
    return 0;
    }
  H5Aclose( attribute );
  return 1;
}

//----------------------------------------------------------------------------
int vtkChomboReader2Internals::GetStringAttribute( 
  hid_t locID, const char* attrName, vtkstd::string& val)
{
  hid_t attribute = H5Aopen_name( locID, attrName );
  if ( attribute < 0 )
    {
    vtkGenericWarningMacro( << "Failed to open attribute " << attrName);
    return 0;
    }

  hid_t atype  = H5Aget_type(attribute);
  size_t size = H5Tget_size(atype);
  char* str = new char[size+1];
  if ( H5Aread(attribute, atype, str ) != 0)
    {
    vtkGenericWarningMacro( << "Failed to read attribute " << attrName );
    H5Aclose( attribute );
    delete [] str;
    str = 0;
    return 0;
    }
  str[size] = 0;
  val = str;
  delete[] str;
  return 1;
}


//----------------------------------------------------------------------------
int vtkChomboReader2Internals::ReadBoxes( 
  hid_t locID , vtkstd::vector<vtkAMRBox>& boxes )
{
  hid_t boxdataset, boxdataspace, memdataspace;
  hsize_t dims[1], maxdims[1];

  boxdataset = H5Dopen( locID, "boxes" );
  if( boxdataset < 0 )
    {
    vtkGenericWarningMacro( << "Failed to open dataset boxes"  );
    return 0;
    }
  boxdataspace = H5Dget_space( boxdataset );
  if( boxdataspace < 0 )
    {
    vtkGenericWarningMacro( << "Failed to open dataspace for dataset boxes" );
    H5Dclose( boxdataset );
    return 0;
    }

  H5Sget_simple_extent_dims( boxdataspace, dims, maxdims );
  memdataspace = H5Screate_simple( 1, dims, 0 );

  int numBoxes = dims[0];

  boxes.resize(numBoxes);
  vtkstd::vector< int > buf( numBoxes*2*this->Reader->Dimensionality );
  this->CreateBoxDataType();
  if ( H5Dread(boxdataset, this->BoxDataType, memdataspace,
               boxdataspace, H5P_DEFAULT, &buf[0] ) < 0 )
    {
    vtkGenericWarningMacro( << "Could not read box information." );
    H5Dclose( boxdataset );
    H5Sclose( boxdataspace );
    H5Sclose( memdataspace );
    return 0;
    }
    
  for(int i=0; i<numBoxes; i++)
    {
    int idx = i*2*this->Reader->Dimensionality;
    boxes[i] = vtkAMRBox(this->Reader->Dimensionality, 
                         &buf[idx], 
                         &buf[idx+this->Reader->Dimensionality]);
    }

  H5Dclose( boxdataset );
  H5Sclose( boxdataspace );
  H5Sclose( memdataspace );

  return 1;
}

//----------------------------------------------------------------------------
vtkChomboReader2::vtkChomboReader2()
{
  this->FileName = 0;

  this->Dimensionality = 0;
  this->NumberOfLevels = 0;
  this->NumberOfComponents = 0;

  this->LevelRange[0] =  0;
  this->LevelRange[1] = -1;
  this->UnderSampleScaleFactor = 2;
  this->TotNumBoxes = 0;
  
  this->PointDataArraySelection = vtkDataArraySelection::New();
  
  this->SetNumberOfInputPorts(0);

  this->SetNumberOfOutputPorts(2);
  vtkPolyData *pd = vtkPolyData::New();
  pd->ReleaseData();
  this->GetExecutive()->SetOutputData(1, pd);
  pd->Delete();

  this->Internals = new vtkChomboReader2Internals(this);
}

//----------------------------------------------------------------------------
vtkChomboReader2::~vtkChomboReader2()
{
  this->SetFileName(0);
  this->PointDataArraySelection->Delete();
  delete this->Internals;
}

//----------------------------------------------------------------------------
int vtkChomboReader2::RequestInformation(
  vtkInformation* request,
  vtkInformationVector** inputVector, 
  vtkInformationVector* outputVector)
{
  if (!this->Superclass::RequestInformation(request, inputVector, outputVector))
    {
    return 0;
    }
    
  vtkInformation* info = outputVector->GetInformationObject(0);
  vtkHierarchicalBoxDataSet *output = static_cast<vtkHierarchicalBoxDataSet *>(
    info->Get(vtkDataObject::DATA_OBJECT()));

  this->TotNumBoxes = 0;
  this->Internals->Boxes.clear();
  this->Internals->Offsets.clear();
  this->Internals->DXs.clear();

  hid_t fileHID = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if( fileHID < 0 )
    {
    vtkErrorMacro( << "Failed to open file " << this->FileName );
    return 0;
    }

  hid_t globalGroupId = H5Gopen(fileHID, "Chombo_global");
  if( globalGroupId < 0 )
    {
    vtkErrorMacro( << "Failed to open group Chombo_global");
    H5Fclose( fileHID );
    return 0;
    }

  if (!this->Internals->GetIntAttribute(globalGroupId, 
                                        "SpaceDim", 
                                        this->Dimensionality))
    {
    return 0;
    }

  hid_t attribute = H5Aopen_name(globalGroupId, "testReal");
  if ( attribute < 0 )
    {
    vtkErrorMacro( << "Failed to open attribute testReal");
    H5Gclose(globalGroupId);
    H5Fclose(fileHID);
    return 0;
    }

  hid_t datatype = H5Aget_type( attribute );
  if(H5Tget_precision(datatype) == H5Tget_precision(H5T_NATIVE_FLOAT))
    {
    this->RealType = vtkChomboReader2::FLOAT;
    } 
  else if(H5Tget_precision(datatype) == H5Tget_precision(H5T_NATIVE_DOUBLE))
    {
    this->RealType = vtkChomboReader2::DOUBLE;
    } 
  else
    {
    vtkErrorMacro( << "Illegal datatype. Only float and double are supported." );
    H5Aclose(attribute);
    H5Gclose(globalGroupId);
    H5Fclose(fileHID);
    return 0;
    }

  H5Aclose(attribute);

  H5Gclose(globalGroupId);

  hid_t rootGroupId = H5Gopen(fileHID, "/");
  if( rootGroupId < 0 )
    {
    vtkErrorMacro( << "Failed to open group /");
    H5Fclose( fileHID );
    return 0;
    }
  
  if (!this->Internals->GetIntAttribute( 
        rootGroupId, "num_components", this->NumberOfComponents))
    {
    H5Gclose(rootGroupId);
    H5Fclose(fileHID);
    return 0;
    }

  this->Internals->ComponentNames.resize(this->NumberOfComponents);
  for(int i=0; i<this->NumberOfComponents; i++)
    {
    ostrstream attrName;
    attrName << "component_" << i << ends;
    if (!this->Internals->GetStringAttribute( 
          rootGroupId, attrName.str(), this->Internals->ComponentNames[i] ) )
      {
      delete[] attrName.str();
      H5Gclose(rootGroupId);
      H5Fclose(fileHID);
      return 0;
      }
    this->PointDataArraySelection->AddArray((const char *)this->Internals->ComponentNames[i].c_str());

    delete[] attrName.str();
    }

  int numLevels;

  if (!this->Internals->GetIntAttribute( rootGroupId, "num_levels", numLevels))
    {
    H5Gclose(rootGroupId);
    H5Fclose(fileHID);
    return 0;
    }
  this->NumberOfLevels = numLevels;
  if (this->LevelRange[1] == -1)
    {
    this->LevelRange[1] = numLevels;
    }

  H5Gclose( rootGroupId );

  int lastLevel = this->LevelRange[1];
  if (lastLevel > this->NumberOfLevels)
    {
    lastLevel = this->NumberOfLevels;
    }
  int firstLevel = this->LevelRange[0];
  if (firstLevel < 0)
    {
    firstLevel = 0;
    }
  if (firstLevel > this->LevelRange[1])
    {
    firstLevel = this->LevelRange[1];
    }
  output->SetNumberOfLevels(this->NumberOfLevels);

  if ( numLevels > 0 )
    {
    typedef vtkstd::vector<vtkAMRBox> LevelBoxesType;
    typedef vtkstd::vector<LevelBoxesType> BoxesType;
    
    this->Internals->Boxes.resize(numLevels);
    this->Internals->DXs.resize(numLevels);
    this->Internals->Offsets.resize(numLevels);
    
    int levelId;
    for (levelId=0; levelId<numLevels; levelId++)
      {
      ostrstream levelName;
      levelName << "/level_" << levelId << ends;
      hid_t levelGroupId = H5Gopen( fileHID, levelName.str() );
      delete[] levelName.str();
      if( levelGroupId < 0 )
        {
        vtkErrorMacro( << "Failed to open group /level_" << levelId);
        H5Fclose(fileHID);
        return 0;
        }
      
      if (!this->Internals->GetRealTypeAttribute(
            levelGroupId, "dx", this->Internals->DXs[levelId]))
        {
        vtkErrorMacro( << "Failed to read dx for level " << levelId );
        H5Gclose(levelGroupId);
        H5Fclose(fileHID);
        return 0;
        }
      LevelBoxesType& boxes = this->Internals->Boxes[levelId];
      if (!this->Internals->ReadBoxes(levelGroupId, boxes))
        {
        vtkErrorMacro( << "Failed to read level " << levelId << " boxes");
        H5Gclose(levelGroupId);
        H5Fclose(fileHID);
        return 0;
        }
      
      int numBoxes = boxes.size();
      if ( levelId >= firstLevel && levelId < lastLevel )
        {
        output->SetNumberOfDataSets(levelId, numBoxes);
        }
      else
        {
        output->SetNumberOfDataSets(levelId, 0);
        }

      if ( numBoxes > 0 )
        {
        vtkChomboReader2Internals::LevelOffsetsType& offsets = 
          this->Internals->Offsets[levelId];
        offsets.resize(numBoxes);
        offsets[0] = 0;
        int boxId;
        for(boxId=1; boxId<numBoxes; boxId++)
          {
          vtkAMRBox& box = boxes[boxId-1];
          offsets[boxId] = 
            offsets[boxId-1] + this->NumberOfComponents*box.GetNumberOfCells();
          }
        for(boxId=0; boxId<numBoxes; boxId++)
          {
          if ( levelId >= firstLevel && levelId < lastLevel )
            {
            vtkAMRBox& box = boxes[boxId];
            output->SetDataSet(levelId, boxId, box, 0);
            this->TotNumBoxes++;
            }
          }
        }
      H5Gclose(levelGroupId);
      }
    }
  hid_t particlesGroupId = H5Gopen(fileHID, "particles");
  if( particlesGroupId < 0 )
    {
    vtkErrorMacro( << "Failed to open group particles");
    H5Fclose( fileHID );
    return 0;
    }
  if (!this->Internals->GetIntAttribute(particlesGroupId, 
                                        "num_particles", 
                                        this->NumberOfParticles))
    {
    return 0;
    }
  H5Gclose(particlesGroupId);
  H5Fclose(fileHID);
  return 1;
}

//----------------------------------------------------------------------------
int vtkChomboReader2::RequestData(
  vtkInformation *request, 
  vtkInformationVector **inputVector, 
  vtkInformationVector *outputVector)
{
  vtkInformation* info = outputVector->GetInformationObject(0);
  vtkHierarchicalBoxDataSet *output = vtkHierarchicalBoxDataSet::SafeDownCast(
    info->Get(vtkDataObject::DATA_OBJECT()));

  hid_t fileHID = H5Fopen(this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if( fileHID < 0 )
    {
    vtkErrorMacro( << "Failed to open file " << this->FileName );
    return 0;
    }

  int piece, numberOfPieces, ghostLevel;
  piece = info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numberOfPieces = 
    info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());
  ghostLevel = 
    info->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS());
  vtkDebugMacro(<< "piece: " << piece 
                  << " number of pieces: " << numberOfPieces
                  << " ghost levels: " << ghostLevel);

  int numBoxesPerPiece = this->TotNumBoxes / numberOfPieces;
  int beginBox = piece*numBoxesPerPiece;
  int endBox = (piece+1)*numBoxesPerPiece;
  if (piece == numberOfPieces - 1)
    {
    endBox = this->TotNumBoxes;
    }
  if (numBoxesPerPiece == 0)
    {
    if (piece == 0)
      {
      beginBox = 0;
      endBox = this->TotNumBoxes;
      }
    else
      {
      beginBox = 1;
      endBox = 0;
      }
    }

  int lastLevel = this->LevelRange[1]+1;
  if (lastLevel > this->NumberOfLevels)
    {
    lastLevel = this->NumberOfLevels;
    }
  int firstLevel = this->LevelRange[0];
  if (firstLevel < 0)
    {
    firstLevel = 0;
    }
  if (firstLevel > this->LevelRange[1])
    {
    firstLevel = this->LevelRange[1];
    }
  output->SetNumberOfLevels(this->NumberOfLevels);

  int globalBoxId = 0;

  for (int levelId=0; levelId<this->NumberOfLevels; levelId++)
    {
    double dx = this->Internals->DXs[levelId];
    vtkChomboReader2Internals::LevelBoxesType& boxes = this->Internals->Boxes[levelId];
    vtkChomboReader2Internals::LevelOffsetsType& offsets = this->Internals->Offsets[levelId];
    int numBoxes = boxes.size();

    ostrstream levelName;
    levelName << "/level_" << levelId << ends;
    hid_t levelGroupId = H5Gopen( fileHID, levelName.str() );
    delete[] levelName.str();
    if( levelGroupId < 0 )
      {
      vtkErrorMacro( << "Failed to open group /level_" << levelId);
      H5Fclose(fileHID);
      return 0;
      }
    
    int ref;
    if (!this->Internals->GetIntAttribute(levelGroupId, "ref_ratio", ref))
      {
      vtkErrorMacro( << "Failed to read ref_ratio for level " << levelId );
      H5Gclose(levelGroupId);
      H5Fclose(fileHID);
      return 0;
      }
    output->SetRefinementRatio(levelId, ref);
      
    hid_t dataset = H5Dopen( levelGroupId, "data:datatype=0" );
    if ( dataset < 0 )
      {
      vtkErrorMacro( << "Failed to open group /level_" << levelId);
      H5Gclose( levelGroupId );
      H5Fclose(fileHID);
      return 0;
      }

    if ( levelId >= firstLevel && levelId < lastLevel )
      {
      output->SetNumberOfDataSets(levelId, numBoxes);
      
      hid_t fileSpace = H5Dget_space( dataset );
      
      for (int boxId=0; boxId<numBoxes; boxId++, globalBoxId++)
        {
        vtkAMRBox& box = boxes[boxId];

        if (globalBoxId < beginBox || globalBoxId >= endBox)
          {
          output->SetDataSet(levelId, boxId, box, 0);
          continue;
          }
        
        int numPts[3];
        numPts[2] = 1;
        for(int i=0; i<this->Dimensionality; i++)
          {
          numPts[i] = box.HiCorner[i] - box.LoCorner[i] + 2;
          }
        vtkUniformGrid* ug = vtkUniformGrid::New();;
        ug->SetSpacing(dx, dx, dx);
        ug->SetWholeExtent(0, numPts[0]-1, 0, numPts[1]-1, 0, numPts[2]-1);
        ug->SetOrigin(box.LoCorner[0]*dx, 
                      box.LoCorner[1]*dx, 
                      box.LoCorner[2]*dx);
        ug->SetDimensions(numPts[0], numPts[1], numPts[2]);
        
        output->SetDataSet(levelId, boxId, box, ug);
        ug->Delete();
        
        hsize_t numItems = box.GetNumberOfCells();
        hid_t memSpace = H5Screate_simple( 1, &numItems, 0 );

        for(int component=0; component < this->NumberOfComponents; component++)
          {
          if(!this->GetPointArrayStatus((const char *)this->Internals->ComponentNames[component].c_str())) continue;

          hsize_t offset = offsets[boxId] + component * numItems;
          H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, &offset, 0, &numItems, 0);
        
          vtkDataArray* array;
          if ( this->RealType == vtkChomboReader2::FLOAT )
            {
            array = vtkFloatArray::New();
            array->SetNumberOfTuples(numItems);
            float* buf = reinterpret_cast<float*>(array->GetVoidPointer(0));
            if (H5Dread(dataset, 
                    H5T_NATIVE_FLOAT, 
                    memSpace, 
                    fileSpace, 
                    H5P_DEFAULT, 
                    buf) < 0)
              {
              vtkErrorMacro( << "Error reading FArray." );
              array->Delete();
              H5Sclose(fileSpace);
              H5Dclose(dataset);
              H5Sclose(memSpace);
              H5Gclose(levelGroupId);
              H5Fclose(fileHID);
              return 0;
              }
            }
          else
            {
            array = vtkDoubleArray::New();
            array->SetNumberOfTuples(numItems);
            double* buf = reinterpret_cast<double*>(array->GetVoidPointer(0));
            if (H5Dread(dataset, 
                  H5T_NATIVE_DOUBLE, 
                  memSpace, 
                  fileSpace, 
                  H5P_DEFAULT, 
                  buf ) < 0)
              {
              vtkErrorMacro( << "Error reading FArray." );
              array->Delete();
              H5Sclose(fileSpace);
              H5Dclose(dataset);
              H5Sclose(memSpace);
              H5Gclose(levelGroupId);
              H5Fclose(fileHID);
              return 0;
              }
            }
          array->SetName(this->Internals->ComponentNames[component].c_str());
          ug->GetCellData()->AddArray(array);
          array->Delete();
        }// loop over all components

        H5Sclose(memSpace);
        } // loop over all boxes
        H5Sclose(fileSpace);
      } // loop over all levels
    else
      {
      output->SetNumberOfDataSets(levelId, 0);
      }
      
    H5Dclose(dataset);
    H5Gclose(levelGroupId);
    }
  output->GenerateVisibilityArrays();
  
// to mimic the feature of ChomboVis, we have been asked to provide a slider in logarithmic scale
// to undersample the particles. If 1, undersample by a factor of 10, if 2, by a factor of 100, etc.

  info = outputVector->GetInformationObject(1);
  vtkPolyData *output2 = vtkPolyData::SafeDownCast(info->Get(vtkDataObject::DATA_OBJECT()));
  
  hid_t particlesGroupId = H5Gopen(fileHID, "particles");
  
  if( particlesGroupId < 0 )
    {
    vtkErrorMacro( << "Failed to open group particles");
    H5Fclose( fileHID );
    return 0;
    }
  if (!this->Internals->GetIntAttribute(particlesGroupId, 
                                        "num_particles", 
                                        this->NumberOfParticles))
    {
    return 0;
    }

  int mystride = (int)pow(10, this->UnderSampleScaleFactor);
  hsize_t NumberOfParticlesToRead = this->NumberOfParticles / mystride;

  NumberOfParticlesToRead /= numberOfPieces;

  // to save memory, we make all particles single-length "float" length. See data type below.
  // if double, HDF5 does on-the-fly conversion to float
  vtkPoints *newPoints = vtkPoints::New();
  newPoints->SetDataTypeToFloat();
  newPoints->Allocate(NumberOfParticlesToRead);
  vtkCellArray *newVerts = vtkCellArray::New();
  newVerts->Allocate(newVerts->EstimateSize(1, NumberOfParticlesToRead));

  hid_t posxDataSetId = H5Dopen(particlesGroupId, "position_x");
  hid_t posyDataSetId = H5Dopen(particlesGroupId, "position_y");
  hid_t poszDataSetId = H5Dopen(particlesGroupId, "position_z");
  
  hid_t fileSpaceX = H5Dget_space(posxDataSetId);
  hid_t fileSpaceY = H5Dget_space(posyDataSetId);
  hid_t fileSpaceZ = H5Dget_space(poszDataSetId);

  hid_t memSpace = H5Screate_simple( 1, &NumberOfParticlesToRead, 0 );
  hsize_t stride[1], start[1], count[1], block[1];
    
  start[0]  = (this->NumberOfParticles / numberOfPieces ) * piece;
  stride[0] = mystride;
  count[0]  = NumberOfParticlesToRead;
  block[0]  = 1;
    
  if ( this->RealType == vtkChomboReader2::DOUBLE || this->RealType == vtkChomboReader2::FLOAT)
    {
    float* bufx = new float[NumberOfParticlesToRead];
    float* bufy = new float[NumberOfParticlesToRead];
    float* bufz = new float[NumberOfParticlesToRead];

    H5Sselect_hyperslab(fileSpaceX, H5S_SELECT_SET, start, stride, count, block);
    H5Sselect_hyperslab(fileSpaceY, H5S_SELECT_SET, start, stride, count, block);
    H5Sselect_hyperslab(fileSpaceZ, H5S_SELECT_SET, start, stride, count, block);
    if (
       (H5Dread(posxDataSetId, H5T_NATIVE_FLOAT, memSpace, fileSpaceX, H5P_DEFAULT, bufx) < 0) || 
       (H5Dread(posyDataSetId, H5T_NATIVE_FLOAT, memSpace, fileSpaceY, H5P_DEFAULT, bufy) < 0) ||
       (H5Dread(poszDataSetId, H5T_NATIVE_FLOAT, memSpace, fileSpaceZ, H5P_DEFAULT, bufz) < 0))	    
      {
            vtkErrorMacro( << "Error reading position_?." );

            H5Dclose(posxDataSetId);
            H5Dclose(posyDataSetId);
            H5Dclose(poszDataSetId);
            H5Fclose(fileHID);
            return 0;
      }

    for (int i=0; i < NumberOfParticlesToRead; i++)
      {
      newVerts->InsertNextCell(1);
      newVerts->InsertCellPoint(newPoints->InsertNextPoint(bufx[i], bufy[i], bufz[i]));
      }
    delete [] bufx;
    delete [] bufy;           
    delete [] bufz;
    }
  H5Dclose(posxDataSetId); //close the datasets
  H5Dclose(posyDataSetId);
  H5Dclose(poszDataSetId);
  H5Sclose(fileSpaceX);    // close the filespaces
  H5Sclose(fileSpaceY);
  H5Sclose(fileSpaceZ);
//
// now load the scalars.
  hid_t velxDataSetId = H5Dopen(particlesGroupId, "velocity_x");
  hid_t velyDataSetId = H5Dopen(particlesGroupId, "velocity_y");
  hid_t velzDataSetId = H5Dopen(particlesGroupId, "velocity_z");
  hid_t massDataSetId = H5Dopen(particlesGroupId, "mass");
  
  vtkDataArray* array;
  if ( this->RealType == vtkChomboReader2::DOUBLE || this->RealType == vtkChomboReader2::FLOAT )
    {
    array = vtkFloatArray::New();
    array->SetNumberOfTuples(NumberOfParticlesToRead);
    float* buf = reinterpret_cast<float*>(array->GetVoidPointer(0));
    hid_t fileSpace = H5Dget_space(velxDataSetId);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start, stride, count, block);
    
    if (H5Dread(velxDataSetId,  H5T_NATIVE_FLOAT,  memSpace, fileSpace,  H5P_DEFAULT,  buf) < 0)
      {
      vtkErrorMacro( << "Error reading velocity_x" );
      array->Delete();
      H5Dclose(velxDataSetId);
      H5Sclose(fileSpace);
      H5Fclose(fileHID);
      return 0;
      }
    H5Sclose(fileSpace);
    H5Dclose(velxDataSetId);
    }

  array->SetName("velocity_x");
  output2->GetPointData()->AddArray(array);
  array->Delete();
  
  if ( this->RealType == vtkChomboReader2::DOUBLE || this->RealType == vtkChomboReader2::FLOAT )
    {
    array = vtkFloatArray::New();
    array->SetNumberOfTuples(NumberOfParticlesToRead);
    float* buf = reinterpret_cast<float*>(array->GetVoidPointer(0));
    hid_t fileSpace = H5Dget_space(velyDataSetId);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start, stride, count, block);
    
    if (H5Dread(velyDataSetId,  H5T_NATIVE_FLOAT,  memSpace, fileSpace,  H5P_DEFAULT,  buf) < 0)
      {
      vtkErrorMacro( << "Error reading velocity_y" );
      array->Delete();
      H5Dclose(velyDataSetId);
      H5Sclose(fileSpace);
      H5Fclose(fileHID);
      return 0;
      }
    H5Sclose(fileSpace);
    H5Dclose(velyDataSetId);
    }

  array->SetName("velocity_y");
  output2->GetPointData()->AddArray(array);
  array->Delete();

  if ( this->RealType == vtkChomboReader2::DOUBLE || this->RealType == vtkChomboReader2::FLOAT )
    {
    array = vtkFloatArray::New();
    array->SetNumberOfTuples(NumberOfParticlesToRead);
    float* buf = reinterpret_cast<float*>(array->GetVoidPointer(0));
    hid_t fileSpace = H5Dget_space(velzDataSetId);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start, stride, count, block);
    
    if (H5Dread(velzDataSetId,  H5T_NATIVE_FLOAT,  memSpace, fileSpace,  H5P_DEFAULT,  buf) < 0)
      {
      vtkErrorMacro( << "Error reading velocity_z" );
      array->Delete();
      H5Dclose(velzDataSetId);
      H5Sclose(fileSpace);
      H5Fclose(fileHID);
      return 0;
      }
    H5Sclose(fileSpace);
    H5Dclose(velzDataSetId);
    }

  array->SetName("velocity_z");
  output2->GetPointData()->AddArray(array);
  array->Delete();

  if ( this->RealType == vtkChomboReader2::DOUBLE || this->RealType == vtkChomboReader2::FLOAT )
    {
    array = vtkFloatArray::New();
    array->SetNumberOfTuples(NumberOfParticlesToRead);
    float* buf = reinterpret_cast<float*>(array->GetVoidPointer(0));
    hid_t fileSpace = H5Dget_space(massDataSetId);
    H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, start, stride, count, block);
    
    if (H5Dread(massDataSetId,  H5T_NATIVE_FLOAT,  memSpace, fileSpace,  H5P_DEFAULT,  buf) < 0)
      {
      vtkErrorMacro( << "Error reading mass" );
      array->Delete();
      H5Dclose(massDataSetId);
      H5Sclose(fileSpace);
      H5Fclose(fileHID);
      return 0;
      }
    H5Sclose(fileSpace);
    H5Dclose(massDataSetId);
    }

  array->SetName("mass");
  output2->GetPointData()->AddArray(array);
  array->Delete();

  output2->SetPoints(newPoints);
  newPoints->Delete();

  output2->SetVerts(newVerts);
  newVerts->Delete();

  H5Gclose(particlesGroupId);

  H5Fclose(fileHID);

  return 1;
}

//----------------------------------------------------------------------------
int vtkChomboReader2::FillOutputPortInformation(
  int port, vtkInformation* info)
{
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkHierarchicalBoxDataSet");
  if(port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPolyData");

  return 1;
}

//----------------------------------------------------------------------------
void vtkChomboReader2::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: "
     << (this->FileName? this->FileName:"(none)") << "\n";
  os << indent << "RealType: " << this->RealType << endl;
  os << indent << "NumberOfComponents: " << this->NumberOfComponents << endl;
  os << indent << "Dimensionality: " << this->Dimensionality << endl;
  os << indent << "LevelRange: " 
     << this->LevelRange[0] << " " << this->LevelRange[1] << endl;
  os << indent << "NumberOfLevels: " << this->NumberOfLevels << endl;
}

const char* vtkChomboReader2::GetPointArrayName(int index)
{
  return this->PointDataArraySelection->GetArrayName(index);
}

int vtkChomboReader2::GetPointArrayStatus(const char* name)
{
  return this->PointDataArraySelection->ArrayIsEnabled(name);
}

void vtkChomboReader2::SetPointArrayStatus(const char* name, int status)
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

int vtkChomboReader2::GetNumberOfPointArrays()
{
  return this->PointDataArraySelection->GetNumberOfArrays();
}

void vtkChomboReader2::Enable(const char* name)
{
  this->SetPointArrayStatus(name, 1);
}

void vtkChomboReader2::Disable(const char* name)
{
  this->SetPointArrayStatus(name, 0);
}

void vtkChomboReader2::EnableAll()
{
  this->PointDataArraySelection->EnableAllArrays();
}

void vtkChomboReader2::DisableAll()
{
  this->PointDataArraySelection->DisableAllArrays();
}
