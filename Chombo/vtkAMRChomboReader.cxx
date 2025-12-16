#include "vtkAMRChomboReader.h"
#include "vtkObjectFactory.h"
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

#include "vtkAMRUtilities.h"
#define PARALLEL_DEBUG 1

#include <vector>
#include <string>

#include <hdf5.h>

vtkStandardNewMacro(vtkAMRChomboReader);

//----------------------------------------------------------------------------
vtkAMRChomboReader::vtkAMRChomboReader()
{
  this->FileName = 0;

  this->Dimensionality = 0;
  this->NumberOfLevels = 0;
  this->NumberOfComponents = 0;

  this->LevelRange[0] =  0;
  this->LevelRange[1] = -1;
  this->UnderSampleScaleFactor = 2;
  this->TotNumBoxes = 0;

  this->Internal = new vtkAMRChomboReaderInternal();
}

//----------------------------------------------------------------------------
vtkAMRChomboReader::~vtkAMRChomboReader()
{
  delete this->Internal;
  this->Internal=NULL;

  delete [] this->FileName;
  this->FileName = NULL;
}

//----------------------------------------------------------------------------
int vtkAMRChomboReader::RequestInformation(
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
  this->Internal->Boxes.clear();
  this->Internal->Offsets.clear();
  this->Internal->DXs.clear();

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

  if (!this->Internal->GetIntAttribute(globalGroupId, 
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
    this->RealType = vtkAMRChomboReader::FLOAT;
    } 
  else if(H5Tget_precision(datatype) == H5Tget_precision(H5T_NATIVE_DOUBLE))
    {
    this->RealType = vtkAMRChomboReader::DOUBLE;
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
  
  if (!this->Internal->GetIntAttribute( 
        rootGroupId, "num_components", this->NumberOfComponents))
    {
    H5Gclose(rootGroupId);
    H5Fclose(fileHID);
    return 0;
    }

  this->Internal->ComponentNames.resize(this->NumberOfComponents);
  for(int i=0; i<this->NumberOfComponents; i++)
    {
    ostrstream attrName;
    attrName << "component_" << i << ends;
    if (!this->Internal->GetStringAttribute( 
          rootGroupId, attrName.str(), this->Internal->ComponentNames[i] ) )
      {
      delete[] attrName.str();
      H5Gclose(rootGroupId);
      H5Fclose(fileHID);
      return 0;
      }
    this->PointDataArraySelection->AddArray((const char *)this->Internal->ComponentNames[i].c_str());

    delete[] attrName.str();
    }

  int numLevels;

  if (!this->Internal->GetIntAttribute( rootGroupId, "num_levels", numLevels))
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
    
    this->Internal->Boxes.resize(numLevels);
    this->Internal->DXs.resize(numLevels);
    this->Internal->Offsets.resize(numLevels);
    
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
      
      if (!this->Internal->GetRealTypeAttribute(
            levelGroupId, "dx", this->Internal->DXs[levelId]))
        {
        vtkErrorMacro( << "Failed to read dx for level " << levelId );
        H5Gclose(levelGroupId);
        H5Fclose(fileHID);
        return 0;
        }
      LevelBoxesType& boxes = this->Internal->Boxes[levelId];
      if (!this->Internal->ReadBoxes(levelGroupId, boxes))
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
        vtkAMRChomboReaderInternal::LevelOffsetsType& offsets = 
          this->Internal->Offsets[levelId];
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
  if (!this->Internal->GetIntAttribute(particlesGroupId, 
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
int vtkAMRChomboReader::RequestData(
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
    double dx = this->Internal->DXs[levelId];
    vtkAMRChomboReaderInternal::LevelBoxesType& boxes = this->Internal->Boxes[levelId];
    vtkAMRChomboReaderInternal::LevelOffsetsType& offsets = this->Internal->Offsets[levelId];
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
    if (!this->Internal->GetIntAttribute(levelGroupId, "ref_ratio", ref))
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
          if(!this->GetPointArrayStatus((const char *)this->Internal->ComponentNames[component].c_str())) continue;

          hsize_t offset = offsets[boxId] + component * numItems;
          H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, &offset, 0, &numItems, 0);
        
          vtkDataArray* array;
          if ( this->RealType == vtkAMRChomboReader::FLOAT )
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
          array->SetName(this->Internal->ComponentNames[component].c_str());
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

  H5Fclose(fileHID);

  return 1;
}


//----------------------------------------------------------------------------
void vtkAMRChomboReader::PrintSelf(ostream& os, vtkIndent indent)
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


