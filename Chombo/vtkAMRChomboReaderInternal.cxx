#include "vtkAMRChomboReaderInternal.h"

#include <hdf5.h>

struct vtkAMRChomboReaderInternal
{
  vtkAMRChomboReaderInternal();
  ~vtkAMRChomboReaderInternal();

  int GetIntAttribute( hid_t locID, const char* attrName, int& val );
  int GetRealTypeAttribute( hid_t locID, const char* attrName, double& val );
  int GetStringAttribute( hid_t locID, const char* attrName, vtkstd::string& val );
  int ReadBoxes( hid_t locID , vtkstd::vector<vtkAMRBox>& boxes );
  void CreateBoxDataType();

  hid_t BoxDataType;
  vtkAMRChomboReader* Reader;
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
vtkAMRChomboReaderInternal::vtkAMRChomboReaderInternal()
{
  this->BoxDataType = -1;
  this->Reader = reader;
}

//----------------------------------------------------------------------------
vtkAMRChomboReaderInternal::~vtkAMRChomboReaderInternal()
{
  if ( this->BoxDataType >= 0 )
    {
    H5Tclose( this->BoxDataType );
    }
}

//----------------------------------------------------------------------------
void vtkAMRChomboReaderInternal::CreateBoxDataType()
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
int vtkAMRChomboReaderInternal::GetRealTypeAttribute( 
  hid_t locID, const char* attrName, double& val )
{
  hid_t attribute = H5Aopen_name( locID, attrName );
  if ( attribute < 0 )
    {
    vtkGenericWarningMacro( << "Failed to open attribute " << attrName );
    return 0;
    }

  if ( this->Reader->RealType == vtkAMRChomboReader::FLOAT )
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
int vtkAMRChomboReaderInternal::GetIntAttribute( 
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
int vtkAMRChomboReaderInternal::GetStringAttribute( 
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
int vtkAMRChomboReaderInternal::ReadBoxes( 
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

