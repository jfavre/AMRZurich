#ifndef __vtkAMRChomboReader_h
#define __vtkAMRChomboReader_h

#include "vtkCompositeDataSetAlgorithm.h"
#include "vtkDataArraySelection.h"

#include "vtkIOAMRModule.h" // For export macro
#include "vtkAMRBaseReader.h"

class vtkAMRChomboReaderInternal;
class vtkOverlappingAMR;

class VTK_EXPORT vtkAMRChomboReader : public vtkAMRBaseReader
{
public:
  vtkTypeMacro(vtkAMRChomboReader,vtkAMRBaseReader);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkAMRChomboReader* New();

  // Description:
  // Get/Set the name of the input file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  // Description:
  vtkGetMacro(UnderSampleScaleFactor, double);
  vtkSetMacro(UnderSampleScaleFactor, double);

  // Description:
  vtkGetMacro(NumberOfComponents, int);

  // Description:
  vtkGetMacro(NumberOfLevels, int);

  // Description:
  vtkGetMacro(Dimensionality, int);

  // Description:
  vtkGetMacro(RealType, int);

  // Description:
  vtkSetVector2Macro(LevelRange, int);
  vtkGetVector2Macro(LevelRange, int);

  vtkGetObjectMacro(PointDataArraySelection, vtkDataArraySelection);
  int GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int GetPointArrayStatus(const char* name);
  void SetPointArrayStatus(const char* name, int status);
  void DisableAll();  
  void EnableAll();
  void Disable(const char* name);  
  void Enable(const char* name);

protected:
  vtkAMRChomboReader();
  ~vtkAMRChomboReader();

  // Description:
  // See vtkAMRBaseReader::ReadMetaData
  void ReadMetaData();

  // Description:
  // See vtkAMRBaseReader::GetBlockLevel
  int GetBlockLevel( const int blockIdx );
  // Description:
  // See vtkAMRBaseReader::FillMetaData
  int FillMetaData( );

  // Description:
  // See vtkAMRBaseReader::GetAMRGrid
  vtkUniformGrid* GetAMRGrid( const int blockIdx );

  // Description:
  // See vtkAMRBaseReader::GetAMRGridData
  void GetAMRGridData(
      const int blockIdx, vtkUniformGrid *block, const char *field);

  // Description:
  // See vtkAMRBaseReader::GetAMRGridData
  void GetAMRGridPointData(
      const int vtkNotUsed(blockIdx), vtkUniformGrid *vtkNotUsed(block), const char *vtkNotUsed(field)) {;};

  // Description:
  // See vtkAMRBaseReader::SetUpDataArraySelections
  void SetUpDataArraySelections();

  vtkAMRChomboReaderInternal *Internal;

  enum Real_T
  {
    FLOAT,
    DOUBLE
  };

  int LevelRange[2];
  int RealType;
  int Dimensionality;
  int NumberOfLevels;
  int NumberOfComponents;
  int NumberOfParticles;
  int TotNumBoxes;
  double UnderSampleScaleFactor; // take the log of this to undersmaple particlesam
  vtkAMRChomboReaderInternals* Internals;

private:
  vtkAMRChomboReader(const vtkAMRChomboReader&);  // Not implemented.
  void operator=(const vtkAMRChomboReader&);  // Not implemented.
};

#endif
