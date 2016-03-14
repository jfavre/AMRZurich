// .NAME vtkAMRAmazeReader - Reads A-MAZE AMR files (in development)
// .SECTION Description

#ifndef __vtkAMRAmazeReader_h
#define __vtkAMRAmazeReader_h

class vtkMultiBlockDataSet;
class vtkPolyData;
class vtkAMRBox;
class vtkUniformGrid;
class vtkDoubleArray;

#include "vtkIOAMRModule.h" // For export macro
#include "vtkAMRBaseReader.h"
#include "vtkAMRAmazeReaderInternal.h"
#include <vector> // Needed for vector ivar
#include <map>

class vtkOverlappingAMR;

class VTKIOAMR_EXPORT vtkAMRAmazeReader : public vtkAMRBaseReader
{
public:
  static vtkAMRAmazeReader* New();
  vtkTypeMacro(vtkAMRAmazeReader, vtkAMRBaseReader);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // See vtkAMRBaseReader::SetFileName
  void SetFileName( const char* fileName );

  vtkSetStringMacro(HDF5SaveFileName);
  vtkGetStringMacro(HDF5SaveFileName);

  int GetNumberOfComponents(){return myreader->NumberOfComponents;};
  int GetNumberOfLevels(){return myreader->NumberOfLevels;};

  int GetNumberOfBlocks(){return myreader->NumberOfGrids;};

  double GetTime(){return myreader->Time;};
  void SetTime(double T) {myreader->Time = T;};

  vtkSetMacro(MaxLevelWrite, int);
  vtkGetMacro(MaxLevelWrite, int);

  // The range of valid levels values.
  int LevelRange[2];
  vtkGetVector2Macro(LevelRange, int);

  int LevelRead[2];
  vtkSetVector2Macro(LevelRead, int);
  vtkGetVector2Macro(LevelRead, int);	

  vtkSetMacro(ShiftedGrid, int);
  vtkGetMacro(ShiftedGrid, int);
  vtkBooleanMacro(ShiftedGrid, int);

  void SetLogData(int _arg);
  vtkGetMacro(LogData, int);
  vtkBooleanMacro(LogData, int);

  void SetScaleChoice(int _arg);
  vtkGetMacro(ScaleChoice, int);

  void SetLengthScale(int _arg);
  vtkGetMacro(LengthScale, int);
  vtkBooleanMacro(LengthScale, int);

  vtkSetMacro(DataScale, int);
  vtkGetMacro(DataScale, int);
  vtkBooleanMacro(DataScale, int);

  //vtkGetObjectMacro(PointDataArraySelection, vtkDataArraySelection);
/*
  int GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int GetPointArrayStatus(const char* name);
  void SetPointArrayStatus(const char* name, int status);
  void DisableAll();  
  void EnableAll();
  void Disable(const char* name);  
  void Enable(const char* name);
*/ 
  vtkMultiBlockDataSet* GetStarsOutput();
  int LoadStars(hid_t root_id, vtkMultiBlockDataSet*);
  int GridsPerLevels(int l){return myreader->Levels[l].GridsPerLevel;};
  vtkSetMacro(MaximumLevelsToReadByDefault, unsigned int);
  vtkGetMacro(MaximumLevelsToReadByDefault, unsigned int);
  int CanReadFile(const char* fname);

protected:
  vtkAMRAmazeReader();
  ~vtkAMRAmazeReader();

   virtual int RequestInformation(
      vtkInformation* rqst,
      vtkInformationVector** inputVector,
      vtkInformationVector* outputVector );
  virtual int FillOutputPortInformation(int port, vtkInformation* info);

  // The input file's name.
  hid_t file_id;
  //char* FileName;
  char* HDF5SaveFileName;
  int LogData; // will automatically calculate log10() for Density, Temperature and Pressure
  int LengthScale; // will automatically scale the grids to real length
  int ShiftedGrid; //will use the shifter grid to provide stationary slab animation
  int ScaleChoice;
  int DataScale;
  double TimeScalor;
  int MaxLevelWrite;
  int MaxLevelRead;
  int MinLevelRead;
  unsigned int MaximumLevelsToReadByDefault;

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

  // This is for cell-data. Not used in AMAZE
  void GetAMRGridData(
      const int vtkNotUsed(blockIdx), vtkUniformGrid *vtkNotUsed(block), const char *vtkNotUsed(field)) {;};

  // Description:
  // See vtkAMRBaseReader::GetAMRGridData. This is for node-data used in AMAZE
  void GetAMRGridPointData(
      const int blockIdx, vtkUniformGrid *block, const char *field);

  // Description:
  // See vtkAMRBaseReader::SetUpDataArraySelections
  void SetUpDataArraySelections();

 private:
  vtkAMRAmazeReader(const vtkAMRAmazeReader&);  // Not implemented.
  void operator=(const vtkAMRAmazeReader&);  // Not implemented.

  vtkAMRAmazeReaderInternal *myreader;
};

#endif
