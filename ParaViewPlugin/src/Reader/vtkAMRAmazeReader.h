// .NAME vtkAMRAmazeReader - Reads AMAZE AMR files
// .SECTION Description
// This is an AMAZE AMR file reader.

#ifndef __vtkAMRAmazeReader_h
#define __vtkAMRAmazeReader_h

class vtkDataArraySelection;
class vtkMultiBlockDataSet;
class vtkPolyData;
class vtkAMRBox;
class vtkUniformGrid;
class vtkDoubleArray;

#include "vtkOverlappingAMRAlgorithm.h"
#include "AMAZEReaderModule.h" // for export macro
#include "vtkAMRAmazeReaderInternal.h"

#include <vector> // Needed for vector ivar
#include <map>

class AMAZEREADER_EXPORT vtkAMRAmazeReader : public vtkOverlappingAMRAlgorithm
{
public:
  vtkTypeMacro(vtkAMRAmazeReader,vtkOverlappingAMRAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkAMRAmazeReader* New();

  // Description:
  // Get/Set the name of the input file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  int GetNumberOfComponents(){return Internal->NumberOfComponents;};
  int GetNumberOfLevels(){return Internal->NumberOfLevels;};
  int GetNumberOfBlocks(){return Internal->NumberOfGrids;};

  double GetTime(){return Internal->GetAMAZETime();};
  void SetTime(double T) {Internal->SetAMAZETime(T);};

  vtkSetMacro(MaxLevelWrite, int);
  vtkGetMacro(MaxLevelWrite, int);

  // The range of valid levels values.
  int LevelRange[2];
  vtkGetVector2Macro(LevelRange, int);

  int LevelRead[2];
  vtkSetVector2Macro(LevelRead, int);
  vtkGetVector2Macro(LevelRead, int);	

  vtkSetMacro(ShiftedGrid, vtkTypeBool);
  vtkGetMacro(ShiftedGrid, vtkTypeBool);
  vtkBooleanMacro(ShiftedGrid, vtkTypeBool);

  vtkSetMacro(LogData, vtkTypeBool);
  vtkGetMacro(LogData, vtkTypeBool);
  vtkBooleanMacro(LogData, vtkTypeBool);

  vtkSetMacro(ScaleChoice, int);
  vtkGetMacro(ScaleChoice, int);

  vtkSetMacro(LengthScale, vtkTypeBool);
  vtkGetMacro(LengthScale, vtkTypeBool);
  vtkBooleanMacro(LengthScale, vtkTypeBool);

  vtkSetMacro(LengthScaleFactor, double);
  vtkGetMacro(LengthScaleFactor, double);

  vtkSetMacro(DataScale, vtkTypeBool);
  vtkGetMacro(DataScale, vtkTypeBool);
  vtkBooleanMacro(DataScale, vtkTypeBool)

  vtkGetObjectMacro(PointDataArraySelection, vtkDataArraySelection);

  int GetNumberOfPointArrays();
  const char* GetPointArrayName(int index);
  int GetPointArrayStatus(const char* name);
  void SetPointArrayStatus(const char* name, int status);

  void DisableAll();  
  void EnableAll();
  void Disable(const char* name);  
  void Enable(const char* name); 
  vtkMultiBlockDataSet* GetStarsOutput();
  int LoadStars(vtkMultiBlockDataSet*);
  int GridsPerLevels(int l){return Internal->Levels[l].GridsPerLevel;};
  vtkSetMacro(MaximumLevelsToReadByDefault, unsigned int);
  vtkGetMacro(MaximumLevelsToReadByDefault, unsigned int);
  vtkTypeBool CanReadFile(VTK_FILEPATH const char* fname);

protected:
  vtkAMRAmazeReader();
  ~vtkAMRAmazeReader();

  // New pipeline execution methods.
  virtual int RequestData(vtkInformation*, 
                  vtkInformationVector**, 
                  vtkInformationVector*);
  virtual int RequestInformation(vtkInformation*, 
                         vtkInformationVector**, 
                         vtkInformationVector*);
 
  virtual int FillOutputPortInformation(int port, vtkInformation* info);

  /**
   * See vtkAMRBaseReader::GetAMRGrid
   */
  vtkUniformGrid* GetAMRGrid(int blockIdx);

  /**
   * See vtkAMRBaseReader::GetAMRGridData   for CellData NOT USED here
   */
  void GetAMRGridData(int vtkNotUsed(blockIdx), vtkUniformGrid* vtkNotUsed(block),
    const char* vtkNotUsed(field)) {};

  /**
   * See vtkAMRBaseReader::GetAMRGridPointData
   */
  void GetAMRGridPointData(int blockIdx, vtkUniformGrid* block, const char* field);
  
  int GetBlockLevel(int blockIdx);
  void SetUpDataArraySelections();
  void FillMetaData(vtkOverlappingAMR*);
    
  // The input file's name.
  char* FileName;

  vtkTypeBool LogData; // will automatically calculate log10() for Density, Temperature and Pressure
  vtkTypeBool LengthScale; // will automatically scale the grids to real length
  vtkTypeBool ShiftedGrid; //will use the shifter grid to provide stationary slab animation
  int ScaleChoice;
  vtkTypeBool DataScale;
  double LengthScaleFactor;

  int MaxLevelWrite;
  int MaxLevelRead;
  int MinLevelRead;
  unsigned int MaximumLevelsToReadByDefault;
  vtkAMRAmazeReaderInternal *Internal;

  int nbstars;
  vtkDataArraySelection* PointDataArraySelection;

 private:
  vtkAMRAmazeReader(const vtkAMRAmazeReader&) = delete;  // Not implemented.
  void operator=(const vtkAMRAmazeReader&) = delete;  // Not implemented.
};

#endif
