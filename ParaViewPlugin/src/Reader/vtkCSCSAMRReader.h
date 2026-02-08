// .NAME vtkCSCSAMRReader - Reads CSCSAMR files (in development)
// .SECTION Description
// This is an experimental CSCSAMR file reader. It is mainly used
// for development and does not support all features of CSCSAMR
// format. Use at your own risk.

#ifndef __vtkCSCSAMRReader_h
#define __vtkCSCSAMRReader_h

class vtkDataArraySelection;
class vtkMultiBlockDataSet;
class vtkPolyData;
class vtkAMRBox;
class vtkUniformGrid;
class vtkDoubleArray;

#include "vtkOverlappingAMRAlgorithm.h"
#include "AMAZEReaderModule.h" // for export macro
#include "vtkAMAZEReader.h"

#include <vector> // Needed for vector ivar
#include <map>

class AMAZEREADER_EXPORT vtkCSCSAMRReader : public vtkOverlappingAMRAlgorithm
{
public:
  vtkTypeMacro(vtkCSCSAMRReader,vtkOverlappingAMRAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkCSCSAMRReader* New();

  // Description:
  // Get/Set the name of the input file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);
  vtkSetStringMacro(HDF5SaveFileName);
  vtkGetStringMacro(HDF5SaveFileName);

  int GetNumberOfComponents(){return myreader->NumberOfComponents;};
  int GetNumberOfLevels(){return myreader->NumberOfLevels;};
  int GetNumberOfGrids(){return myreader->NumberOfGrids;};

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

  vtkSetMacro(LogData, int);
  vtkGetMacro(LogData, int);
  vtkBooleanMacro(LogData, int);

  vtkSetMacro(ScaleChoice, int);
  vtkGetMacro(ScaleChoice, int);

  vtkSetMacro(LengthScale, int);
  vtkGetMacro(LengthScale, int);
  vtkBooleanMacro(LengthScale, int);

  vtkSetMacro(LengthScaleFactor, double);
  vtkGetMacro(LengthScaleFactor, double);

  vtkSetMacro(DataScale, int);
  vtkGetMacro(DataScale, int);
  vtkBooleanMacro(DataScale, int)

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
  int LoadStars(hid_t root_id, vtkMultiBlockDataSet*);
  int GridsPerLevels(int l){return myreader->Levels[l].GridsPerLevel;};
  vtkSetMacro(MaximumLevelsToReadByDefault, unsigned int);
  vtkGetMacro(MaximumLevelsToReadByDefault, unsigned int);
  int CanReadFile(const char* fname);

protected:
  vtkCSCSAMRReader();
  ~vtkCSCSAMRReader();

  // New pipeline execution methods.
  virtual int RequestData(vtkInformation*, 
                  vtkInformationVector**, 
                  vtkInformationVector*);
  virtual int RequestInformation(vtkInformation*, 
                         vtkInformationVector**, 
                         vtkInformationVector*);
 
  virtual int FillOutputPortInformation(int port, vtkInformation* info);

  // The input file's name.
  hid_t file_id;
  char* FileName;
  char* HDF5SaveFileName;
  int LogData; // will automatically calculate log10() for Density, Temperature and Pressure
  int LengthScale; // will automatically scale the grids to real length
  int ShiftedGrid; //will use the shifter grid to provide stationary slab animation
  int ScaleChoice;
  int DataScale;
  double LengthScaleFactor;
  double TimeScalor;
  int MaxLevelWrite;
  int MaxLevelRead;
  int MinLevelRead;
  unsigned int MaximumLevelsToReadByDefault;
  vtkAMAZEReader *myreader;

  vtkDataArraySelection* PointDataArraySelection;

 private:
  vtkCSCSAMRReader(const vtkCSCSAMRReader&) = delete;  // Not implemented.
  void operator=(const vtkCSCSAMRReader&) = delete;  // Not implemented.
};

#endif
