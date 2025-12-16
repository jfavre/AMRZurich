#ifndef __vtkAMRChomboReaderInternal_h
#define __vtkAMRChomboReaderInternal_h

class vtkDataArray;
class vtkDataSet;

class vtkAMRChomboReaderInternal
{
public:
  vtkAMRChomboReaderInternal();
  ~vtkAMRChomboReaderInternal();

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

  char* FileName;

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

#endif
