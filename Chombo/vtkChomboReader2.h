/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile: vtkChomboReader2.h,v $
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
// .NAME vtkChomboReader2 -
// .SECTION Description

#ifndef __vtkChomboReader2_h
#define __vtkChomboReader2_h

#include "vtkCompositeDataSetAlgorithm.h"
#include "vtkDataArraySelection.h"
class vtkHierarchicalBoxDataSet;
//BTX
struct vtkChomboReader2Internals;
//ETX

class VTK_EXPORT vtkChomboReader2 : public vtkCompositeDataSetAlgorithm
{
public:
  vtkTypeRevisionMacro(vtkChomboReader2,vtkCompositeDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);
  static vtkChomboReader2* New();

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
  vtkChomboReader2();
  ~vtkChomboReader2();
  
  // Standard pipeline exectution methods.
  virtual int RequestInformation(vtkInformation *request, 
                                 vtkInformationVector **inputVector, 
                                 vtkInformationVector *outputVector);
  virtual int RequestData(vtkInformation *request, 
                          vtkInformationVector **inputVector, 
                          vtkInformationVector *outputVector);
  
  // The input file's name.
  char* FileName;

  // see algorithm for more info
  virtual int FillOutputPortInformation(int port, vtkInformation* info);

//BTX
  friend struct vtkChomboReader2Internals;

  enum Real_T
  {
    FLOAT,
    DOUBLE
  };

//ETX

  int LevelRange[2];
  int RealType;
  int Dimensionality;
  int NumberOfLevels;
  int NumberOfComponents;
  int NumberOfParticles;
  int TotNumBoxes;
  double UnderSampleScaleFactor; // take the log of this to undersmaple particlesam
  vtkChomboReader2Internals* Internals;

  vtkDataArraySelection *PointDataArraySelection;

 
private:
  vtkChomboReader2(const vtkChomboReader2&);  // Not implemented.
  void operator=(const vtkChomboReader2&);  // Not implemented.
};

#endif
