#include <iostream>
#include "vtkAMRAmazeReaderInternal.h"
#include "vtkAMRAmazeReader.h"

#include "vtkRectilinearGrid.h"
#include "vtkSmartPointer.h"

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New(

int main(int argc, char **argv)
{

  if(argc != 2)
    {
    cerr << "Usage: TestAMAZEReader   *amr5\n"; exit(1);
    }


  vtkAMRAmazeReaderInternal *reader = vtkAMRAmazeReaderInternal::New();
  reader->SetFileName(argv[1]);
  reader->ReadMetaData();
  cerr << "reader->GetScaleChoice() = " << reader->GetScaleChoice() << endl;
  //reader->WriteXMLFile("/tmp/amr.vthb");

  vtkRectilinearGrid *rg;
  int dims[3];
/*
  for(int i=0; i < reader->GetNumberOfGrids(); i++)
    {
    rg = reader->ReadRectilinearGrid(i);
    rg->GetDimensions(dims);
    cerr << "dims[] = " << dims[0] << ", " << dims[1] << ", " << dims[2] << endl;
    }
*/
  reader->Delete();
/*
  vtkAMRAmazeReader *Reader = vtkAMRAmazeReader::New();
  Reader->SetFileName(argv[1]);
  Reader->CanReadFile(argv[1]);
  Reader->DataScaleOn();
  Reader->SetLengthScale(1);
  Reader->SetScaleChoice(3);
  Reader->LogDataOn();
  //Reader->SetPointArrayStatus("Density", 1);
  Reader->SetMaxLevel(2);
  //Reader->UpdateInformation();
  Reader->Update();
  Reader->Delete();
*/
  return 0;
}
