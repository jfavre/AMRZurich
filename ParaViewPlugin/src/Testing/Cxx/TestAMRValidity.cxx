// A simple Validity check for the vtkAMROverlapping datastructure
// with a second argument (filename), it can produce the same data in Kitware's own format
// Can be used also with valgrind and other diagnostics

#include "vtkAMRAmazeReader.h"
#include "vtkOverlappingAMR.h"

#include "vtkSmartPointer.h"

#include "vtkXMLUniformGridAMRWriter.h"

#include <iostream>
using namespace std;

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New();

int main(int argc, char **argv)
{
  VTK_CREATE(vtkAMRAmazeReader, reader);
  if(argc < 2)
    {
    std::cerr << "missing a filename argument: Syntax ./bin/TestAMRValidity <path-to>/file.amr5 [optional output.vtkhb]\n";
    exit(1);
    }
  std::cout << "Opening " <<  argv[1] << std::endl;
  reader->SetFileName(argv[1]);
  reader->DebugOn();
  reader->DataScaleOn();
  reader->LogDataOff();
  reader->UpdateInformation();
  reader->DisableAll();
  reader->Enable("Density");
  // if we do not specify levels, should read them all
  //reader->SetLevelRead(0, 19);
  reader->DebugOff();
  reader->Update();
  
  vtkSmartPointer<vtkOverlappingAMR> amr;
  amr = vtkOverlappingAMR::SafeDownCast(reader->GetOutputDataObject(0));
  if (!amr->CheckValidity())
  {
    std::cerr << "amr->CheckValidity() failure\n";
    return -1;
  }
  else if(argc == 3)
  {
  vtkNew<vtkXMLUniformGridAMRWriter> writer;
  writer->SetFileName(argv[2]);
  writer->SetInputData(amr);
  writer->Write();
  }

  return 0;
}
