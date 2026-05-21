// A simple Validity check for the vtkAMROverlapping datastructure
// with a second argument (filename), it can produce the same data in Kitware's own format
// Can be used also with valgrind and other diagnostics

#include "vtkAMRAmazeReader.h"
#include "vtkOverlappingAMR.h"

#include "vtkNew.h"

#include "vtkXMLUniformGridAMRWriter.h"

#include <iostream>

int main(int argc, char **argv)
{
  vtkNew<vtkAMRAmazeReader> reader;
  if(argc < 2)
    {
    std::cerr << "missing a filename argument: Syntax ./bin/TestAMRValidity <path-to>/file.amr5 [optional output.vtkhb]\n";
    exit(1);
    }
  std::cout << "Opening " <<  argv[1] << std::endl;
  reader->SetFileName(argv[1]);
  reader->DebugOn();
  reader->DataScaleOn();
  reader->SetScaleChoice(0);
  reader->LogDataOff();
  reader->UpdateInformation();
  reader->DisableAll();
  reader->Enable("Density");
  // if we do not specify levels, should read them all
  //reader->SetLevelRead(0, 19);
  reader->DebugOff();
  reader->Update();
  
  vtkOverlappingAMR* amr = nullptr;
  amr = reader->GetOutput();
  if (amr != nullptr)
  {
    if (!amr->CheckValidity())
    {
      std::cerr << "ERROR: output AMR dataset is not valid!\n";
      return 1;
    }
  }
  else
  {
    std::cerr << "ERROR: output AMR dataset is nullptr!";
    return 1;
  }

  if(argc == 3)
  {
  vtkNew<vtkXMLUniformGridAMRWriter> writer;
  writer->SetFileName(argv[2]);
  writer->SetInputData(amr);
  writer->Write();
  }
  if (!amr->CheckValidity())
  {
    std::cerr << "amr->CheckValidity() failure\n";
    return -1;
  }
  else
    return 0;
}
