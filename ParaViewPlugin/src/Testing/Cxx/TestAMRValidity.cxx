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
  std::cout << "TestAMRValidity opening " <<  argv[1] << std::endl;
  reader->SetFileName(argv[1]);
  reader->DebugOn();
  reader->DataScaleOn();
  reader->LogDataOff();
  reader->SetPointArrayStatus("Density", 1);

  reader->DebugOn();
  reader->UpdateInformation();
  reader->SetMaxLevel(reader->GetNumberOfLevels());
  reader->Update();
  vtkOverlappingAMR* amr = nullptr;
  amr = reader->GetOutput();
  if (amr != nullptr)
  {
    if (!amr->CheckValidity())
    {
      std::cerr << "ERROR: output AMR dataset is not valid!\n";
      //return 1;
    }
    else
      {
      std::cerr << "this AMR dataset passed the validity check\n";
      if(argc == 3)
        {
        vtkNew<vtkXMLUniformGridAMRWriter> writer;
        writer->SetFileName(argv[2]);
        std::cout << "TestAMRValidity writing " <<  argv[2] << std::endl;
        writer->SetInputData(amr);
        writer->Write();
        }
     } 
  }
  else
  {
    std::cerr << "ERROR: output AMR dataset is nullptr!";
    //return 1;
  }


}
