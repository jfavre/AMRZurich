#include "vtkAMAZEFileSeriesReader.h"

#include "vtkObjectFactory.h"

vtkStandardNewMacro(vtkAMAZEFileSeriesReader);

vtkAMAZEFileSeriesReader::vtkAMAZEFileSeriesReader()
{
  // Override base class's SetNumberOfOutputPorts(1) to support 2 outputs
  this->SetNumberOfOutputPorts(2);
}

void vtkAMAZEFileSeriesReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}
