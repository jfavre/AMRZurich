// .NAME vtkAMAZEFileSeriesReader - Reads series of AMAZE HDF5 files
// .SECTION Description

#ifndef vtkAMAZEFileSeriesReader_h
#define vtkAMAZEFileSeriesReader_h

#include "vtkFileSeriesReader.h"
#include "AMAZEReaderModule.h" // needed for export

class AMAZEREADER_EXPORT vtkAMAZEFileSeriesReader : public vtkFileSeriesReader
{
public:
  static vtkAMAZEFileSeriesReader* New();
  vtkTypeMacro(vtkAMAZEFileSeriesReader, vtkFileSeriesReader);
  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  vtkAMAZEFileSeriesReader();
  ~vtkAMAZEFileSeriesReader() override = default;

private:
  vtkAMAZEFileSeriesReader(const vtkAMAZEFileSeriesReader&) = delete;
  void operator=(const vtkAMAZEFileSeriesReader&) = delete;
};

#endif
