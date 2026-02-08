// .NAME vtkCSCSAMRReader - Reads CSCSAMR files (in development)
// .SECTION Description
// This is an experimental CSCSAMR file reader. It is mainly used
// for development and does not support all features of CSCSAMR
// format. Use at your own risk.

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
