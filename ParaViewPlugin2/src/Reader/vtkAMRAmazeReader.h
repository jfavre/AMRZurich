/**
 * @class   vtkAMRAmazeReader
 * @brief   A concrete instance of vtkAMRBaseReader that implements functionality
 * for reading Amaze AMR datasets.
 */

#ifndef vtkAMRAmazeReader_h
#define vtkAMRAmazeReader_h

#include "vtkAMRBaseReader.h"
#include "vtkIOAMRModule.h" // For export macro
#include <memory> // For std::unique_ptr

class vtkOverlappingAMR;
class vtkAMRAmazeReaderInternal;
class vtkMultiBlockDataSet;

class VTKIOAMR_EXPORT vtkAMRAmazeReader : public vtkAMRBaseReader
{
public:
  static vtkAMRAmazeReader* New();
  vtkTypeMacro(vtkAMRAmazeReader, vtkAMRBaseReader);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //--- Grid Shift Settings ---
  vtkSetMacro(ShiftedGrid, vtkTypeBool);
  vtkGetMacro(ShiftedGrid, vtkTypeBool);
  vtkBooleanMacro(ShiftedGrid, vtkTypeBool);

  //--- Data Logarithmic Scaling Settings ---
  vtkSetMacro(LogData, vtkTypeBool);
  vtkGetMacro(LogData, vtkTypeBool);
  vtkBooleanMacro(LogData, vtkTypeBool);

  //--- Scale Choice Settings ---
  vtkSetMacro(ScaleChoice, int);
  vtkGetMacro(ScaleChoice, int);

  //--- Length Scaling Settings ---
  vtkSetMacro(LengthScale, vtkTypeBool);
  vtkGetMacro(LengthScale, vtkTypeBool);
  vtkBooleanMacro(LengthScale, vtkTypeBool);

  //--- Data Scaling Settings ---
  vtkSetMacro(DataScale, vtkTypeBool);
  vtkGetMacro(DataScale, vtkTypeBool);
  vtkBooleanMacro(DataScale, vtkTypeBool);

  //--- Overridden Base Reader API ---
  int GetNumberOfBlocks() override;
  int GetNumberOfLevels() override;
  void SetFileName(VTK_FILEPATH const char* fileName) override;

  //--- Amaze Specific API ---
  vtkMultiBlockDataSet* GetStarsOutput();
  int LoadStars(vtkMultiBlockDataSet* starsDataSet);
  vtkTypeBool CanReadFile(VTK_FILEPATH const char* fname) const;

protected:
  vtkAMRAmazeReader();
  ~vtkAMRAmazeReader() override;

  //--- Pipeline Methods ---
  int RequestInformation(vtkInformation* request,
                         vtkInformationVector** inputVector,
                         vtkInformationVector* outputVector) override;
                  
  int RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector) override;
                  
  int FillOutputPortInformation(int port, vtkInformation* info) override;
  
  //--- Metadata and Grid Parsing ---
  void ReadMetaData() override;
  int GetBlockLevel(int blockIdx) override;
  int FillMetaData() override;

  vtkUniformGrid* GetAMRGrid(int blockIdx) override;
  void GetAMRGridData(int vtkNotUsed(blockIdx), vtkUniformGrid* vtkNotUsed(block), const char* vtkNotUsed(field)) override {}
  void GetAMRGridPointData(int blockIdx, vtkUniformGrid* block, const char* field) override;
  void SetUpDataArraySelections() override;

  //--- Internal State Variables ---
  vtkTypeBool LogData{0};          // Calculates log10() for Density, Temperature and Pressure
  vtkTypeBool LengthScale;      // Scales the grids to real length
  vtkTypeBool ShiftedGrid{0};      // Uses shifted grid to provide stationary slab animation
  vtkTypeBool DataScale;
  int ScaleChoice;

  bool IsReady;
  unsigned int MaximumLevelsToReadByDefault;
  int nbstars;

private:
  vtkAMRAmazeReader(const vtkAMRAmazeReader&) = delete;
  void operator=(const vtkAMRAmazeReader&) = delete;

  std::unique_ptr<vtkAMRAmazeReaderInternal> Internal;
};

#endif /* vtkAMRAmazeReader_h */
