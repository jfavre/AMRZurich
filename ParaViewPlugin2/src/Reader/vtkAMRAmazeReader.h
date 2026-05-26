
/**
 * @class   vtkAMRAmazeReader
 * @brief   A concrete instance of vtkAMRBaseReader that implements functionality
 * for reading Amaze AMR datasets.
 */

#ifndef vtkAMRAmazeReader_h
#define vtkAMRAmazeReader_h

#include "vtkAMRBaseReader.h"
#include "vtkIOAMRModule.h" // For export macro

class vtkOverlappingAMR;
class vtkAMRAmazeReaderInternal;
class vtkMultiBlockDataSet;

class VTKIOAMR_EXPORT vtkAMRAmazeReader : public vtkAMRBaseReader
{
public:
  static vtkAMRAmazeReader* New();
  vtkTypeMacro(vtkAMRAmazeReader, vtkAMRBaseReader);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  vtkSetMacro(ShiftedGrid, vtkTypeBool);
  vtkGetMacro(ShiftedGrid, vtkTypeBool);
  vtkBooleanMacro(ShiftedGrid, vtkTypeBool);

  vtkSetMacro(LogData, vtkTypeBool);
  vtkGetMacro(LogData, vtkTypeBool);
  vtkBooleanMacro(LogData, vtkTypeBool);

  vtkSetMacro(ScaleChoice, int);
  vtkGetMacro(ScaleChoice, int);

  vtkSetMacro(LengthScale, vtkTypeBool);
  vtkGetMacro(LengthScale, vtkTypeBool);
  vtkBooleanMacro(LengthScale, vtkTypeBool);

  vtkSetMacro(LengthScaleFactor, double);
  vtkGetMacro(LengthScaleFactor, double);

  vtkSetMacro(DataScale, vtkTypeBool);
  vtkGetMacro(DataScale, vtkTypeBool);
  vtkBooleanMacro(DataScale, vtkTypeBool)

  /**
   * See vtkAMRBaseReader::GetNumberOfBlocks
   */
  int GetNumberOfBlocks() override;

  /**
   * See vtkAMRBaseReader::GetNumberOfLevels
   */
  int GetNumberOfLevels() override;

  /**
   * See vtkAMRBaseReader::SetFileName
   */
  void SetFileName(VTK_FILEPATH const char* fileName) override;
  vtkMultiBlockDataSet* GetStarsOutput();
  int LoadStars(vtkMultiBlockDataSet*);
  vtkTypeBool CanReadFile(VTK_FILEPATH const char* fname);
  
  vtkTypeBool LogData; // will automatically calculate log10() for Density, Temperature and Pressure
  vtkTypeBool LengthScale; // will automatically scale the grids to real length
  vtkTypeBool ShiftedGrid; //will use the shifter grid to provide stationary slab animation
  int ScaleChoice;
  vtkTypeBool DataScale;
  double LengthScaleFactor;
  
protected:
  vtkAMRAmazeReader();
  ~vtkAMRAmazeReader() override;

  int RequestInformation(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector) override;
                  
  int RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector) override;
                  
  int FillOutputPortInformation(int port, vtkInformation* info) override;
  
  
  /**
   * See vtkAMRBaseReader::ReadMetaData
   */
  void ReadMetaData() override;

  /**
   * See vtkAMRBaseReader::GetBlockLevel
   */
  int GetBlockLevel(int blockIdx) override;

  /**
   * See vtkAMRBaseReader::FillMetaData
   */
  int FillMetaData() override;

  /**
   * See vtkAMRBaseReader::GetAMRGrid
   */
  vtkUniformGrid* GetAMRGrid(int blockIdx) override;

  /**
   * See vtkAMRBaseReader::GetAMRGridData
   */
    void GetAMRGridData(int vtkNotUsed(blockIdx), vtkUniformGrid* vtkNotUsed(block),const char* vtkNotUsed(field)) override {};
 
  /**
   * See vtkAMRBaseReader::GetAMRGridData
   */
  void GetAMRGridPointData(int blockIdx, vtkUniformGrid* block, const char* field) override;
  /**
   * See vtkAMRBaseReader::SetUpDataArraySelections
   */
  void SetUpDataArraySelections() override;

  bool IsReady;

  unsigned int MaximumLevelsToReadByDefault;

  int nbstars;

private:
  vtkAMRAmazeReader(const vtkAMRAmazeReader&) = delete;
  void operator=(const vtkAMRAmazeReader&) = delete;

  vtkAMRAmazeReaderInternal* Internal;
};

#endif /* vtkAMRAmazeReader_h */
