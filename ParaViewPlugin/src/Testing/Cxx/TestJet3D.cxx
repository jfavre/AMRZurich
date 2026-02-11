#include <iostream>
#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
#include "vtkAMRAmazeReader.h"
#include "vtkOverlappingAMR.h"
#include "vtkArrayCalculator.h"
#include "vtkCamera.h"
#include "vtkClipDataSet.h"
#include "vtkCompositeDataPipeline.h"
#include "vtkCompositeDataSet.h"
#include "vtkContourFilter.h"
#include "vtkCutter.h"
#include "vtkDataSet.h"
#include "vtkExecutive.h"
#include "vtkCompositeDataGeometryFilter.h"
#include "vtkCompositePolyDataMapper.h"
#include "vtkLookupTable.h"
#include "vtkOutlineCornerFilter.h"
#include "vtkOutlineFilter.h"
#include "vtkPlane.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkTIFFWriter.h"
#include "vtkTimerLog.h"
#include "vtkWindowToImageFilter.h"
#include "vtkTestUtilities.h"
#include "vtkXMLHierarchicalBoxDataReader.h"
#include "vtkXMLHierarchicalBoxDataWriter.h"

#include <iostream>
using namespace std;

#define VTK_CREATE(type, var) \
  vtkSmartPointer<type> var = vtkSmartPointer<type>::New();

//#define ALL 1

int main(int argc, char **argv)
{
  //vtkCompositeDataPipeline* prototype = vtkCompositeDataPipeline::New();
  //vtkAlgorithm::SetDefaultExecutivePrototype(prototype);
  //prototype->Delete();

  VTK_CREATE(vtkTimerLog, timer);
  timer->StartTimer();

  VTK_CREATE(vtkAMRAmazeReader, reader);
  reader->SetFileName("/local/data/Walder/jet3d.amr5");
  reader->DebugOff();
  reader->DataScaleOn();
  reader->LogDataOn();
  reader->UpdateInformation();
  reader->DisableAll();
  reader->Enable("Density");
  //reader->SetPhysicalSpaceScale(1.0);
  //reader->SetPointArrayStatus("Magnetic field", 0);
  reader->SetLevelRead(0, 4);
  reader->Update();
  //static_cast<vtkOverlappingAMR*>(reader->GetOutput())->Audit();
  //reader->GetOutput()->Print(cout);

#ifdef ALL
  VTK_CREATE(vtkContourFilter, contour);
  contour->SetInputConnection(0, reader->GetOutputPort(0));
  contour->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "Log10(Density) [N/cm^3]");
  //contour->GenerateValues(16, 6876, 6876349);
  contour->GenerateValues(9, 5.46, 7.63);
  contour->ComputeScalarsOn();
  contour->DebugOff();

  VTK_CREATE(vtkLookupTable, lut);
  lut->SetTableRange(5.46, 7.63); // 3.83, 6.83);
  lut->SetHueRange(0.66,0.0);
  lut->SetNumberOfTableValues(256);
  lut->Build();

  VTK_CREATE(vtkPlane, plane);
  plane->SetNormal(0,-1,0);
  plane->SetOrigin(.25,.25,.16037);

  VTK_CREATE(vtkClipDataSet, clip);
  clip->SetInputConnection(0, contour->GetOutputPort(0));
  clip->SetClipFunction(plane);

  VTK_CREATE(vtkCompositeDataGeometryFilter, clipGeom);
  clipGeom->SetInputConnection(0, clip->GetOutputPort(0)); // or contour

  // Rendering objects
  VTK_CREATE(vtkPolyDataMapper, clipMapper);
  clipMapper->SetInputConnection(0, clipGeom->GetOutputPort(0));
  clipMapper->SetLookupTable(lut);
  clipMapper->UseLookupTableScalarRangeOn();
  clipMapper->SetColorModeToMapScalars();
  clipMapper->ScalarVisibilityOn();
  clipMapper->SetScalarModeToUsePointFieldData();
  clipMapper->SelectColorArray(0);

  VTK_CREATE(vtkActor, clipActor);
  clipActor->SetMapper(clipMapper);
  //clipActor->GetProperty()->SetOpacity(.3);
  VTK_CREATE(vtkCompositePolyDataMapper, clipMapper2);
  clipMapper2->SetInputConnection(0, clipGeom->GetOutputPort(0));
  clipMapper2->ScalarVisibilityOff();

  VTK_CREATE(vtkActor, clipActor2);
  clipActor2->SetMapper(clipMapper2);

  VTK_CREATE(vtkProperty, backP);
  backP->SetDiffuseColor(.6, .6, .6);
  clipActor2->SetProperty(backP);
  clipActor2->GetProperty()->BackfaceCullingOn();
  //clipActor2->SetBackfaceProperty(backP);
  //clipActor2->GetProperty()->FrontfaceCullingOn();

  VTK_CREATE(vtkCutter, cutter);
  cutter->SetInputConnection(0, contour->GetOutputPort(0));
  cutter->SetCutFunction(plane);

  VTK_CREATE(vtkCompositeDataGeometryFilter, cutterGeom);
  cutterGeom->SetInputConnection(0, cutter->GetOutputPort(0));

  VTK_CREATE(vtkCompositePolyDataMapper, cutterMapper);
  cutterMapper->SetInputConnection(0, cutter->GetOutputPort(0));
  //cutterMapper->SetResolveCoincidentTopologyToPolygonOffset();
  cutterMapper->SetResolveCoincidentTopologyToShiftZBuffer();
  cutterMapper->SetLookupTable(lut);
  cutterMapper->UseLookupTableScalarRangeOn();
  cutterMapper->SetScalarModeToUsePointData();
  cutterMapper->SelectColorArray("Log10(Density)");
  //cutterMapper->ScalarVisibilityOff();

  VTK_CREATE(vtkActor, cutterActor);
  cutterActor->SetMapper(cutterMapper);
  //cutterActor->GetProperty()->SetColor(.3,.3,.3);

  // corner outline
  VTK_CREATE(vtkOutlineCornerFilter, ocf);
  ocf->SetInputConnection(0, reader->GetOutputPort(0));

  // Rendering objects
  // This one is actually just a vtkPolyData so it doesn't need a hierarchical
  // mapper, but we use this one to test hierarchical mapper with polydata input
  VTK_CREATE(vtkCompositePolyDataMapper, ocMapper);
  ocMapper->SetInputConnection(0, ocf->GetOutputPort(0));

  VTK_CREATE(vtkActor, ocActor);
  ocActor->SetMapper(ocMapper);
  ocActor->GetProperty()->SetColor(0, 0, 0);

  VTK_CREATE(vtkRenderer, Ren1);
  Ren1->SetBackground(0.33, 0.35, 0.43);
  Ren1->GetActiveCamera()->ParallelProjectionOn();
  Ren1->GetActiveCamera()->SetParallelScale(0.049615);
  Ren1->GetActiveCamera()->SetFocalPoint(0.257375, 0.274539, 0.107061);
  Ren1->GetActiveCamera()->SetPosition(0.347757, 0.547007, 0.306317);
  Ren1->GetActiveCamera()->SetViewUp(0.964893, -0.23644, -0.11436);

  //Ren1->SetUseDepthPeeling(1);
  //Ren1->SetMaximumNumberOfPeels(4);
  //Ren1->SetOcclusionRatio(0.0);

  VTK_CREATE(vtkRenderWindow, RenWin1);
  RenWin1->AddRenderer(Ren1);
  RenWin1->SetSize(1024,640);
  RenWin1->SetMultiSamples(1);
  RenWin1->SetAlphaBitPlanes(1);

  VTK_CREATE(vtkRenderWindowInteractor, iren);
  iren->SetRenderWindow(RenWin1);

  Ren1->AddActor(clipActor);
  Ren1->AddActor(clipActor2);
  Ren1->AddActor(cutterActor);
  Ren1->AddActor(ocActor);
  RenWin1->Render();
/*
  cerr  << *reader;
  VTK_CREATE(vtkWindowToImageFilter, wToImg);
  wToImg->SetInput(RenWin1);

  VTK_CREATE(vtkTIFFWriter, writer);
  writer->SetInputConnection(0,wToImg->GetOutputPort(0));
  writer->SetFileName("/tmp/jet3d.tif");
  writer->Write();
*/
  iren->Initialize();
  iren->Start();

  timer->StopTimer();
  cerr << "time elapsed " << timer->GetElapsedTime() << endl;

  if(Ren1->GetLastRenderingUsedDepthPeeling())
    {
    cout<<"depth peeling was used"<<endl;
    }
  else
    {
    cout<<"depth peeling was not used (alpha blending instead)"<<endl;
    }
#endif

  return 0;
}
