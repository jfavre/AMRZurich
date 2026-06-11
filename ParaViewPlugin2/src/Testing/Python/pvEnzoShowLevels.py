# state file generated using paraview version 6.1.1-617-g6391ede1d0
import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.ViewSize = [1295, 658]

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.Set(
    ViewSize=[971, 1071],
    CenterOfRotation=[0.7509765625, 0.75048828125, 0.75048828125],
    CameraPosition=[0.5224628769618425, 0.75048828125, 0.75048828125],
    CameraFocalPoint=[0.5377894756118452, 0.75048828125, 0.75048828125],
    CameraViewUp=[0.0, 0.0, 1.0],
)

# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.Set(
    ViewSize=[970, 1071],
    CenterOfRotation=[0.7509765625, 0.75048828125, 0.75048828125],
    CameraPosition=[0.5224628769618425, 0.75048828125, 0.75048828125],
    CameraFocalPoint=[0.5377894756118452, 0.75048828125, 0.75048828125],
    CameraViewUp=[0.0, 0.0, 1.0],
)

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.500000)
layout1.AssignView(1, renderView2)
layout1.AssignView(2, renderView3)
layout1.SetSize(1942, 1071)

# create new layout object 'Layout #1'
layout1_1 = CreateLayout(name='Layout #1')
layout1_1.AssignView(0, renderView1)
layout1_1.SetSize(1295, 658)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView3)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Enzo Reader'
reader = EnzoReader(registrationName='moving7_0010.hierarchy', FileNames=['/local/data/Walder/Enzo/DD0010/moving7_0010.hierarchy'])
reader.Set(
    Level=7,
    CellArrayStatus=['Density'],
)

# create a new 'Axis-Aligned Slice'
axisAlignedSlice1 = AxisAlignedSlice(registrationName='AxisAlignedSlice1', Input=reader)
axisAlignedSlice1.Level = 7

# init the 'Axis Aligned Plane' selected for 'CutFunction'
axisAlignedSlice1.CutFunction.Origin = [0.7498756825440642, 0.5, 0.5]

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=reader)
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Origin = [0.7498756825440642, 0.5, 0.5]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from axisAlignedSlice1
axisAlignedSlice1Display = Show(axisAlignedSlice1, renderView2, 'AMRRepresentation')

# get color transfer function/color map for 'vtkAMRLevel'
vtkAMRLevelLUT = GetColorTransferFunction('vtkAMRLevel')
vtkAMRLevelLUT.Set(
    InterpretValuesAsCategories=1,
    AnnotationsInitialized=1,
    ScalarRangeInitialized=1.0,
    Annotations=['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7'],
    ActiveAnnotatedValues=['0', '1', '2', '3', '4', '5', '6', '7'],
    IndexedColors=[1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.6299992203712463, 0.6299992203712463, 1.0],
    IndexedOpacities=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
)

# trace defaults for the display properties.
axisAlignedSlice1Display.Set(
    Representation='Surface With Edges',
    ColorArrayName=['CELLS', 'vtkAMRLevel'],
    LookupTable=vtkAMRLevelLUT,
    Assembly='Hierarchy',
)

# setup the color legend parameters for each legend in this view

# get color legend/bar for vtkAMRLevelLUT in view renderView2
vtkAMRLevelLUTColorBar = GetScalarBar(vtkAMRLevelLUT, renderView2)
vtkAMRLevelLUTColorBar.Set(
    WindowLocation='Upper Right Corner',
    Title='vtkAMRLevel',
    ComponentTitle='',
)

# set color bar visibility
vtkAMRLevelLUTColorBar.Visibility = 1

# show color legend
axisAlignedSlice1Display.SetScalarBarVisibility(renderView2, True)

slice1Display = Show(slice1, renderView3, 'GeometryRepresentation')

# get color transfer function/color map for 'Density'
densityLUT = GetColorTransferFunction('Density')
densityLUT.Set(
    RGBPoints=[
        # scalar, red, green, blue
        0.09997177124023439, 0.0564, 0.0564, 0.47,
        0.9256784561960354, 0.243, 0.46035, 0.81,
        4.800494160958115, 0.356814, 0.745025, 0.954368,
        27.16894158082622, 0.6882, 0.93, 0.91791,
        65.52274055026149, 0.899496, 0.944646, 0.768657,
        205.7650925873603, 0.957108, 0.833819, 0.508916,
        949.7085632135038, 0.927521, 0.621439, 0.315357,
        5951.937631029598, 0.8, 0.352, 0.16,
        42944.41796875, 0.59, 0.0767, 0.119475,
    ],
    UseLogScale=1,
    ScalarRangeInitialized=1.0,
    EnableOpacityMapping=1,
)

# trace defaults for the display properties.
slice1Display.Set(
    Representation='Surface With Edges',
    ColorArrayName=['CELLS', 'Density'],
    LookupTable=densityLUT,
    Assembly='Hierarchy',
)

# setup the color legend parameters for each legend in this view

# get color legend/bar for densityLUT in view renderView3
densityLUTColorBar = GetScalarBar(densityLUT, renderView3)
densityLUTColorBar.Set(
    Title='Density',
    ComponentTitle='',
)

# set color bar visibility
densityLUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView3, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'vtkAMRLevel'
vtkAMRLevelPWF = GetOpacityTransferFunction('vtkAMRLevel')

# get opacity transfer function/opacity map for 'Density'
densityPWF = GetOpacityTransferFunction('Density')
densityPWF.Set(
    Points=[0.09997177124023438, 0.0, 0.5, 0.0, 50.0, 0.0, 0.5, 0.0, 60.0, 1.0, 0.5, 0.0, 42944.41796875, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)


AddCameraLink(renderView2, renderView3, 'CameraLink0')

# save screenshot
SaveScreenshot(filename='/local/data/Walder/foo.png', viewOrLayout=layout1, location=16, SaveAllViews=1)

