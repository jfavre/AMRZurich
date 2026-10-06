# state file generated using paraview version 6.2.0
import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.Set(
    ViewSize=[635, 566],
    CenterOfRotation=[48575107039232.0, 50906936639488.0, 50000000188416.0],
    CameraPosition=[48573864061551.234, 50905597061954.17, 49999009834585.65],
    CameraFocalPoint=[48574512006126.11, 50906295362502.305, 49999526090349.22],
    CameraViewUp=[0.2845263357202557, 0.38505479403026266, -0.8779393884978145],
)

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.Set(
    ViewSize=[635, 566],
    CenterOfRotation=[48575107039232.0, 50906936639488.0, 50000000188416.0],
    CameraPosition=[48573864061551.234, 50905597061954.17, 49999009834585.65],
    CameraFocalPoint=[48574512006126.11, 50906295362502.305, 49999526090349.22],
    CameraViewUp=[0.2845263357202557, 0.38505479403026266, -0.8779393884978145],
)

# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.Set(
    ViewSize=[633, 566],
    CenterOfRotation=[48575107039232.0, 50906936639488.0, 50000000188416.0],
    CameraPosition=[48573864061551.234, 50905597061954.17, 49999009834585.65],
    CameraFocalPoint=[48574512006126.11, 50906295362502.305, 49999526090349.22],
    CameraViewUp=[0.2845263357202557, 0.38505479403026266, -0.8779393884978145],
)

# Create a new 'Render View'
renderView4 = CreateView('RenderView')
renderView4.Set(
    ViewSize=[633, 566],
    CenterOfRotation=[48575107039232.0, 50906936639488.0, 50000000188416.0],
    CameraPosition=[48573864061551.234, 50905597061954.17, 49999009834585.65],
    CameraFocalPoint=[48574512006126.11, 50906295362502.305, 49999526090349.22],
    CameraViewUp=[0.2845263357202557, 0.38505479403026266, -0.8779393884978145],
)

# Create a new 'Render View'
renderView5 = CreateView('RenderView')
renderView5.Set(
    ViewSize=[633, 566],
    CenterOfRotation=[48575107039232.0, 50906936639488.0, 50000000188416.0],
    CameraPosition=[48573864061551.234, 50905597061954.17, 49999009834585.65],
    CameraFocalPoint=[48574512006126.11, 50906295362502.305, 49999526090349.22],
    CameraViewUp=[0.2845263357202557, 0.38505479403026266, -0.8779393884978145],
)

# Create a new 'Render View'
renderView6 = CreateView('RenderView')
renderView6.Set(
    ViewSize=[633, 566],
    CenterOfRotation=[48575107039232.0, 50906936639488.0, 50000000188416.0],
    CameraPosition=[48573864061551.234, 50905597061954.17, 49999009834585.65],
    CameraFocalPoint=[48574512006126.11, 50906295362502.305, 49999526090349.22],
    CameraViewUp=[0.2845263357202557, 0.38505479403026266, -0.8779393884978145],
)

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.333333)
layout1.SplitVertical(1, 0.500000)
layout1.AssignView(3, renderView1)
layout1.AssignView(4, renderView2)
layout1.SplitHorizontal(2, 0.500000)
layout1.SplitVertical(5, 0.500000)
layout1.AssignView(11, renderView4)
layout1.AssignView(12, renderView3)
layout1.SplitVertical(6, 0.500000)
layout1.AssignView(13, renderView5)
layout1.AssignView(14, renderView6)
layout1.SetSize(1903, 1133)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView6)
# ----------------------------------------------------------------

directory_JeanLaptop = '/local/data/Walder/CygX-1_M1=14.8_M2=19.2_V=750_ML=1-6_Gamma=1.01_RBH=20_HR/DATA/'
directory_r740gpu02 = '/Xnfs/highenergy/rwalder/Simulations/CygX-1_M1=14.8_M2=19.2_V=750_ML=1-6_Gamma=1.01_RBH=20_HR/'

fname = 'DATA.CX-1_M1=14.8_M2=19.2_V=750_ML=1-6_G=1.01_RBH=20RG_.NT=0000008600.Time=189.4754123165H.amr5'

# create a new directory name for Rolf and Doris's laptops, or use Jean's laptop, or use r740gpu02
# watch for the last "/" in the directory name

#filename = directory_r740gpu02 + fname
filename = directory_JeanLaptop + fname

reader = AMAZEReader(registrationName='reader', FileNames=filename)
reader.Set(
    PointArrayStatus=['Density', 'Energy Density', 'Thermal Energy Density', 'Velocity'],
    Level=100,
)
reader.UpdatePipeline()

# create a new 'Extract Block'
extractBlock1 = ExtractBlock(registrationName='ExtractBlock1', Input=reader)
extractBlock1.Set(
    Assembly='Hierarchy',
    Selectors=['/Root/Level17', '/Root/Level18'],
)

# temporar solution to get a simple vtkImageData (a single grid) to test with

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=extractBlock1)
resampleToImage1.SamplingDimensions = [512, 512, 256]

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(registrationName='ProgrammableFilter1', Input=resampleToImage1)
programmableFilter1.Set(
    RequestInformationScript='',
    RequestUpdateExtentScript='',
    CopyArrays=1,
    PythonPath='',
    Script="""import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

image = inputs[0]

nx, ny, nz = image.GetDimensions()
dx, dy, dz = image.GetSpacing()
ox, oy, oz = image.GetOrigin()

# Node coordinates
x = ox + np.arange(nx) * dx
y = oy + np.arange(ny) * dy
z = oz + np.arange(nz) * dz

# 1-D trapezoidal weights
wx = np.ones(nx)
wy = np.ones(ny)
wz = np.ones(nz)

if nx > 1:
    wx[[0, -1]] = 0.5
if ny > 1:
    wy[[0, -1]] = 0.5
if nz > 1:
    wz[[0, -1]] = 0.5

# Point-data arrays
rho = image.PointData["Density [gr/cm^3]"]
v = image.PointData["Velocity [cm/s]"]

# Reshape to image dimensions.
# VTK point ordering corresponds to z varying fastest.
rho = rho.reshape((nz, ny, nx))
v = v.reshape((nz, ny, nx, 3))

# SimpleAccretor at Position (4.857511e+13, 5.090694e+13, 5.000000e+13), Radius = 2.982200e-05
x0 = np.asarray([4.857511e+13, 5.090694e+13, 5.000000e+13])

# Coordinate differences, using broadcasting.
rx = x[None, None, :] - x0[0]
ry = y[None, :, None] - x0[1]
rz = z[:, None, None] - x0[2]

# Momentum density
px = rho * v[..., 0]
py = rho * v[..., 1]
pz = rho * v[..., 2]

Mdens = make_vector(px.ravel(), py.ravel(), pz.ravel())
output.PointData.append(Mdens, "Momentum density")

# Angular momentum density
Lx = ry * pz - rz * py
Ly = rz * px - rx * pz
Lz = rx * py - ry * px

Amd = make_vector(Lx.ravel(), Ly.ravel(), Lz.ravel())

output.PointData.append(Amd, "Angular momentum density")

# Tensor-product integration weights
W = (
    wz[:, None, None] *
    wy[None, :, None] *
    wx[None, None, :]
)

factor = dx * dy * dz

L = np.array([
    np.sum(Lx * W) * factor,
    np.sum(Ly * W) * factor,
    np.sum(Lz * W) * factor
])

print("Angular momentum:", L)""",
)

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=programmableFilter1)
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Set(
    Origin=[48575107574462.89, 50906936645507.81, 50000000000000.0],
    Normal=[0.0, 0.0, 1.0],
)

# init the 'Plane' selected for 'HyperTreeGridSlicer'
slice1.HyperTreeGridSlicer.Origin = [48575107574462.89, 50906936645507.81, 50000000000000.0]

# create a new 'Extract Block'
accretor = ExtractBlock(registrationName='Accretor', Input=OutputPort(reader,1))
accretor.Set(
    Assembly='Hierarchy',
    Selectors=['/Root/SimpleAccretor'],
)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from accretor
accretorDisplay = Show(accretor, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
accretorDisplay.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
accretorDisplay.ScaleTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
accretorDisplay.OpacityTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'Densitygrcm3'
densitygrcm3LUT = GetColorTransferFunction('Densitygrcm3')
densitygrcm3LUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=[
        # scalar, red, green, blue
        9.999996163704522e-20, 0.05639999999999999, 0.05639999999999999, 0.47,
        1.668618794746745e-17, 0.24300000000000013, 0.4603500000000004, 0.81,
        7.343078662909453e-16, 0.3568143826543521, 0.7450246485363142, 0.954367702893722,
        3.9506252452387274e-14, 0.6882, 0.93, 0.9179099999999999,
        2.9901484269261274e-13, 0.8994959551205902, 0.944646394975174, 0.7686567142818399,
        4.152811035857781e-12, 0.957107977357604, 0.8338185108985666, 0.5089156299842102,
        1.3980018529988463e-10, 0.9275207599610714, 0.6214389091739178, 0.31535705838676426,
        9.508378022069595e-09, 0.8, 0.3520000000000001, 0.15999999999999998,
        8.940991045077145e-07, 0.59, 0.07670000000000013, 0.11947499999999994,
    ],
    UseLogScale=1,
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
slice1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Density [gr/cm^3]'],
    LookupTable=densitygrcm3LUT,
    Assembly='Assembly',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [9.999996163704499e-20, 0.0, 0.5, 0.0, 8.940991045077148e-07, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [9.999996163704499e-20, 0.0, 0.5, 0.0, 8.940991045077148e-07, 1.0, 0.5, 0.0]

# show data from programmableFilter1
programmableFilter1Display = Show(programmableFilter1, renderView1, 'UniformGridRepresentation')

# trace defaults for the display properties.
programmableFilter1Display.Set(
    Representation='Outline',
    ColorArrayName=[None, ''],
)

# init the 'Plane' selected for 'SliceFunction'
programmableFilter1Display.SliceFunction.Origin = [48575111389160.16, 50906936645507.81, 50000000000000.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for densitygrcm3LUT in view renderView1
densitygrcm3LUTColorBar = GetScalarBar(densitygrcm3LUT, renderView1)
densitygrcm3LUTColorBar.Set(
    WindowLocation='Any Location',
    Orientation = 'Horizontal',
    Position=[0.5593700787401575, 0.03936395759717308],
    Title='Density [gr/cm^3]',
    ComponentTitle='',
    ScalarBarLength=0.36464566929133857,
    AllowOverlappingLabels=1,
)

# set color bar visibility
densitygrcm3LUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from slice1
slice1Display_1 = Show(slice1, renderView2, 'GeometryRepresentation')

# get color transfer function/color map for 'Velocitycms'
velocitycmsLUT = GetColorTransferFunction('Velocitycms')
velocitycmsLUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=[
        # scalar, red, green, blue
        0.0, 1.11641e-07, 0.0, 1.62551e-06,
        515819775.56828755, 0.0413146, 0.0619808, 0.209857,
        1031639551.1365751, 0.0185557, 0.101341, 0.350684,
        1547463437.150122, 0.00486405, 0.149847, 0.461054,
        2063283212.7184093, 0.0836345, 0.210845, 0.517906,
        2579102988.286697, 0.173222, 0.276134, 0.541793,
        3094922763.8549843, 0.259857, 0.343877, 0.535869,
        3610742539.423272, 0.362299, 0.408124, 0.504293,
        4126564657.945357, 0.468266, 0.468276, 0.468257,
        4642386201.005107, 0.582781, 0.527545, 0.374914,
        5158205976.573394, 0.691591, 0.585251, 0.274266,
        5674025752.141682, 0.784454, 0.645091, 0.247332,
        6189845527.70997, 0.862299, 0.710383, 0.27518,
        6705669413.7235155, 0.920863, 0.782923, 0.351563,
        7221489189.291802, 0.955792, 0.859699, 0.533541,
        7737308964.860091, 0.976162, 0.93433, 0.780671,
        8220890518.261018, 1.0, 1.0, 0.999983,
    ],
    NanColor=[1.0, 0.0, 0.0],
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
slice1Display_1.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity [cm/s]'],
    LookupTable=velocitycmsLUT,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display_1.ScaleTransferFunction.Points = [9.999996473106875e-20, 0.0, 0.5, 0.0, 8.551440235296542e-07, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display_1.OpacityTransferFunction.Points = [9.999996473106875e-20, 0.0, 0.5, 0.0, 8.551440235296542e-07, 1.0, 0.5, 0.0]

# show data from accretor
accretorDisplay_1 = Show(accretor, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
accretorDisplay_1.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
accretorDisplay_1.ScaleTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
accretorDisplay_1.OpacityTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# show data from programmableFilter1
programmableFilter1Display_1 = Show(programmableFilter1, renderView2, 'UniformGridRepresentation')

# trace defaults for the display properties.
programmableFilter1Display_1.Set(
    Representation='Outline',
    ColorArrayName=['POINTS', ''],
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
programmableFilter1Display_1.ScaleTransferFunction.Points = [9.999996473104078e-20, 0.0, 0.5, 0.0, 2.961973623283542e-06, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
programmableFilter1Display_1.OpacityTransferFunction.Points = [9.999996473104078e-20, 0.0, 0.5, 0.0, 2.961973623283542e-06, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
programmableFilter1Display_1.SliceFunction.Origin = [48575111389160.16, 50906936645507.81, 50000000000000.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for velocitycmsLUT in view renderView2
velocitycmsLUTColorBar = GetScalarBar(velocitycmsLUT, renderView2)
velocitycmsLUTColorBar.Set(
    WindowLocation='Any Location',
    Orientation = 'Horizontal',
    Position=[0.6223622047244095, 0.04240282685512368],
    Title='Velocity [cm/s]',
    ComponentTitle='Magnitude',
    ScalarBarLength=0.32999999999999996,
)

# set color bar visibility
velocitycmsLUTColorBar.Visibility = 1

# show color legend
slice1Display_1.SetScalarBarVisibility(renderView2, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView3'
# ----------------------------------------------------------------

# show data from programmableFilter1
programmableFilter1Display_2 = Show(programmableFilter1, renderView3, 'UniformGridRepresentation')

# trace defaults for the display properties.
programmableFilter1Display_2.Set(
    Representation='Outline',
    ColorArrayName=['POINTS', ''],
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
programmableFilter1Display_2.ScaleTransferFunction.Points = [9.999996473104078e-20, 0.0, 0.5, 0.0, 2.961973623283542e-06, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
programmableFilter1Display_2.OpacityTransferFunction.Points = [9.999996473104078e-20, 0.0, 0.5, 0.0, 2.961973623283542e-06, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
programmableFilter1Display_2.SliceFunction.Origin = [48575111389160.16, 50906936645507.81, 50000000000000.0]

# show data from slice1
slice1Display_2 = Show(slice1, renderView3, 'GeometryRepresentation')

# get color transfer function/color map for 'Momentumdensity'
momentumdensityLUT = GetColorTransferFunction('Momentumdensity')
momentumdensityLUT.Set(
    RGBPoints=GenerateRGBPoints(
        range_min=0.0,
        range_max=4908.336132557782,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
slice1Display_2.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Momentum density'],
    LookupTable=momentumdensityLUT,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display_2.ScaleTransferFunction.Points = [9.999996473106875e-20, 0.0, 0.5, 0.0, 8.551440235296542e-07, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display_2.OpacityTransferFunction.Points = [9.999996473106875e-20, 0.0, 0.5, 0.0, 8.551440235296542e-07, 1.0, 0.5, 0.0]

# show data from accretor
accretorDisplay_2 = Show(accretor, renderView3, 'GeometryRepresentation')

# trace defaults for the display properties.
accretorDisplay_2.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
accretorDisplay_2.ScaleTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
accretorDisplay_2.OpacityTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for momentumdensityLUT in view renderView3
momentumdensityLUTColorBar = GetScalarBar(momentumdensityLUT, renderView3)
momentumdensityLUTColorBar.Set(
    WindowLocation='Any Location',
    Orientation = 'Horizontal',
    Position=[0.594565560821485, 0.08000000000000007],
    Title='Momentum density',
    ComponentTitle='Magnitude',
    ScalarBarLength=0.32999999999999985,
)

# set color bar visibility
momentumdensityLUTColorBar.Visibility = 1

# show color legend
slice1Display_2.SetScalarBarVisibility(renderView3, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView4'
# ----------------------------------------------------------------

# show data from programmableFilter1
programmableFilter1Display_3 = Show(programmableFilter1, renderView4, 'UniformGridRepresentation')

# trace defaults for the display properties.
programmableFilter1Display_3.Set(
    Representation='Outline',
    ColorArrayName=['POINTS', ''],
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
programmableFilter1Display_3.ScaleTransferFunction.Points = [9.999996473104078e-20, 0.0, 0.5, 0.0, 2.961973623283542e-06, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
programmableFilter1Display_3.OpacityTransferFunction.Points = [9.999996473104078e-20, 0.0, 0.5, 0.0, 2.961973623283542e-06, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
programmableFilter1Display_3.SliceFunction.Origin = [48575111389160.16, 50906936645507.81, 50000000000000.0]

# show data from slice1
slice1Display_3 = Show(slice1, renderView4, 'GeometryRepresentation')

# get color transfer function/color map for 'EnergyDensityergcm3'
energyDensityergcm3LUT = GetColorTransferFunction('EnergyDensityergcm3')
energyDensityergcm3LUT.Set(
    RGBPoints=GenerateRGBPoints(
        range_min=1.3423333959252284e-05,
        range_max=20881997766670.5,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
slice1Display_3.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Energy Density [erg/cm^3]'],
    LookupTable=energyDensityergcm3LUT,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display_3.ScaleTransferFunction.Points = [9.999996473106875e-20, 0.0, 0.5, 0.0, 8.551440235296542e-07, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display_3.OpacityTransferFunction.Points = [9.999996473106875e-20, 0.0, 0.5, 0.0, 8.551440235296542e-07, 1.0, 0.5, 0.0]

# show data from accretor
accretorDisplay_3 = Show(accretor, renderView4, 'GeometryRepresentation')

# trace defaults for the display properties.
accretorDisplay_3.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
accretorDisplay_3.ScaleTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
accretorDisplay_3.OpacityTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for energyDensityergcm3LUT in view renderView4
energyDensityergcm3LUTColorBar = GetScalarBar(energyDensityergcm3LUT, renderView4)
energyDensityergcm3LUTColorBar.Set(
    WindowLocation='Any Location',
    Orientation = 'Horizontal',
    Position=[0.5961453396524486, 0.046431095406360395],
    Title='Energy Density [erg/cm^3]',
    ComponentTitle='',
    ScalarBarLength=0.3300000000000002,
)

# set color bar visibility
energyDensityergcm3LUTColorBar.Visibility = 1

# show color legend
slice1Display_3.SetScalarBarVisibility(renderView4, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView5'
# ----------------------------------------------------------------

# show data from programmableFilter1
programmableFilter1Display_4 = Show(programmableFilter1, renderView5, 'UniformGridRepresentation')

# trace defaults for the display properties.
programmableFilter1Display_4.Set(
    Representation='Outline',
    ColorArrayName=['POINTS', ''],
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
programmableFilter1Display_4.ScaleTransferFunction.Points = [9.999996473104078e-20, 0.0, 0.5, 0.0, 2.961973623283542e-06, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
programmableFilter1Display_4.OpacityTransferFunction.Points = [9.999996473104078e-20, 0.0, 0.5, 0.0, 2.961973623283542e-06, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
programmableFilter1Display_4.SliceFunction.Origin = [48575111389160.16, 50906936645507.81, 50000000000000.0]

# show data from accretor
accretorDisplay_4 = Show(accretor, renderView5, 'GeometryRepresentation')

# trace defaults for the display properties.
accretorDisplay_4.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
accretorDisplay_4.ScaleTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
accretorDisplay_4.OpacityTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# show data from slice1
slice1Display_4 = Show(slice1, renderView5, 'GeometryRepresentation')

# get color transfer function/color map for 'ThermalEnergyDensityergcm3'
thermalEnergyDensityergcm3LUT = GetColorTransferFunction('ThermalEnergyDensityergcm3')
thermalEnergyDensityergcm3LUT.Set(
    RGBPoints=GenerateRGBPoints(
        range_min=1.342333395925226e-05,
        range_max=5735108516751.374,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
slice1Display_4.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Thermal Energy Density [erg/cm^3]'],
    LookupTable=thermalEnergyDensityergcm3LUT,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display_4.ScaleTransferFunction.Points = [9.999996473106875e-20, 0.0, 0.5, 0.0, 8.551440235296542e-07, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display_4.OpacityTransferFunction.Points = [9.999996473106875e-20, 0.0, 0.5, 0.0, 8.551440235296542e-07, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for thermalEnergyDensityergcm3LUT in view renderView5
thermalEnergyDensityergcm3LUTColorBar = GetScalarBar(thermalEnergyDensityergcm3LUT, renderView5)
thermalEnergyDensityergcm3LUTColorBar.Set(
    WindowLocation='Any Location',
    Orientation = 'Horizontal',
    Position=[0.5977251184834125, 0.0340636042402826],
    Title='Thermal Energy Density [erg/cm^3]',
    ComponentTitle='',
    ScalarBarLength=0.32999999999999996,
)

# set color bar visibility
thermalEnergyDensityergcm3LUTColorBar.Visibility = 1

# show color legend
slice1Display_4.SetScalarBarVisibility(renderView5, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView6'
# ----------------------------------------------------------------

# show data from accretor
accretorDisplay_5 = Show(accretor, renderView6, 'GeometryRepresentation')

# trace defaults for the display properties.
accretorDisplay_5.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
accretorDisplay_5.ScaleTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
accretorDisplay_5.OpacityTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 14.801953315734863, 1.0, 0.5, 0.0]

# show data from slice1
slice1Display_5 = Show(slice1, renderView6, 'GeometryRepresentation')

# get color transfer function/color map for 'Angularmomentumdensity'
angularmomentumdensityLUT = GetColorTransferFunction('Angularmomentumdensity')
angularmomentumdensityLUT.Set(
    RGBPoints=GenerateRGBPoints(
        range_min=0.0,
        range_max=520158334040.5685,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
slice1Display_5.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Angular momentum density'],
    LookupTable=angularmomentumdensityLUT,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display_5.ScaleTransferFunction.Points = [9.999996473106875e-20, 0.0, 0.5, 0.0, 8.551440235296542e-07, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display_5.OpacityTransferFunction.Points = [9.999996473106875e-20, 0.0, 0.5, 0.0, 8.551440235296542e-07, 1.0, 0.5, 0.0]

# show data from programmableFilter1
programmableFilter1Display_5 = Show(programmableFilter1, renderView6, 'UniformGridRepresentation')

# trace defaults for the display properties.
programmableFilter1Display_5.Set(
    Representation='Outline',
    ColorArrayName=['POINTS', ''],
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
programmableFilter1Display_5.ScaleTransferFunction.Points = [9.999996473104078e-20, 0.0, 0.5, 0.0, 2.961973623283542e-06, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
programmableFilter1Display_5.OpacityTransferFunction.Points = [9.999996473104078e-20, 0.0, 0.5, 0.0, 2.961973623283542e-06, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
programmableFilter1Display_5.SliceFunction.Origin = [48575111389160.16, 50906936645507.81, 50000000000000.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for angularmomentumdensityLUT in view renderView6
angularmomentumdensityLUTColorBar = GetScalarBar(angularmomentumdensityLUT, renderView6)
angularmomentumdensityLUTColorBar.Set(
    WindowLocation='Any Location',
    Orientation = 'Horizontal',
    Position=[0.580347551342812, 0.04289752650176673],
    Title='Angular momentum density',
    ComponentTitle='Magnitude',
    ScalarBarLength=0.33000000000000007,
)

# set color bar visibility
angularmomentumdensityLUTColorBar.Visibility = 1

# show color legend
slice1Display_5.SetScalarBarVisibility(renderView6, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'ThermalEnergyDensityergcm3'
thermalEnergyDensityergcm3PWF = GetOpacityTransferFunction('ThermalEnergyDensityergcm3')
thermalEnergyDensityergcm3PWF.Set(
    Points=[1.342333395925226e-05, 0.0, 0.5, 0.0, 5735108516751.374, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'Velocitycms'
velocitycmsPWF = GetOpacityTransferFunction('Velocitycms')
velocitycmsPWF.Set(
    Points=[0.0, 0.0, 0.5, 0.0, 8220890518.261018, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'Densitygrcm3'
densitygrcm3PWF = GetOpacityTransferFunction('Densitygrcm3')
densitygrcm3PWF.Set(
    Points=[9.999996163704499e-20, 0.0, 0.5, 0.0, 7.383330380613588e-18, 0.0, 0.5, 0.0, 8.192588101331325e-18, 1.0, 0.5, 0.0, 8.940991045077148e-07, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'EnergyDensityergcm3'
energyDensityergcm3PWF = GetOpacityTransferFunction('EnergyDensityergcm3')
energyDensityergcm3PWF.Set(
    Points=[1.3423333959252284e-05, 0.0, 0.5, 0.0, 20881997766670.5, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'Momentumdensity'
momentumdensityPWF = GetOpacityTransferFunction('Momentumdensity')
momentumdensityPWF.Set(
    Points=[0.0, 0.0, 0.5, 0.0, 4908.336132557782, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'Angularmomentumdensity'
angularmomentumdensityPWF = GetOpacityTransferFunction('Angularmomentumdensity')
angularmomentumdensityPWF.Set(
    Points=[0.0, 0.0, 0.5, 0.0, 520158334040.5685, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper
timeKeeper1.SuppressedTimeSources = [reader, resampleToImage1]

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.Set(
    ViewModules=[renderView1, renderView2, renderView3, renderView4, renderView5, renderView6],
    Cues=timeAnimationCue1,
    AnimationTime=682111.484339536,
    StartTime=682111.484339536,
    EndTime=682112.484339536,
    PlayMode='Snap To TimeSteps',
)

# ----------------------------------------------------------------
# restore active source
SetActiveSource(slice1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup camera links
AddCameraLink(renderView1, renderView4, 'CameraLink0')
AddCameraLink(renderView1, renderView5, 'CameraLink1')
AddCameraLink(renderView1, renderView2, 'CameraLink2')
AddCameraLink(renderView1, renderView3, 'CameraLink3')
AddCameraLink(renderView1, renderView6, 'CameraLink4')


