# state file generated using paraview version 6.1.1
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
materialLibrary1 = GetMaterialLibrary()

renderView1 = GetRenderView()
renderView1.Set(
    CenterOfRotation=[48596391034880.0, 50928556179456.0, 50002793594880.0],
    CameraPosition=[46728175275234.4, 47312228148177.38, 48611337824725.445],
    CameraFocalPoint=[46738475360820.055, 47331829696592.125, 48618913426822.555],
    CameraViewUp=[0.1826575486013235, 0.2646634591923189, -0.9468840865213182],
    CameraViewAngle=4.6618150684931505,
    EnableRayTracing=0,
    BackEnd='OSPRay raycaster',
    Shadows=1,
    AmbientSamples=1,
    SamplesPerPixel = 50,
    OSPRayMaterialLibrary=materialLibrary1,
    OrientationAxesVisibility = 0
)

directory_JeanLaptop = '/local/data/Walder/CygX-1_M1=14.8_M2=19.2_V=750_ML=1-6_Gamma=1.01_RBH=20_HR/DATA/'
directory_r740gpu02 = '/Xnfs/highenergy/rwalder/Simulations/CygX-1_M1=14.8_M2=19.2_V=750_ML=1-6_Gamma=1.01_RBH=20_HR/'

fname = 'DATA.CX-1_M1=14.8_M2=19.2_V=750_ML=1-6_G=1.01_RBH=20RG_.NT=0000008600.Time=189.4754123165H.amr5'

# create a new directory name for Rolf and Doris's laptops, or use Jean's laptop, or use r740gpu02
# watch for the last "/" in the directory name

#filename = directory_r740gpu02 + fname
filename = directory_JeanLaptop + fname

reader = AMAZEReader(registrationName='reader', FileNames=filename)
reader.Set(
    PointArrayStatus=['Density', 'Velocity'],
    Level=100,
)

reader.UpdatePipeline()
Global_bounds = reader.GetDataInformation().GetBounds()
mid_plane = [(Global_bounds[1] - Global_bounds[0]) * 0.5,
             (Global_bounds[1] - Global_bounds[0]) * 0.5,
             (Global_bounds[1] - Global_bounds[0]) * 0.5
             ]
             

# create a new 'Contour'
contour1 = Contour(registrationName='Contour1', Input=reader)
contour1.Set(
    ContourBy=['POINTS', 'Density [gr/cm^3]'],
    Isosurfaces=[1e-10],
)

# create a new 'Mask Points'
maskPoints1 = MaskPoints(registrationName='MaskPoints1', Input=contour1)
maskPoints1.Set(
    MaximumNumberofPoints=100,
    ProportionallyDistributeMaximumNumberOfPoints=1,
    RandomSampling=1,
    GenerateVertices=1,
)

# create a new 'Point Source'
pointSource1 = PointSource(registrationName='PointSource1')
pointSource1.Set(
    Center=[48583300000000.0, 50919680000000.0, 50000000000000.0],
    NumberOfPoints=100,
    Radius=100000000000.0,
)

# create a new 'Stream Tracer With Custom Source'
streamTracerWithCustomSource1 = StreamTracerWithCustomSource(registrationName='StreamTracerWithCustomSource1', Input=reader,
    SeedSource=pointSource1)
streamTracerWithCustomSource1.Set(
    Vectors=['POINTS', 'Velocity [cm/s]'],
    MaximumSteps=5000,
    MaximumStreamlineLength=1e+17,
)

# create a new 'Contour'
contour2 = Contour(registrationName='Contour2', Input=streamTracerWithCustomSource1)
contour2.Set(
    ContourBy=['POINTS', 'IntegrationTime'],
    Isosurfaces=[-3000.0, 3000.0, -2000.0, 2000.0, -1000.0, 1000.0],
)

# create a new 'Glyph'
glyph2 = Glyph(registrationName='Glyph2', Input=contour2,
    GlyphType='Arrow')
glyph2.Set(
    OrientationArray=['POINTS', 'Velocity [cm/s]'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=37500000000.0,
    GlyphMode='All Points',
)

# init the 'Arrow' selected for 'GlyphType'
glyph2.GlyphType.Set(
    TipResolution=12,
    ShaftResolution=12,
)

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=maskPoints1,
    GlyphType='Arrow')
glyph1.Set(
    OrientationArray=['POINTS', 'Velocity [cm/s]'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=21552431104.0,
    GlyphMode='All Points',
)

# create a new 'Clip'
halfdomain = Clip(registrationName='half-domain', Input=reader)
halfdomain.Invert = 0

# init the 'Plane' selected for 'ClipType'
halfdomain.ClipType.Set(
    Origin=mid_plane,
    Normal=[0.0, 0.0, 1.0],
)

# init the 'Plane' selected for 'HyperTreeGridClipper'
halfdomain.HyperTreeGridClipper.Origin = [50000000000000.0, 50000000000000.0, 50000000000000.0]

# create a new 'Clip'
clip3 = Clip(registrationName='Clip3', Input=halfdomain)
clip3.Set(
    ClipType='Scalar',
    Scalars=['POINTS', 'Density [gr/cm^3]'],
    Value=4e-14,
    Invert=0,
)

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=halfdomain)
clip1.Set(
    ClipType='Scalar',
    Scalars=['POINTS', 'Density [gr/cm^3]'],
    Value=5e-15,
)

# create a new 'Clip'
clip2 = Clip(registrationName='Clip2', Input=clip1)
clip2.Set(
    ClipType='Scalar',
    Scalars=['POINTS', 'Density [gr/cm^3]'],
    Value=1e-15,
    Invert=0,
)

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=reader)
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Set(
    Origin=[50000000000000.0, 50000000000000.0, 50000000000000.0],
    Normal=[0.0, 0.0, 1.0],
)

# show data from clip3
clip3Display = Show(clip3, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'Densitygrcm3'
densitygrcm3LUT = GetColorTransferFunction('Densitygrcm3')
densitygrcm3LUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=[
        # scalar, red, green, blue
        2.0000000000000002e-16, 0.05639999999999999, 0.05639999999999999, 0.47,
        1.1060900170346836e-14, 0.24300000000000013, 0.4603500000000004, 0.81,
        2.1509520034506664e-13, 0.3568143826543521, 0.7450246485363142, 0.954367702893722,
        4.89671651274604e-12, 0.6882, 0.93, 0.9179099999999999,
        2.3946005290892195e-11, 0.8994959551205902, 0.944646394975174, 0.7686567142818399,
        1.8849182633410137e-10, 0.957107977357604, 0.8338185108985666, 0.5089156299842102,
        2.9708959530143666e-09, 0.9275207599610714, 0.6214389091739178, 0.31535705838676426,
        8.128330351117764e-08, 0.8, 0.3520000000000001, 0.15999999999999998,
        2.8670558469571843e-06, 0.59, 0.07670000000000013, 0.11947499999999994,
    ],
    UseLogScale=1,
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
clip3Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Density [gr/cm^3]'],
    LookupTable=densitygrcm3LUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
clip3Display.ScaleTransferFunction.Points = [9.99999999981972e-16, 0.0, 0.5, 0.0, 3.750342361890499e-07, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
clip3Display.OpacityTransferFunction.Points = [9.99999999981972e-16, 0.0, 0.5, 0.0, 3.750342361890499e-07, 1.0, 0.5, 0.0]

# show data from streamTracerWithCustomSource1
streamTracerWithCustomSource1Display = Show(streamTracerWithCustomSource1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'Velocitycms'
velocitycmsLUT = GetColorTransferFunction('Velocitycms')
velocitycmsLUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=[
        # scalar, red, green, blue
        1000.0, 1.11641e-07, 0.0, 1.62551e-06,
        6275437.255, 0.0413146, 0.0619808, 0.209857,
        12549874.51, 0.0185557, 0.101341, 0.350684,
        18824361.7645, 0.00486405, 0.149847, 0.461054,
        25098799.0195, 0.0836345, 0.210845, 0.517906,
        31373236.274499997, 0.173222, 0.276134, 0.541793,
        37647673.5295, 0.259857, 0.343877, 0.535869,
        43922110.784499995, 0.362299, 0.408124, 0.504293,
        50196576.539215, 0.468266, 0.468276, 0.468257,
        56471035.29400001, 0.582781, 0.527545, 0.374914,
        62745472.548999995, 0.691591, 0.585251, 0.274266,
        69019909.804, 0.784454, 0.645091, 0.247332,
        75294347.05900002, 0.862299, 0.710383, 0.27518,
        81568834.3135, 0.920863, 0.782923, 0.351563,
        87843271.5685, 0.955792, 0.859699, 0.533541,
        94117708.8235, 0.976162, 0.93433, 0.780671,
        100000000.0, 1.0, 1.0, 0.999983,
    ],
    NanColor=[1.0, 0.0, 0.0],
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
streamTracerWithCustomSource1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity [cm/s]'],
    LookupTable=velocitycmsLUT,
    LineWidth=.3,
    RenderLinesAsTubes=0,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
streamTracerWithCustomSource1Display.ScaleTransferFunction.Points = [7.959786682934369e-18, 0.0, 0.5, 0.0, 1.2647800773757307e-06, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
streamTracerWithCustomSource1Display.OpacityTransferFunction.Points = [7.959786682934369e-18, 0.0, 0.5, 0.0, 1.2647800773757307e-06, 1.0, 0.5, 0.0]

# show data from glyph2
glyph2Display = Show(glyph2, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph2Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity [cm/s]'],
    LookupTable=velocitycmsLUT,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
glyph2Display.ScaleTransferFunction.Points = [-5555555.555555556, 0.0, 0.5, 0.0, -1111111.111111112, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
glyph2Display.OpacityTransferFunction.Points = [-5555555.555555556, 0.0, 0.5, 0.0, -1111111.111111112, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for densitygrcm3LUT in view renderView1
densitygrcm3LUTColorBar = GetScalarBar(densitygrcm3LUT, renderView1)
densitygrcm3LUTColorBar.Set(
    WindowLocation='Upper Right Corner',
    Title='Density [gr/cm^3]',
    ComponentTitle='',
)

# set color bar visibility
densitygrcm3LUTColorBar.Visibility = 1

# get color legend/bar for velocitycmsLUT in view renderView1
velocitycmsLUTColorBar = GetScalarBar(velocitycmsLUT, renderView1)
velocitycmsLUTColorBar.Set(
    Title='Velocity [cm/s]',
    ComponentTitle='Magnitude',
)

# set color bar visibility
velocitycmsLUTColorBar.Visibility = 1

# show color legend
clip3Display.SetScalarBarVisibility(renderView1, True)

# show color legend
streamTracerWithCustomSource1Display.SetScalarBarVisibility(renderView1, True)

# show color legend
glyph2Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'Velocitycms'
velocitycmsPWF = GetOpacityTransferFunction('Velocitycms')
velocitycmsPWF.Set(
    Points=[1000.0, 0.0, 0.5, 0.0, 100000000.0, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'Densitygrcm3'
densitygrcm3PWF = GetOpacityTransferFunction('Densitygrcm3')
densitygrcm3PWF.Set(
    Points=[2e-16, 0.0, 0.5, 0.0, 2.233550340855008e-16, 0.0, 0.5, 0.0, 2.259500339775932e-16, 1.0, 0.5, 0.0, 2.8670558469571783e-06, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation scene

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.Set(
    ViewModules=renderView1,
    Cues=timeAnimationCue1,
    AnimationTime=681420.9843395391,
    StartTime=3954516.40220226,
    EndTime=3954517.40220226,
    PlayMode='Snap To TimeSteps',
)

# ----------------------------------------------------------------
# restore active source
SetActiveSource(glyph2)
if __name__ == '__main__':
  print("rendering image")
  SaveScreenshot("StreamlinesAroundDisk.png")
