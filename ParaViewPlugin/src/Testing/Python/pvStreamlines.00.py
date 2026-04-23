# state file generated using paraview version 6.1.0
import paraview
from math import sqrt

paraview.compatibility.major = 6
paraview.compatibility.minor = 1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


renderView1 = GetRenderView()
renderView1.Set(
    ViewSize=[1024,1024],
    CenterOfRotation=[100540323201024.0, 100501953708032.0, 100000000376832.0],
    CameraPosition=[103635548547907.8, 99500414158541.61, 102007840290891.8],
    CameraFocalPoint=[94084105875000.77, 106104842719205.73, 93073213528864.52],
    CameraViewUp=[-0.4714061580718542, 0.3889935826356, 0.7914924047647297],
    OrientationAxesVisibility = 0
)

# select your own filename, local, or remote (if remote, create the connection to server before loading this script)

fname = "/local/data/Walder/DATA.CygX-1_M1=14.8_M2=19.2_D=2e14_V=850_ML=1-6_Gamma=1.6666_R=256.COOL_HR.NT=0000021000.Time=3921312.6324999821S.amr5"

reader = AMAZEReader(registrationName='DATA', FileNames=fname)
reader.Set(
    PointArrayStatus=['Density', 'Velocity'],
    LevelSet=[0, 14],
)
reader.UpdatePipeline()

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=reader)
slice1.SliceType.Set(
    Origin=[1.e14, 1.e14, 1.e14],
    Normal=[0.0, 0.0, 1.0],  # a Z slice
)

# create a new 'Extract Block' to place a disk inside it
# we measure its radius and center. The Selectors get a particular refinement level
extractBlock1 = ExtractBlock(registrationName='ExtractBlock1', Input=slice1)
extractBlock1.Set(
    Assembly='Hierarchy',
    Selectors=['/Root/Block8'],
)
extractBlock1.UpdatePipeline()

eb = FindSource("ExtractBlock1")
bbb = eb.GetDataInformation().GetBounds()
diskCenter = [bbb[0]+0.5*(bbb[1] - bbb[0]),
              bbb[2]+0.5*(bbb[3] - bbb[2]),
              bbb[4]+0.5*(bbb[5] - bbb[4])]
diskRadius = 0.5 * sqrt((bbb[1] - bbb[0])*(bbb[1] - bbb[0]) +
                        (bbb[3] - bbb[2])*(bbb[3] - bbb[2]))
# create a new 'Disk'
disk1 = Disk(registrationName='Disk1')
disk1.Set(
    OuterRadius=diskRadius,
    RadialResolution=128,
    CircumferentialResolution=128,
)
axisAlignedTransform1 = AxisAlignedTransform(registrationName='AxisAlignedTransform1', Input=disk1)
axisAlignedTransform1.Translation=diskCenter

# create a new disk-shaped cut through the slice to use as input to the masking filter below
resampleWithDataset1 = ResampleWithDataset(registrationName='ResampleWithDataset1', SourceDataArrays=slice1,
    DestinationMesh=axisAlignedTransform1)

# create a new 'Mask Points'
maskPoints1 = MaskPoints(registrationName='MaskPoints1', Input=resampleWithDataset1)
maskPoints1.Set(
    MaximumNumberofPoints=500,
    RandomSampling=1,
    RandomSamplingMode='Random Sampling',
    GenerateVertices=1,
)

# display small velocity arrows on the disk, above the slice
glyph1 = Glyph(registrationName='Glyph1', Input=maskPoints1,
    GlyphType='Arrow')
glyph1.Set(
    OrientationArray=['POINTS', 'Velocity [cm/s]'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=60000000000.0,
)

reader_1Display = Show(OutputPort(reader, 1), renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'vtkBlockColors'
vtkBlockColorsLUT = GetColorTransferFunction('vtkBlockColors')
vtkBlockColorsLUT.Set(
    InterpretValuesAsCategories=1,
    Annotations=['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11'],
    ActiveAnnotatedValues=['0', '1'],
    IndexedColors=[1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.63, 0.63, 1.0, 0.67, 0.5, 0.33, 1.0, 0.5, 0.75, 0.53, 0.35, 0.7, 1.0, 0.75, 0.5],
)

# trace defaults for the display properties.
reader_1Display.Set(
    Representation='Surface',
    ColorArrayName=['FIELD', 'vtkBlockColors'],
    LookupTable=vtkBlockColorsLUT,
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)


# show data from slice1
slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'Densitygrcm3'
densitygrcm3LUT = GetColorTransferFunction('Densitygrcm3')
densitygrcm3LUT.Set(
    RGBPoints=[
        # scalar, red, green, blue
        9.999827218728697e-20, 0.0564, 0.0564, 0.47,
        1.0557671474038171e-17, 0.243, 0.46035, 0.81,
        3.311932956042457e-16, 0.356814, 0.745025, 0.954368,
        1.247577721203983e-14, 0.6882, 0.93, 0.91791,
        7.878922518882417e-14, 0.899496, 0.944646, 0.768657,
        8.647917119777687e-13, 0.957108, 0.833819, 0.508916,
        2.1255780705225206e-11, 0.927521, 0.621439, 0.315357,
        9.912141753803464e-10, 0.8, 0.352, 0.16,
        6.207849265864068e-08, 0.59, 0.0767, 0.119475,
    ],
    UseLogScale=1,
    ScalarRangeInitialized=1.0,
)

# show slice with or without LIC
slice1Display.Set(
    Representation='Surface LIC',
    #Representation='Surface',
    ColorArrayName=['POINTS', 'Density [gr/cm^3]'],
    LookupTable=densitygrcm3LUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
slice1Display.ScaleTransferFunction.Points = [9.999827218728682e-20, 0.0, 0.5, 0.0, 6.207849265864068e-08, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [9.999827218728682e-20, 0.0, 0.5, 0.0, 6.207849265864068e-08, 1.0, 0.5, 0.0]


# show data from glyph1
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', ''],
    DisableLighting = 0
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [8.783839289829877e-15, 0.0, 0.5, 0.0, 4.136151274181865e-09, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [8.783839289829877e-15, 0.0, 0.5, 0.0, 4.136151274181865e-09, 1.0, 0.5, 0.0]


# get color legend/bar for densitygrcm3LUT in view renderView1
densitygrcm3LUTColorBar = GetScalarBar(densitygrcm3LUT, renderView1)
densitygrcm3LUTColorBar.Set(
    WindowLocation='Upper Right Corner',
    Title='Density [gr/cm^3]',
    ComponentTitle='',
)

# set color bar visibility
densitygrcm3LUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'Densitygrcm3'
densitygrcm3PWF = GetOpacityTransferFunction('Densitygrcm3')
densitygrcm3PWF.Set(
    Points=[9.999827218728682e-20, 0.0, 0.5, 0.0, 6.207849265864068e-08, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

SetActiveSource(slice1)
HideInteractiveWidgets()

SaveScreenshot("SliceLICwithArrows.png")
