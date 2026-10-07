# state file generated using paraview version 6.2.0
import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

renderView1 = GetRenderView()

if __name__ == '__main__':
    renderView1.ViewSize=[1920,1080]
renderView1.Set(
    ViewSize=[1920, 1080],
    OrientationAxesVisibility=0,
    CenterOfRotation=[49942858301440.0, 50042690338816.0, 49960265449472.0],
    CameraPosition=[49942858301440.0, 50042690338816.0, 243438390883037.03],
    CameraFocalPoint=[49942858301440.0, 50042690338816.0, -43091813706817.64],
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

# create a new 'Axis-Aligned Slice'
axisAlignedSlice1 = AxisAlignedSlice(registrationName='AxisAlignedSlice1', Input=reader)
axisAlignedSlice1.Level = 18

# init the 'Axis Aligned Plane' selected for 'CutFunction'
axisAlignedSlice1.CutFunction.Set(
    Origin=[50000000000000.0, 50000000000000.0, 50000000000000.0],
    Normal=[0.0, 0.0, 1.0],
)

# create a new 'Programmable Filter'
programmableFilter2 = ProgrammableFilter(registrationName='Co-Rotating Velocity', Input=axisAlignedSlice1)
programmableFilter2.Set(
    CopyArrays=1,
    PythonPath='',
    Script="""import numpy as np
from vtkmodules.vtkCommonDataModel import vtkOverlappingAMR
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk

######################################################
# define my computation for a single cartesian grid
def compute_corotating_velocity(grid, output_grid):
        # --------------------------------------------------------------
        # Get velocity
        # --------------------------------------------------------------

        velocity_vtk = grid.GetPointData().GetArray("Velocity [cm/s]")

        if velocity_vtk is None:
            raise RuntimeError(
                f"No \'velocity\' array in AMR level {level}, block {block}"
            )

        velocity = vtk_to_numpy(velocity_vtk)

        # Check dimensions
        if velocity.ndim != 2 or velocity.shape[1] != 3:
            raise RuntimeError(
                f"Velocity must be N x 3, got {velocity.shape}"
            )

        # --------------------------------------------------------------
        # Coordinates
        # --------------------------------------------------------------

        npoints = grid.GetNumberOfPoints()

        coords = np.empty((npoints, 3), dtype=np.float64)

        for i in range(npoints):
            grid.GetPoint(i, coords[i])

        # --------------------------------------------------------------
        # Corotating velocity
        #
        # Vcorot = V - Omega x (Coords - Pos_CM)
        # --------------------------------------------------------------

        r = coords - Pos_CM

        omega_cross_r = np.cross(Omega, r)

        velocity_corot = velocity - omega_cross_r

        # --------------------------------------------------------------
        # Add result to output
        # --------------------------------------------------------------

        corot_vtk = numpy_to_vtk(
            velocity_corot,
            deep=True
        )

        corot_vtk.SetName("Co-Rotating Velocity")

        output_grid.GetPointData().AddArray(corot_vtk)
        
######################################################

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------

Omega = np.array([0.0, 0.0, 1.2983630952380953e-05])
Pos_CM = np.array([1.0e+14, 1.0e+14, 1.0e+14])

# ----------------------------------------------------------------------
# Input / output of the nested AMR structures
# ----------------------------------------------------------------------

input_amr = self.GetInputDataObject(0, 0)
output_amr = self.GetOutputDataObject(0)

if not isinstance(input_amr, vtkOverlappingAMR):
    raise RuntimeError(
        f"Expected vtkOverlappingAMR, got {input_amr.GetClassName()}"
    )

output_amr.CopyStructure(input_amr)

# ----------------------------------------------------------------------
# Process every AMR block
# ----------------------------------------------------------------------

for level in range(input_amr.GetNumberOfLevels()):

    for block in range(input_amr.GetNumberOfBlocks(level)):

        grid = input_amr.GetDataSetAsImageData(level, block)

        if grid is None:
            continue

        # Make a copy of the block
        output_grid = grid.NewInstance()
        output_grid.ShallowCopy(grid)

        compute_corotating_velocity(grid, output_grid)

        # Put block into AMR output
        output_amr.SetDataSet(level, block, output_grid)""",
)

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=programmableFilter2)

# create a new 'Mask Points'
maskPoints1 = MaskPoints(registrationName='MaskPoints1', Input=resampleToImage1)
maskPoints1.Set(
    OnRatio=4,
    MaximumNumberofPoints=600,
    RandomSampling=1,
    RandomSamplingMode='Uniform Spatial Distribution (Bounds Based)',
    GenerateVertices=1,
)

# create a new 'Glyph'
VelocityGlyphs = Glyph(registrationName='Velocity Glyphs', Input=maskPoints1,
    GlyphType='Arrow')
VelocityGlyphs.Set(
    OrientationArray=['POINTS', 'Velocity [cm/s]'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=2500000000000.0,
    GlyphMode='All Points',
)

# create a new 'Glyph'
coRotatingVelocityGlyphs = Glyph(registrationName='Co-Rotating Velocity Glyphs', Input=maskPoints1,
    GlyphType='Arrow')
coRotatingVelocityGlyphs.Set(
    OrientationArray=['POINTS', 'Co-Rotating Velocity'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=2500000000000.0,
    GlyphMode='All Points',
)


reader_1Display = Show(OutputPort(reader, 1), renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
reader_1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)


# show data from axisAlignedSlice1
axisAlignedSlice1Display = Show(axisAlignedSlice1, renderView1, 'AMRRepresentation')

# get color transfer function/color map for 'Densitygrcm3'
densitygrcm3LUT = GetColorTransferFunction('Densitygrcm3')
densitygrcm3LUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=[
        # scalar, red, green, blue
        3.0993905462995587e-18, 0.05639999999999999, 0.05639999999999999, 0.47,
        1.8420972073731352e-17, 0.24300000000000013, 0.4603500000000004, 0.81,
        6.88244792114201e-17, 0.3568143826543521, 0.7450246485363142, 0.954367702893722,
        2.757832323063568e-16, 0.6882, 0.93, 0.9179099999999999,
        5.581213762189007e-16, 0.8994959551205902, 0.944646394975174, 0.7686567142818399,
        1.3954294817763714e-15, 0.957107977357604, 0.8338185108985666, 0.5089156299842102,
        4.74911010586632e-15, 0.9275207599610714, 0.6214389091739178, 0.31535705838676426,
        2.0648925237405664e-14, 0.8, 0.3520000000000001, 0.15999999999999998,
        1.0050345896692112e-13, 0.59, 0.07670000000000013, 0.11947499999999994,
    ],
    UseLogScale=1,
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
axisAlignedSlice1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Density [gr/cm^3]'],
    LookupTable=densitygrcm3LUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
axisAlignedSlice1Display.ScaleTransferFunction.Points = [9.999996164020771e-20, 0.0, 0.5, 0.0, 8.940991045077148e-07, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
axisAlignedSlice1Display.OpacityTransferFunction.Points = [9.999996164020771e-20, 0.0, 0.5, 0.0, 8.940991045077148e-07, 1.0, 0.5, 0.0]

# show data from VelocityGlyphs
VelocityGlyphsDisplay = Show(VelocityGlyphs, renderView1, 'AMRRepresentation')

# get color transfer function/color map for 'Co-Rotating Velocity'
velocity_corotatingLUT = GetColorTransferFunction('Co-Rotating Velocity')
velocity_corotatingLUT.Set(
    RGBPoints=GenerateRGBPoints(
        range_min=36774534.81499946,
        range_max=8616573489.715328,
    ),
    ScalarRangeInitialized=1.0,
)

# get color transfer function/color map for 'Velocitycms'
velocitycmsLUT = GetColorTransferFunction('Velocity [cm/s]')
velocitycmsLUT.Set(
    RGBPoints=GenerateRGBPoints(
        range_min=1086796.4471111032,
        range_max=83595073.42514606,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
VelocityGlyphsDisplay.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity [cm/s]'],
    LookupTable=velocitycmsLUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
VelocityGlyphsDisplay.ScaleTransferFunction.Points = [3.104755277450931e-18, 0.0, 0.5, 0.0, 1.0050345896692121e-13, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
VelocityGlyphsDisplay.OpacityTransferFunction.Points = [3.104755277450931e-18, 0.0, 0.5, 0.0, 1.0050345896692121e-13, 1.0, 0.5, 0.0]

# show data from coRotatingVelocityGlyphs
coRotatingVelocityGlyphsDisplay = Show(coRotatingVelocityGlyphs, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
coRotatingVelocityGlyphsDisplay.Set(
    Representation='Surface',
    AmbientColor=[0.0, 0.0, 0.0],
    ColorArrayName=['POINTS', ''],
    DiffuseColor=[0.0, 0.0, 0.0],
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
coRotatingVelocityGlyphsDisplay.ScaleTransferFunction.Points = [3.104755277450931e-18, 0.0, 0.5, 0.0, 9.961072422757231e-14, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
coRotatingVelocityGlyphsDisplay.OpacityTransferFunction.Points = [3.104755277450931e-18, 0.0, 0.5, 0.0, 9.961072422757231e-14, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for densitygrcm3LUT in view renderView1
densitygrcm3LUTColorBar = GetScalarBar(densitygrcm3LUT, renderView1)
densitygrcm3LUTColorBar.Set(
    WindowLocation='Upper Right Corner',
    Title='Density [gr/cm^3]',
    ComponentTitle='',
    AllowOverlappingLabels=1,
)

# set color bar visibility
densitygrcm3LUTColorBar.Visibility = 1

# show color legend
axisAlignedSlice1Display.SetScalarBarVisibility(renderView1, True)
VelocityGlyphsDisplay.SetScalarBarVisibility(renderView1, True)

timeAnimationCue1 = GetTimeTrack()

timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper
timeKeeper1.SuppressedTimeSources = reader

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.Set(
    ViewModules=renderView1,
    Cues=timeAnimationCue1,
    AnimationTime=682111.484339536,
    StartTime=682111.484339536,
    EndTime=682112.484339536,
    PlayMode='Snap To TimeSteps',
)

if __name__ == '__main__':
  img_fname = "Co-Rotating-Velocity.png"
  print(f"rendering image {img_fname}")
  SaveScreenshot(img_fname, ImageResolution=[1920,1080])
