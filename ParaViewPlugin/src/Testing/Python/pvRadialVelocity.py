# state file generated using paraview version 6.2.0
import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 1
import h5py
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

renderView1 = CreateView('RenderView')

import sys
import os
    
renderView1.Set(
    OrientationAxesVisibility=0,
    CenterOfRotation=[50000000188416.0, 50000000188416.0, 50000000188416.0],
    CameraPosition=[50000000188416.0, 50000000188416.0, 116792085367294.34],
    CameraFocalPoint=[50000000188416.0, 50000000188416.0, -211323526086314.56],
)

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.Set(
    OrientationAxesVisibility=0,
    CenterOfRotation=[50000000188416.0, 50000000188416.0, 50000000188416.0],
    CameraPosition=[50000000188416.0, 50000000188416.0, 116792085367294.34],
    CameraFocalPoint=[50000000188416.0, 50000000188416.0, -211323526086314.56],
)
executable = os.path.basename(sys.executable)

if executable.startswith("pvbatch"):
    renderView1.ViewSize=[1920//2,1080]
    renderView2.ViewSize=[1920//2,1080]
    
SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitHorizontal(0, 0.5)
layout1.AssignView(1, renderView1)
layout1.AssignView(2, renderView2)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

directory_JeanLaptop = '/local/data/Walder/CygX-1_M1=14.8_M2=19.2_V=750_ML=1-6_Gamma=1.01_RBH=20_HR/DATA/'
directory_r740gpu02 = '/Xnfs/highenergy/rwalder/Simulations/CygX-1_M1=14.8_M2=19.2_V=750_ML=1-6_Gamma=1.01_RBH=20_HR/'

fname = 'DATA.CX-1_M1=14.8_M2=19.2_V=750_ML=1-6_G=1.01_RBH=20RG_.NT=0000008600.Time=189.4754123165H.amr5'

# create a new directory name for Rolf and Doris's laptops, or use Jean's laptop, or use r740gpu02
# watch for the last "/" in the directory name

#filename = directory_r740gpu02 + fname
filename = directory_JeanLaptop + fname

##############################################################################
# Read my FieldData to use later inside the ProgrammableFilter
# These numerical values are visible from inside the ProgrammableFilter (but not from the Python shell)
myfile = h5py.File(filename, "r")

atts_g = myfile["APR_StellarSystems"]["Stars"]["Stellar Orbit"]["Fixed Binary Orbit"]

CM_x, CM_y, CM_z = atts_g.attrs["Center of Mass - x"], atts_g.attrs["Center of Mass - y"], atts_g.attrs["Center of Mass - z"]

Omega = atts_g.attrs["Omega"]

dset = myfile["Thermodynamics"]["Polytropic EOS"][0]
GasConstant, MeanParticleMass = dset['Specific Gasconstant'], dset['Mean Particle Mass']

myfile.close()
##############################################################################

reader = AMAZEReader(registrationName='reader', FileNames=filename)
reader.Set(
    PointArrayStatus=['Density', 'Velocity'],
    Level=100,
)
reader.UpdatePipeline()

Global_bounds = reader.GetDataInformation().GetBounds() # an array of [xmin, xmax, ymin, ymax, zmin, zmax]
center_of_box = [(Global_bounds[1] - Global_bounds[0]) * 0.5,
                 (Global_bounds[3] - Global_bounds[2]) * 0.5,
                 (Global_bounds[5] - Global_bounds[4]) * 0.5
                ]
                
# create a new 'Axis-Aligned Slice'
axisAlignedSlice1 = AxisAlignedSlice(registrationName='AxisAlignedSlice1', Input=reader)
axisAlignedSlice1.Level = 18

# init the 'Axis Aligned Plane' selected for 'CutFunction'
axisAlignedSlice1.CutFunction.Set(
    Origin=center_of_box,
    Normal=[0.0, 0.0, 1.0],
)

# create a new 'Programmable Filter'
RadialVelocity = ProgrammableFilter(registrationName='Radial Velocity', Input=axisAlignedSlice1)
RadialVelocity.Set(
    RequestInformationScript='',
    RequestUpdateExtentScript='',
    CopyArrays=1,
    PythonPath='',
    Script="""import numpy as np
from vtkmodules.vtkCommonDataModel import vtkOverlappingAMR
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk

######################################################
# define my computation for a single cartesian grid
def compute_radial_velocity(grid, output_grid):
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
        # Radial Velocity (with respect to CM or BH respectively):
        # Radius_CM= Coords-Pos_CM (Radius_BH=Coords-Pos_BH)
        # dot(normalize(Radius_CM), Velocity)*normalize(Radius_CM)
        # --------------------------------------------------------------

        Radius_CM = coords - Pos_CM
        r_norm = np.linalg.norm(Radius_CM, axis=1)
        radial_velocity = np.zeros_like(r_norm)
        np.divide(np.sum(velocity * Radius_CM, axis=1),
                  r_norm,
                  out=radial_velocity,
                  where=r_norm > 0)

        radial_vtk = numpy_to_vtk(radial_velocity, deep=True)

        radial_vtk.SetName("Radial Velocity")
        output_grid.GetPointData().AddArray(radial_vtk)
        
######################################################

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------

Omega_v = np.array([0.0, 0.0, Omega])
Pos_CM = np.array([CM_x, CM_y, CM_z])

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
        
        compute_radial_velocity(grid, output_grid)
        # Put block into AMR output
        output_amr.SetDataSet(level, block, output_grid)""",
)

# show data from reader_1
reader_1Display = Show(OutputPort(reader, 1), renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
reader_1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
reader_1Display.ScaleTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 19.2, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
reader_1Display.OpacityTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 19.2, 1.0, 0.5, 0.0]

# show data from RadialVelocity
RadialVelocityDisplay = Show(RadialVelocity, renderView1, 'AMRRepresentation')

# get color transfer function/color map for 'RadialVelocity'
radialVelocityLUT = GetColorTransferFunction('RadialVelocity')
radialVelocityLUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=GenerateRGBPoints(
        range_min=-1000000.0,
        range_max=1000000.0,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
RadialVelocityDisplay.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Radial Velocity'],
    LookupTable=radialVelocityLUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
RadialVelocityDisplay.ScaleTransferFunction.Points = [-6660371613.680384, 0.0, 0.5, 0.0, 7294142954.603252, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
RadialVelocityDisplay.OpacityTransferFunction.Points = [-6660371613.680384, 0.0, 0.5, 0.0, 7294142954.603252, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for radialVelocityLUT in view renderView1
radialVelocityLUTColorBar = GetScalarBar(radialVelocityLUT, renderView1)
radialVelocityLUTColorBar.Set(
    WindowLocation='Upper Right Corner',
    Title='Radial Velocity',
    ComponentTitle='',
)

# set color bar visibility
radialVelocityLUTColorBar.Visibility = 1

# show color legend
RadialVelocityDisplay.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from RadialVelocity
RadialVelocityDisplay_1 = Show(RadialVelocity, renderView2, 'AMRRepresentation')

# get color transfer function/color map for 'Velocitycms'
velocitycmsLUT = GetColorTransferFunction('Velocitycms')
velocitycmsLUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=[
        # scalar, red, green, blue
        1000000.0, 0.0564, 0.0564, 0.47,
        3271719.0293607833, 0.243, 0.46035, 0.81,
        7860913.184908469, 0.356814, 0.745025, 0.954368,
        19787321.037596498, 0.6882, 0.93, 0.91791,
        31622776.60168379, 0.899496, 0.944646, 0.768657,
        58167178.75430428, 0.957108, 0.833819, 0.508916,
        131347859.58042039, 0.927521, 0.621439, 0.315357,
        349076007.3789005, 0.8, 0.352, 0.16,
        1000000000.0, 0.59, 0.0767, 0.119475,
    ],
    UseLogScale=1,
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
RadialVelocityDisplay_1.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity [cm/s]'],
    LookupTable=velocitycmsLUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
RadialVelocityDisplay_1.ScaleTransferFunction.Points = [-6660371613.680384, 0.0, 0.5, 0.0, 7294142954.603252, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
RadialVelocityDisplay_1.OpacityTransferFunction.Points = [-6660371613.680384, 0.0, 0.5, 0.0, 7294142954.603252, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for velocitycmsLUT in view renderView2
velocitycmsLUTColorBar = GetScalarBar(velocitycmsLUT, renderView2)
velocitycmsLUTColorBar.Set(
    Title='Velocity [cm/s]',
)

# set color bar visibility
velocitycmsLUTColorBar.Visibility = 1

# show color legend
RadialVelocityDisplay_1.SetScalarBarVisibility(renderView2, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'Velocitycms'
velocitycmsPWF = GetOpacityTransferFunction('Velocitycms')
velocitycmsPWF.Set(
    Points=[1000000.0, 0.0, 0.5, 0.0, 1000000000.0, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'RadialVelocity'
radialVelocityPWF = GetOpacityTransferFunction('RadialVelocity')
radialVelocityPWF.Set(
    Points=[-1000000.0, 0.0, 0.5, 0.0, 1000000.0, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper
timeKeeper1.SuppressedTimeSources = reader

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.Set(
    ViewModules=[renderView1, renderView2],
    Cues=timeAnimationCue1,
    AnimationTime=682111.484339536,
    StartTime=682111.484339536,
    EndTime=682112.484339536,
    PlayMode='Snap To TimeSteps',
)

# initialize the animation scene

# ----------------------------------------------------------------
# restore active source
SetActiveSource(RadialVelocity)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup camera links
AddCameraLink(renderView2, renderView1, 'CameraLink0')

if executable.startswith("pvbatch"):
  img_fname = "ThermalPressure.png"
  print(f"rendering image {img_fname}")
  SaveScreenshot(img_fname, layout1, SaveAllViews=1)
