# state file generated using paraview version 6.2.0
import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

renderView1 = GetRenderView()

import sys
import os

executable = os.path.basename(sys.executable)

if executable.startswith("pvbatch"):
    renderView1.ViewSize=[1920,1080]
    
renderView1.Set(
    OrientationAxesVisibility=0,
    CenterOfRotation=[50000000188416.0, 50000000188416.0, 50000000188416.0],
    CameraPosition=[50000000188416.0, 50000000188416.0, 70513427149559.25],
    CameraFocalPoint=[50000000188416.0, 50000000188416.0, -202726256550417.97],
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
    PointArrayStatus=['Density', 'Temperature'],
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
ThermalPressure = ProgrammableFilter(registrationName='Thermal Pressure', Input=axisAlignedSlice1)
ThermalPressure.Set(
    RequestInformationScript='',
    RequestUpdateExtentScript='',
    CopyArrays=1,
    PythonPath='',
    Script="""import numpy as np
from vtkmodules.vtkCommonDataModel import vtkOverlappingAMR
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk

######################################################
# define my computation for a single cartesian grid
def compute_thermal_pressure_from_density_and_temperature(grid, output_grid):
        # --------------------------------------------------------------
        # derive thermal pressure from density and temperature:
        # MeanParticleMass=0.61
        # GasConstant_Norm= GasConstant/MeanParticleMass
        # Density*GasConstant_Norm*Temperature
        # --------------------------------------------------------------

        density_vtk = grid.GetPointData().GetArray("Density [gr/cm^3]")
        if density_vtk is None:
            raise RuntimeError(
                f"No \'Density\' array in AMR level {level}, block {block}"
            )
        density = vtk_to_numpy(density_vtk)
        
        temperature_vtk = grid.GetPointData().GetArray("Temperature [Kelvin]")
        if temperature_vtk is None:
            raise RuntimeError(
                f"No \'Temperature\' array in AMR level {level}, block {block}"
            )
        temperature = vtk_to_numpy(temperature_vtk)

        thermal_pressure_D = density * Specific_Gasconstant * temperature
        
        # --------------------------------------------------------------
        # Add result to output
        # --------------------------------------------------------------

        thermal_pressure_D_vtk = numpy_to_vtk(thermal_pressure_D, deep=True)
        thermal_pressure_D_vtk.SetName("thermal_pressure_D")

        output_grid.GetPointData().AddArray(thermal_pressure_D_vtk)
        
######################################################

# ----------------------------------------------------------------------
# Input / output of the nested AMR structures
# ----------------------------------------------------------------------

input_amr = self.GetInputDataObject(0, 0)
output_amr = self.GetOutputDataObject(0)

if not isinstance(input_amr, vtkOverlappingAMR):
    raise RuntimeError(
        f"Expected vtkOverlappingAMR, got {input_amr.GetClassName()}"
    )

# Thermodynamics_EOS data are available in that order:
#typedef struct Thermodynamics_EOS
#  {
#  double Adiabatic_Exponent;
#  double Specific_Gasconstant;
#  double Scaled_Specific_Gasconstant;
#  double Scale_Factor;
#  double Mean_Particle_Mass;
#} Polytropic_EOS;
####################

MeanParticleMass          = input_amr.GetFieldData().GetArray("Thermodynamics_EOS").GetValue(4)
Specific_Gasconstant      = input_amr.GetFieldData().GetArray("Thermodynamics_EOS").GetValue(1)

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
        
        compute_thermal_pressure_from_density_and_temperature(grid, output_grid)
        
        # Put block into AMR output
        output_amr.SetDataSet(level, block, output_grid)""",
)

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

# show data from ThermalPressure
thermalPressureDisplay = Show(ThermalPressure, renderView1, 'AMRRepresentation')

# get color transfer function/color map for 'thermal_pressure'
thermal_pressureLUT = GetColorTransferFunction('thermal_pressure_D')
thermal_pressureLUT.Set(
    RGBPoints=[
        # scalar, red, green, blue
        2.200546489449541e-07, 0.0564, 0.0564, 0.47,
        0.00023504030455570138, 0.243, 0.46035, 0.81,
        0.04082435201591318, 0.356814, 0.745025, 0.954368,
        9.324812745307433, 0.6882, 0.93, 0.91791,
        147.0850080586625, 0.899496, 0.944646, 0.768657,
        5305.960214501072, 0.957108, 0.833819, 0.508916,
        639694.5753154006, 0.927521, 0.621439, 0.315357,
        201113272.2522087, 0.8, 0.352, 0.16,
        98311940689.91692, 0.59, 0.0767, 0.119475,
    ],
    UseLogScale=1,
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
thermalPressureDisplay.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'thermal_pressure_D'],
    LookupTable=thermal_pressureLUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
thermalPressureDisplay.ScaleTransferFunction.Points = [-6660371613.680384, 0.0, 0.5, 0.0, 7294142954.603252, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
thermalPressureDisplay.OpacityTransferFunction.Points = [-6660371613.680384, 0.0, 0.5, 0.0, 7294142954.603252, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for thermal_pressureLUT in view renderView1
thermal_pressureLUTColorBar = GetScalarBar(thermal_pressureLUT, renderView1)
thermal_pressureLUTColorBar.Set(
    Title='thermal_pressure_D',
    ComponentTitle='',
)

# set color bar visibility
thermal_pressureLUTColorBar.Visibility = 1

# show color legend
thermalPressureDisplay.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'thermal_pressure'
thermal_pressurePWF = GetOpacityTransferFunction('thermal_pressure_D')
thermal_pressurePWF.Set(
    Points=[2.200546489449543e-07, 0.0, 0.5, 0.0, 98311940689.91743, 1.0, 0.5, 0.0],
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
    ViewModules=[renderView1],
    Cues=timeAnimationCue1,
    AnimationTime=682111.484339536,
    StartTime=682111.484339536,
    EndTime=682112.484339536,
    PlayMode='Snap To TimeSteps',
)

# initialize the animation scene

# ----------------------------------------------------------------
# restore active source
SetActiveSource(ThermalPressure)

if executable.startswith("pvbatch"):
  img_fname = "ThermalPressure_D.png"
  print(f"rendering image {img_fname}")
  SaveScreenshot(img_fname, ImageResolution=[1920,1080])

