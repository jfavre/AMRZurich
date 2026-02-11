# state file generated using paraview version 6.
import paraview
import paraview.util

#### import the simple module from the paraview
from paraview.simple import *

Version = GetParaViewSourceVersion().split(" ")[-1][0:3]
#returns "6.1"

LoadPlugin("/home/jfavre/Projects/AMRZurich/ParaViewPlugin/build" + Version + "/lib/paraview-" + Version + "/plugins/pvAMAZEReader/pvAMAZEReader.so", ns=globals())
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# Create a new 'Render View'
renderView1 = GetRenderView()
renderView1.Set(
    ViewSize=[1024,1024],
    CenterOfRotation=[100000000376832.0, 100000000376832.0, 100000000376832.0],
    CameraPosition=[100000000376832.0, 100000000376832.0, 134164839639387.86],
    CameraFocalPoint=[100000000376832.0, 100000000376832.0, -535048205872668.75],
)

basename = '/local/data/Walder/DATA.CygX-1_M1/'
fnames = [basename + 'DATA.CygX-1_M1=14.8_M2=19.2_D=2e14_V=1200_ML=2-6_Gamma=1.6666_R=64.COOL.NT=0000016530.Time=2278461.3937737681S.amr5'
,
basename +
'DATA.CygX-1_M1=14.8_M2=19.2_D=2e14_V=1200_ML=2-6_Gamma=1.6666_R=64.COOL.NT=0000016525.Time=2278308.5024622753S.amr5'
,
basename +
'DATA.CygX-1_M1=14.8_M2=19.2_D=2e14_V=1200_ML=2-6_Gamma=1.6666_R=64.COOL.NT=0000016553.Time=2279142.6533078365S.amr5']

#fnames = paraview.util.Glob(path = "/local/data/Walder/DATA.CygX-1_M1/*.amr5")
print("fnames = ", fnames)

reader = AMAZEReader(registrationName='DATA.CygX-1_M1', FileNames=fnames)
reader.PointArrayStatus=['Density']
reader.LevelSet=[0, 16]
reader.UpdatePipelineInformation()

print("TimestepValues = ", reader.TimestepValues)
annotateTime1 = AnnotateTime(registrationName='AnnotateTime1')
annotateTime1Display = Show(annotateTime1, renderView1, 'TextSourceRepresentation')

# show data from reader
Display = Show(reader, renderView1, 'AMRRepresentation')
Display.Set(
    Representation='Outline',
    ColorArrayName=['POINTS', ''],
    Assembly='Hierarchy',
)

slice1 = Slice(registrationName='Slice1', Input=reader)
slice1.SliceOffsetValues = [0.0]
slice1.SliceType.Set(
    Origin=[100000000000000.0, 100000000000000.0, 100000000000000.0],
    Normal=[0.0, 0.0, 1.0],
)

slice1Display = Show(slice1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'Densitygrcm3'
densitygrcm3LUT = GetColorTransferFunction('Densitygrcm3')
densitygrcm3LUT.Set(
    RGBPoints=[
        # scalar, red, green, blue
        9.9999571622828e-20, 0.0564, 0.0564, 0.47,
        1.561072104220104e-17, 0.243, 0.46035, 0.81,
        6.539528989823782e-16, 0.356814, 0.745025, 0.954368,
        3.340513337178177e-14, 0.6882, 0.93, 0.91791,
        2.4626011470992334e-13, 0.899496, 0.944646, 0.768657,
        3.304978765295026e-12, 0.957108, 0.833819, 0.508916,
        1.0628000776895798e-10, 0.927521, 0.621439, 0.315357,
        6.842327307194751e-09, 0.8, 0.352, 0.16,
        6.064430388329856e-07, 0.59, 0.0767, 0.119475,
    ],
    UseLogScale=1,
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
slice1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Density [gr/cm^3]'],
    LookupTable=densitygrcm3LUT,
    Assembly='Hierarchy',
)

densitygrcm3LUTColorBar = GetScalarBar(densitygrcm3LUT, renderView1)
densitygrcm3LUTColorBar.Set(
    Title='Density [gr/cm^3]',
    ComponentTitle='',
)

densitygrcm3LUTColorBar.Visibility = 1

timeAnimationCue1 = GetTimeTrack()

timeKeeper1 = GetTimeKeeper()

animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.Set(
    ViewModules=renderView1,
    Cues=timeAnimationCue1,
    AnimationTime=reader.TimestepValues[0],
    StartTime=reader.TimestepValues[0],
    EndTime=reader.TimestepValues[-1],
    PlayMode='Snap To TimeSteps',
)

SaveAnimation("animation.png")
