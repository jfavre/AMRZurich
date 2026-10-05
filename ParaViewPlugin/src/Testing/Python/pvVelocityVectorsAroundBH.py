##############################################################################
#
# The purpose of this script is to show that using the default scenario to show
# velocity glyphs creates an over-crowded image (the left-side), because ParaView
# replicates the sequence Mask-Points + Glyphs for all levels, such that finer
# levels of resolution get too much coverage
#
# The trick I use here to create the right-side images is to "Merge" all levels
# of the iso-contour level into a single object, and then apply the Mask-Points +
# Glyphs operation.
# IMHO, it makes a better plot.
#
##############################################################################
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
renderView1.Set(
    OrientationAxesVisibility=0,
    CenterOfRotation=[48586775592960.0, 50911124652032.0, 49998796423168.0],
    CameraPosition=[49026158811564.02, 50841115317235.516, 50469980164650.04],
    CameraFocalPoint=[48500874450982.086, 50924811751490.83, 49906678154264.46],
    CameraViewUp=[-0.7075214121952297, 0.17228693346660492, 0.6853689983081681],
)

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.Set(
    CenterOfRotation=[48586775592960.0, 50911124652032.0, 49998796423168.0],
    CameraPosition=[49026158811564.02, 50841115317235.516, 50469980164650.04],
    CameraFocalPoint=[48500874450982.086, 50924811751490.83, 49906678154264.46],
    CameraViewUp=[-0.7075214121952297, 0.17228693346660492, 0.6853689983081681],
)

if __name__ == '__main__':
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
layout1.SetSize(1920,1080)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView2)
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
    PointArrayStatus=['Density', 'Velocity'],
    Level=100,
)
reader.UpdatePipeline()

Global_bounds = reader.GetDataInformation().GetBounds() # an array of [xmin, xmax, ymin, ymax, zmin, zmax]
center_of_box = [(Global_bounds[1] - Global_bounds[0]) * 0.5,
                 (Global_bounds[3] - Global_bounds[2]) * 0.5,
                 (Global_bounds[5] - Global_bounds[4]) * 0.5
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
    RandomSamplingMode='Uniform Spatial Distribution (Bounds Based)',
    GenerateVertices=1,
)

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=maskPoints1,
    GlyphType='Arrow')
glyph1.Set(
    OrientationArray=['POINTS', 'Velocity [cm/s]'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=21500000000.0,
    GlyphMode='All Points',
)

# create a new 'Merge Blocks'
mergeBlocks1 = MergeBlocks(registrationName='MergeBlocks1', Input=contour1)
mergeBlocks1.MergePoints = 0

# create a new 'Mask Points'
maskPoints2 = MaskPoints(registrationName='MaskPoints2', Input=mergeBlocks1)
maskPoints2.Set(
    MaximumNumberofPoints=5000,
    RandomSampling=1,
    RandomSamplingMode='Uniform Spatial Distribution (Bounds Based)',
    GenerateVertices=1,
)

# create a new 'Glyph'
glyph2 = Glyph(registrationName='Glyph2', Input=maskPoints2,
    GlyphType='Arrow')
glyph2.Set(
    OrientationArray=['POINTS', 'Velocity [cm/s]'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=21500000000.0,
    GlyphMode='All Points',
)

# create a new 'Slice'
midboxslice = Slice(registrationName='mid-box-slice', Input=reader)
midboxslice.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
midboxslice.SliceType.Set(
    Origin=center_of_box,
    Normal=[0.0, 0.0, 1.0],
)


# show data from glyph1
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'Velocitycms'
velocitycmsLUT = GetColorTransferFunction('Velocitycms')
velocitycmsLUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=GenerateRGBPoints(
        preset_name='Fast',
        range_min=49097944.71770766,
        range_max=535998529.4941732,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
glyph1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity [cm/s]'],
    LookupTable=velocitycmsLUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [1.000000013351432e-10, 0.0, 0.5, 0.0, 1.000142121898584e-10, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [1.000000013351432e-10, 0.0, 0.5, 0.0, 1.000142121898584e-10, 1.0, 0.5, 0.0]

# show data from midboxslice
midboxsliceDisplay = Show(midboxslice, renderView1, 'GeometryRepresentation')

# get color transfer function/color map for 'Densitygrcm3'
densitygrcm3LUT = GetColorTransferFunction('Densitygrcm3')
densitygrcm3LUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=[
        # scalar, red, green, blue
        2.0000000000000002e-16, 0.301961, 0.047059, 0.090196,
        2.9342424284933736e-16, 0.396078431372549, 0.0392156862745098, 0.058823529411764705,
        4.3048893145853273e-16, 0.49411764705882355, 0.054901960784313725, 0.03529411764705882,
        6.315794438412012e-16, 0.5882352941176471, 0.11372549019607843, 0.023529411764705882,
        9.266036005415464e-16, 0.6627450980392157, 0.16862745098039217, 0.01568627450980392,
        1.3594397995518598e-15, 0.7411764705882353, 0.22745098039215686, 0.00392156862745098,
        1.994462969413797e-15, 0.788235294117647, 0.2901960784313726, 0.0,
        2.92611893345641e-15, 0.8627450980392157, 0.3803921568627451, 0.011764705882352941,
        4.292971162682788e-15, 0.9019607843137255, 0.4588235294117647, 0.027450980392156862,
        6.298309064921156e-15, 0.9176470588235294, 0.5215686274509804, 0.047058823529411764,
        9.24038284302804e-15, 0.9254901960784314, 0.5803921568627451, 0.0784313725490196,
        1.3556761696767493e-14, 0.9372549019607843, 0.6431372549019608, 0.12156862745098039,
        1.9889412681814496e-14, 0.9450980392156862, 0.7098039215686275, 0.1843137254901961,
        2.918017928439701e-14, 0.9529411764705882, 0.7686274509803922, 0.24705882352941178,
        4.2810860063660557e-14, 0.9647058823529412, 0.8274509803921568, 0.3254901960784314,
        6.280872099954241e-14, 0.9686274509803922, 0.8784313725490196, 0.4235294117647059,
        9.214800701813002e-14, 0.9725490196078431, 0.9176470588235294, 0.5137254901960784,
        1.3519229594685056e-13, 0.9803921568627451, 0.9490196078431372, 0.596078431372549,
        1.9834348538634e-13, 0.9803921568627451, 0.9725490196078431, 0.6705882352941176,
        2.909939351179271e-13, 0.9882352941176471, 0.9882352941176471, 0.7568627450980392,
        4.204276735542484e-13, 0.984313725490196, 0.9882352941176471, 0.8549019607843137,
        4.2692337542863475e-13, 0.9882352941176471, 0.9882352941176471, 0.8588235294117647,
        4.270232273592505e-13, 0.9529411764705882, 0.9529411764705882, 0.8941176470588236,
        4.270232273592505e-13, 0.9529411764705882, 0.9529411764705882, 0.8941176470588236,
        6.358013135809509e-13, 0.8901960784313725, 0.8901960784313725, 0.807843137254902,
        9.466541500590945e-13, 0.8274509803921568, 0.8235294117647058, 0.7372549019607844,
        1.4094876192954062e-12, 0.7764705882352941, 0.7647058823529411, 0.6784313725490196,
        2.098607341258701e-12, 0.7254901960784313, 0.7137254901960784, 0.6274509803921569,
        3.1246480724580894e-12, 0.6784313725490196, 0.6627450980392157, 0.5803921568627451,
        4.6523355678628835e-12, 0.6313725490196078, 0.6078431372549019, 0.5333333333333333,
        6.926932484583816e-12, 0.5803921568627451, 0.5568627450980392, 0.48627450980392156,
        1.0313614086101644e-11, 0.5372549019607843, 0.5058823529411764, 0.44313725490196076,
        1.535609532123586e-11, 0.4980392156862745, 0.4588235294117647, 0.40784313725490196,
        2.286392156486181e-11, 0.4627450980392157, 0.4196078431372549, 0.37254901960784315,
        3.404243711623974e-11, 0.43137254901960786, 0.38823529411764707, 0.34509803921568627,
        5.068629725331795e-11, 0.403921568627451, 0.3568627450980392, 0.3176470588235294,
        7.546759124440429e-11, 0.37254901960784315, 0.3215686274509804, 0.29411764705882354,
        1.1236483303896581e-10, 0.34509803921568627, 0.29411764705882354, 0.26666666666666666,
        1.6730169196715729e-10, 0.3176470588235294, 0.2627450980392157, 0.23921568627450981,
        2.490980084967267e-10, 0.28627450980392155, 0.23137254901960785, 0.21176470588235294,
        3.708857759144254e-10, 0.2549019607843137, 0.2, 0.1843137254901961,
        5.52217416774143e-10, 0.23137254901960785, 0.17254901960784313, 0.16470588235294117,
        8.222048274481849e-10, 0.2, 0.1450980392156863, 0.13725490196078433,
        1.2244794719879247e-09, 0.14902, 0.196078, 0.278431,
        2.660112310971431e-09, 0.2, 0.2549019607843137, 0.34509803921568627,
        5.778943354186021e-09, 0.24705882352941178, 0.3176470588235294, 0.41568627450980394,
        1.2554427177059686e-08, 0.3058823529411765, 0.38823529411764707, 0.49411764705882355,
        2.7273782088541532e-08, 0.37254901960784315, 0.4588235294117647, 0.5686274509803921,
        5.925074708087675e-08, 0.44313725490196076, 0.5333333333333333, 0.6431372549019608,
        1.2871889268034198e-07, 0.5176470588235295, 0.615686274509804, 0.7254901960784313,
        2.796345050339624e-07, 0.6, 0.6980392156862745, 0.8,
        6.074901265642326e-07, 0.6862745098039216, 0.7843137254901961, 0.8705882352941177,
        1.3197378979686622e-06, 0.7607843137254902, 0.8588235294117647, 0.9294117647058824,
        1.945189516942253e-06, 0.807843137254902, 0.9019607843137255, 0.9607843137254902,
        2.8670558469571843e-06, 0.8901960784313725, 0.9568627450980393, 0.984313725490196,
    ],
    UseLogScale=1,
    NanColor=[0.25, 0.0, 0.0],
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
midboxsliceDisplay.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Density [gr/cm^3]'],
    LookupTable=densitygrcm3LUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
midboxsliceDisplay.ScaleTransferFunction.Points = [9.999996163704499e-20, 0.0, 0.5, 0.0, 8.940991045077148e-07, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
midboxsliceDisplay.OpacityTransferFunction.Points = [9.999996163704499e-20, 0.0, 0.5, 0.0, 8.940991045077148e-07, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for velocitycmsLUT in view renderView1
velocitycmsLUTColorBar = GetScalarBar(velocitycmsLUT, renderView1)
velocitycmsLUTColorBar.Set(
    Orientation='Horizontal',
    WindowLocation='Any Location',
    Position=[0.40, 0.04],
    Title='Velocity [cm/s]',
    ComponentTitle='Magnitude',
    ScalarBarLength=0.33,
    AllowOverlappingLabels=1,
)

# set color bar visibility
velocitycmsLUTColorBar.Visibility = 1

# show color legend
glyph1Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from midboxslice
midboxsliceDisplay_1 = Show(midboxslice, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
midboxsliceDisplay_1.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Density [gr/cm^3]'],
    LookupTable=densitygrcm3LUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
midboxsliceDisplay_1.ScaleTransferFunction.Points = [9.999996163704499e-20, 0.0, 0.5, 0.0, 8.940991045077148e-07, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
midboxsliceDisplay_1.OpacityTransferFunction.Points = [9.999996163704499e-20, 0.0, 0.5, 0.0, 8.940991045077148e-07, 1.0, 0.5, 0.0]

# show data from glyph2
glyph2Display = Show(glyph2, renderView2, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph2Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity [cm/s]'],
    LookupTable=velocitycmsLUT,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
glyph2Display.ScaleTransferFunction.Points = [1.000000013351432e-10, 0.0, 0.5, 0.0, 1.000142121898584e-10, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
glyph2Display.OpacityTransferFunction.Points = [1.000000013351432e-10, 0.0, 0.5, 0.0, 1.000142121898584e-10, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for velocitycmsLUT in view renderView2
velocitycmsLUTColorBar_1 = GetScalarBar(velocitycmsLUT, renderView2)
velocitycmsLUTColorBar_1.Set(
    Orientation='Horizontal',
    WindowLocation='Any Location',
    Position=[0.40, 0.04],
    Title='Velocity [cm/s]',
    ComponentTitle='Magnitude',
    ScalarBarLength=0.33,
)

# set color bar visibility
velocitycmsLUTColorBar_1.Visibility = 1

# show color legend
glyph2Display.SetScalarBarVisibility(renderView2, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'Velocitycms'
velocitycmsPWF = GetOpacityTransferFunction('Velocitycms')
velocitycmsPWF.Set(
    Points=[49097944.71770766, 0.0, 0.5, 0.0, 535998529.4941732, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'Densitygrcm3'
densitygrcm3PWF = GetOpacityTransferFunction('Densitygrcm3')
densitygrcm3PWF.Set(
    Points=[2e-16, 0.0, 0.5, 0.0, 2.233550340855008e-16, 0.0, 0.5, 0.0, 2.259500339775932e-16, 1.0, 0.5, 0.0, 2.8670558469571783e-06, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

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
    ViewModules=[renderView1, renderView2],
    Cues=timeAnimationCue1,
    AnimationTime=3954516.40220226,
    StartTime=3954516.40220226,
    EndTime=3954517.40220226,
    PlayMode='Snap To TimeSteps',
)

# ----------------------------------------------------------------
# restore active source
SetActiveSource(contour1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup camera links
AddCameraLink(renderView1, renderView2, 'CameraLink0')

if __name__ == '__main__':
  img_fname = "VelocityVectorsAroundBH.png"
  print(f"rendering image {img_fname}")
  SaveScreenshot(img_fname, viewOrLayout=layout1, location=16, SaveAllViews=1, ImageResolution=[1920, 1080])
