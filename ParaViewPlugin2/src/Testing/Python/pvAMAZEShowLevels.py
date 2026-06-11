# state file generated using paraview version 6.1.1-617-g6391ede1d0
import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 1

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


if __name__ == '__main__':
  LoadPlugin("pvAMAZEReader2", ns=globals())
  
# Create a new 'Render View'
renderView3 = CreateView('RenderView')
renderView3.Set(
    ViewSize=[1095, 1167],
    CenterOfRotation=[50000000188416.0, 50000000188416.0, 50000000188416.0],
    CameraPosition=[678853698522586.5, -361671823723627.0, 8862195700943370.0],
    CameraFocalPoint=[678366494102997.5, -361352880967050.1, 8855383830609538.0],
    CameraViewAngle=0.5872725310212507,
)

# Create a new 'Render View'
renderView4 = CreateView('RenderView')
renderView4.Set(
    ViewSize=[1095, 1167],
    CenterOfRotation=[50000000188416.0, 50000000188416.0, 50000000188416.0],
    CameraPosition=[678853698522586.5, -361671823723627.0, 8862195700943370.0],
    CameraFocalPoint=[678366494102997.5, -361352880967050.1, 8855383830609538.0],
    CameraViewAngle=0.5872725310212507,
)

SetActiveView(None)


# create new layout object 'Layout #1'
layout1_1 = CreateLayout(name='Layout #1')
layout1_1.SplitHorizontal(0, 0.5)
layout1_1.AssignView(1, renderView3)
layout1_1.AssignView(2, renderView4)
layout1_1.SetSize(2191, 1167)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView4)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'AMAZE Reader'
reader = AMAZEReader(registrationName='DATA.CX-1_M1=14.8_M2=19.2_V=750_ML=1-6_G=1.01_RBH=20RG_.NT=0000008200.Time=189.2836067610H.amr5', FileNames=['/local/data/Walder/CygX-1_M1=14.8_M2=19.2_V=750_ML=1-6_Gamma=1.01_RBH=20_HR/DATA/DATA.CX-1_M1=14.8_M2=19.2_V=750_ML=1-6_G=1.01_RBH=20RG_.NT=0000008200.Time=189.2836067610H.amr5'])
reader.Set(
    PointArrayStatus=['Density'],
    Level=5,
)

# create a new 'Slice'
slice1 = Slice(registrationName='Slice1', Input=reader)
slice1.SliceOffsetValues = [0.0]

# init the 'Plane' selected for 'SliceType'
slice1.SliceType.Set(
    Origin=[50000000000000.0, 50000000000000.0, 50000000000000.0],
    Normal=[0.0, 0.0, 1.0],
)

# create a new 'Axis-Aligned Slice'
axisAlignedSlice1 = AxisAlignedSlice(registrationName='AxisAlignedSlice1', Input=reader)
axisAlignedSlice1.Level = 18

# init the 'Axis Aligned Plane' selected for 'CutFunction'
axisAlignedSlice1.CutFunction.Set(
    Origin=[50000000000000.0, 50000000000000.0, 50000000000000.0],
    Normal=[0.0, 0.0, 1.0],
)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView3'
# ----------------------------------------------------------------

# show data from axisAlignedSlice1
axisAlignedSlice1Display = Show(axisAlignedSlice1, renderView3, 'AMRRepresentation')

# get color transfer function/color map for 'vtkAMRLevel'
vtkAMRLevelLUT = GetColorTransferFunction('vtkAMRLevel')
vtkAMRLevelLUT.Set(
    InterpretValuesAsCategories=1,
    AnnotationsInitialized=1,
    ScalarRangeInitialized=1.0,
    Annotations=['0', '0', '1', '1', '2', '2', '3', '3', '4', '4', '5', '5', '6', '6', '7', '7', '8', '8', '9', '9', '10', '10', '11', '11', '12', '12', '13', '13', '14', '14', '15', '15', '16', '16', '17', '17', '18', '18'],
    ActiveAnnotatedValues=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18'],
    IndexedColors=[1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.6299992203712463, 0.6299992203712463, 1.0, 0.712463, 0.3712463, 1.0, 0.463, 0.12463, 1.0, 0.5, 1.0, 1.0, 1.0, 0.5, 0.0, 0.0, 1.0, 0.5, 0.0, 0.4, 0.4, 0.5, 1.0, 0.0, 1.0, 0.8, 1.0, 0.0, 1.0, 0.8, 0.3, 0.2, 1.0, 0.7, 0.3, 1.0],
)

# trace defaults for the display properties.
axisAlignedSlice1Display.Set(
    Representation='Surface',
    ColorArrayName=['CELLS', 'vtkAMRLevel'],
    LookupTable=vtkAMRLevelLUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
axisAlignedSlice1Display.ScaleTransferFunction.Points = [9.999984656091761e-20, 0.0, 0.5, 0.0, 2.8670558469571787e-06, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
axisAlignedSlice1Display.OpacityTransferFunction.Points = [9.999984656091761e-20, 0.0, 0.5, 0.0, 2.8670558469571787e-06, 1.0, 0.5, 0.0]

# find source
reader_1 = reader

# show data from reader_1
reader_1Display = Show(OutputPort(reader_1, 1), renderView3, 'GeometryRepresentation')

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

# setup the color legend parameters for each legend in this view

# get color legend/bar for vtkAMRLevelLUT in view renderView3
vtkAMRLevelLUTColorBar = GetScalarBar(vtkAMRLevelLUT, renderView3)
vtkAMRLevelLUTColorBar.Set(
    Title='vtkAMRLevel',
    ComponentTitle='',
)

# set color bar visibility
vtkAMRLevelLUTColorBar.Visibility = 1

# show color legend
axisAlignedSlice1Display.SetScalarBarVisibility(renderView3, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView4'
# ----------------------------------------------------------------

# show data from reader_1
reader_1Display_1 = Show(OutputPort(reader_1, 1), renderView4, 'GeometryRepresentation')

# trace defaults for the display properties.
reader_1Display_1.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
reader_1Display_1.ScaleTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 19.2, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
reader_1Display_1.OpacityTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 19.2, 1.0, 0.5, 0.0]

clip1 = Clip(registrationName='Clip1', Input=slice1)
clip1.Set(
    ClipType='Scalar',
    Scalars=['POINTS', 'Density [gr/cm^3]'],
    Value=5e-17,
    Invert=0,
)

slice1Display = Show(clip1, renderView4, 'GeometryRepresentation')

# get color transfer function/color map for 'Densitygrcm3'
densitygrcm3LUT = GetColorTransferFunction('Densitygrcm3')
densitygrcm3LUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=[
        # scalar, red, green, blue
        1e-20, 0.301961, 0.047059, 0.090196,
        1.5728342357220066e-20, 0.396078431372549, 0.0392156862745098, 0.058823529411764705,
        2.4738075330592287e-20, 0.49411764705882355, 0.054901960784313725, 0.03529411764705882,
        3.890889180582618e-20, 0.5882352941176471, 0.11372549019607843, 0.023529411764705882,
        6.119723710620686e-20, 0.6627450980392157, 0.16862745098039217, 0.01568627450980392,
        9.625310965223929e-20, 0.7411764705882353, 0.22745098039215686, 0.00392156862745098,
        1.513901861557463e-19, 0.788235294117647, 0.2901960784313726, 0.0,
        2.381116677380855e-19, 0.8627450980392157, 0.3803921568627451, 0.011764705882352941,
        3.7451018294333025e-19, 0.9019607843137255, 0.4588235294117647, 0.027450980392156862,
        5.890424373597817e-19, 0.9176470588235294, 0.5215686274509804, 0.047058823529411764,
        9.264661117726001e-19, 0.9254901960784314, 0.5803921568627451, 0.0784313725490196,
        1.4571776188321966e-18, 0.9372549019607843, 0.6431372549019608, 0.12156862745098039,
        2.2918988464271515e-18, 0.9450980392156862, 0.7098039215686275, 0.1843137254901961,
        3.604776970472456e-18, 0.9529411764705882, 0.7686274509803922, 0.24705882352941178,
        5.669716631301337e-18, 0.9647058823529412, 0.8274509803921568, 0.3254901960784314,
        8.917524424553188e-18, 0.9686274509803922, 0.8784313725490196, 0.4235294117647059,
        1.402578771282444e-17, 0.9725490196078431, 0.9176470588235294, 0.5137254901960784,
        2.2060239097699517e-17, 0.9803921568627451, 0.9490196078431372, 0.596078431372549,
        3.469709930107495e-17, 0.9803921568627451, 0.9725490196078431, 0.6705882352941176,
        5.457278566097723e-17, 0.9882352941176471, 0.9882352941176471, 0.7568627450980392,
        8.429304810030729e-17, 0.984313725490196, 0.9882352941176471, 0.8549019607843137,
        8.583394562630401e-17, 0.9882352941176471, 0.9882352941176471, 0.8588235294117647,
        8.585766569884266e-17, 0.9529411764705882, 0.9529411764705882, 0.8941176470588236,
        8.585766569884266e-17, 0.9529411764705882, 0.9529411764705882, 0.8941176470588236,
        1.3741318259107095e-16, 0.8901960784313725, 0.8901960784313725, 0.807843137254902,
        2.1992657960256764e-16, 0.8274509803921568, 0.8235294117647058, 0.7372549019607844,
        3.5198733850464677e-16, 0.7764705882352941, 0.7647058823529411, 0.6784313725490196,
        5.633474893824875e-16, 0.7254901960784313, 0.7137254901960784, 0.6274509803921569,
        9.016244593961862e-16, 0.6784313725490196, 0.6627450980392157, 0.5803921568627451,
        1.4430288251973039e-15, 0.6313725490196078, 0.6078431372549019, 0.5333333333333333,
        2.309533829356004e-15, 0.5803921568627451, 0.5568627450980392, 0.48627450980392156,
        3.696354789177927e-15, 0.5372549019607843, 0.5058823529411764, 0.44313725490196076,
        5.915929246764168e-15, 0.4980392156862745, 0.4588235294117647, 0.40784313725490196,
        9.468306168873898e-15, 0.4627450980392157, 0.4196078431372549, 0.37254901960784315,
        1.51538022123186e-14, 0.43137254901960786, 0.38823529411764707, 0.34509803921568627,
        2.4253305437564208e-14, 0.403921568627451, 0.3568627450980392, 0.3176470588235294,
        3.88168471784337e-14, 0.37254901960784315, 0.3215686274509804, 0.29411764705882354,
        6.212545455928589e-14, 0.34509803921568627, 0.29411764705882354, 0.26666666666666666,
        9.943033462909e-14, 0.3176470588235294, 0.2627450980392157, 0.23921568627450981,
        1.591359212513811e-13, 0.28627450980392155, 0.23137254901960785, 0.21176470588235294,
        2.5469331393677965e-13, 0.2549019607843137, 0.2, 0.1843137254901961,
        4.076306823374471e-13, 0.23137254901960785, 0.17254901960784313, 0.16470588235294117,
        6.52403357648163e-13, 0.2, 0.1450980392156863, 0.13725490196078433,
        1.0444448411802018e-12, 0.14902, 0.196078, 0.278431,
        2.6121430589709417e-12, 0.2, 0.2549019607843137, 0.34509803921568627,
        6.5329360551199055e-12, 0.24705882352941178, 0.3176470588235294, 0.41568627450980394,
        1.6338788702139022e-11, 0.3058823529411765, 0.38823529411764707, 0.49411764705882355,
        4.0863099531478304e-11, 0.37254901960784315, 0.4588235294117647, 0.5686274509803921,
        1.0219808418851138e-10, 0.44313725490196076, 0.5333333333333333, 0.6431372549019608,
        2.555960886852512e-10, 0.5176470588235295, 0.615686274509804, 0.7254901960784313,
        6.392425168234573e-10, 0.6, 0.6980392156862745, 0.8,
        1.5987372788712038e-09, 0.6862745098039216, 0.7843137254901961, 0.8705882352941177,
        3.9984212870471895e-09, 0.7607843137254902, 0.8588235294117647, 0.9294117647058824,
        6.323307114989109e-09, 0.807843137254902, 0.9019607843137255, 0.9607843137254902,
        1e-08, 0.8901960784313725, 0.9568627450980393, 0.984313725490196,
    ],
    UseLogScale=1,
    ScalarRangeInitialized=1.0,
    EnableOpacityMapping=0,
)

# trace defaults for the display properties.
slice1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Density [gr/cm^3]'],
    LookupTable=densitygrcm3LUT,
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
slice1Display.OpacityTransferFunction.Points = [1e-20, 0.0,
                                                5e-17, 0.0, 5e-17, 1.0, 1e-8, 1.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for densitygrcm3LUT in view renderView4
densitygrcm3LUTColorBar = GetScalarBar(densitygrcm3LUT, renderView4)
densitygrcm3LUTColorBar.Set(
    Title='Density [gr/cm^3]',
    ComponentTitle='',
)

# set color bar visibility
densitygrcm3LUTColorBar.Visibility = 1

# show color legend
slice1Display.SetScalarBarVisibility(renderView4, True)

# get opacity transfer function/opacity map for 'Densitygrcm3'
densitygrcm3PWF = GetOpacityTransferFunction('Densitygrcm3')
densitygrcm3PWF.Set(
    Points=[1e-20, 0.0, 0.5, 0.0,
            5e-17, 0.0, 0.5, 0.0,
            5.1e-17, 1.0, 0.5, 0.0,
            1e-8, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

SetActiveSource(slice1)

AddCameraLink(renderView4, renderView3, 'CameraLink0')

SaveScreenshot(filename='/local/data/Walder/amaze.png', viewOrLayout=layout1_1, location=16, SaveAllViews=1)

