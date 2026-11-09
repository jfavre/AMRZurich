# state file generated using paraview version 6.2.0
import paraview

paraview.compatibility.major = 6
paraview.compatibility.minor = 2

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# Create a new 'Render View'
renderView2 = CreateView('RenderView')
renderView2.Set(
    #ViewSize=[802, 566],
    OrientationAxesVisibility=0,
    CenterOfRotation=[50398113038336.0, 49542537936896.0, 50000000188416.0],
    CameraPosition=[50264423932789.664, 49678295576621.61, 57649597726338.12],
    CameraFocalPoint=[50084221597285.12, 50019496616327.27, 47401344924361.61],
    CameraViewUp=[-0.0009620086755001139, 0.9994452038721633, 0.033292026015822616],
)

# Create a new 'Render View'
renderView4 = CreateView('RenderView')
renderView4.Set(
    #ViewSize=[801, 566],
    CenterOfRotation=[50398113038336.0, 49542537936896.0, 50000000188416.0],
    CameraPosition=[50264423932789.664, 49678295576621.61, 57649597726338.12],
    CameraFocalPoint=[50084221597285.12, 50019496616327.27, 47401344924361.61],
    CameraViewUp=[-0.0009620086755001139, 0.9994452038721633, 0.033292026015822616],
)

# Create a new 'Render View'
renderView6 = CreateView('RenderView')
renderView6.Set(
    #ViewSize=[632, 566],
    CenterOfRotation=[50398113038336.0, 49542537936896.0, 50000000188416.0],
    CameraPosition=[50264423932789.664, 49678295576621.61, 57649597726338.12],
    CameraFocalPoint=[50084221597285.12, 50019496616327.27, 47401344924361.61],
    CameraViewUp=[-0.0009620086755001139, 0.9994452038721633, 0.033292026015822616],
)

# Create a new 'Render View'
renderView7 = CreateView('RenderView')
renderView7.Set(
    #ViewSize=[633, 566],
    CenterOfRotation=[50398113038336.0, 49542537936896.0, 50000000188416.0],
    CameraPosition=[50264423932789.664, 49678295576621.61, 57649597726338.12],
    CameraFocalPoint=[50084221597285.12, 50019496616327.27, 47401344924361.61],
    CameraViewUp=[-0.0009620086755001139, 0.9994452038721633, 0.033292026015822616],
)

SetActiveView(None)

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.SplitVertical(0, 0.500000)
layout1.SplitHorizontal(1, 0.500000)
layout1.AssignView(3, renderView2)
layout1.AssignView(4, renderView4)
layout1.SplitHorizontal(2, 0.500000)
layout1.AssignView(5, renderView6)
layout1.AssignView(6, renderView7)
layout1.SetSize(1773, 980)
# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView2)
# ----------------------------------------------------------------

directory_JeanLaptop = '/local/data/Walder/CygX-1_M1=14.8_M2=19.2_V=750_ML=1-6_Gamma=1.01_RBH=20_HR/DATA/'
directory_r740gpu02 = '/Xnfs/highenergy/rwalder/Simulations/CygX-1_M1=14.8_M2=19.2_V=750_ML=1-6_Gamma=1.01_RBH=20_HR/'
directory_RolfLaptop ='/home/rolf/Paraview/Paraview/DATA.ToTest.Paraview/'

fname = 'DATA.CX-1_M1=14.8_M2=19.2_V=750_ML=1-6_G=1.01_RBH=20RG_.NT=0000008614.Time=189.4812456499H.amr5'

# create a new directory name for Rolf and Doris's laptops, or use Jean's laptop, or use r740gpu02
# watch for the last "/" in the directory name

#filename = directory_r740gpu02 + fname
filename = directory_JeanLaptop + fname
#filename = directory_RolfLaptop + fname

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
coordinatesFrames = ProgrammableFilter(registrationName='CoordinatesFrames', Input=axisAlignedSlice1)
coordinatesFrames.Set(
    RequestInformationScript='',
    RequestUpdateExtentScript='',
    CopyArrays=1,
    PythonPath='',
    Script="""import numpy as np
from vtkmodules.vtkCommonDataModel import vtkOverlappingAMR
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk

VECTOR_ARRAY = "Velocity [cm/s]"       # input POINTS vector array (3 components)
OUTPUT_PREFIX = ""               # e.g. "Frame_" or "" for short names
# Physical parameters from CoordFramesDemo.py:
m_sun = 1.989e33

# ----------------------------------------------------------------------
# Input / output of the nested AMR structures
# ----------------------------------------------------------------------

input_amr = self.GetInputDataObject(0, 0)
output_amr = self.GetOutputDataObject(0)

if not isinstance(input_amr, vtkOverlappingAMR):
    raise RuntimeError(
        f"Expected vtkOverlappingAMR, got {input_amr.GetClassName()}"
    )


# get coordinates for CM and the stars
CM_x = input_amr.GetFieldData().GetArray("CenterOfMass").GetValue(0)
CM_y = input_amr.GetFieldData().GetArray("CenterOfMass").GetValue(1)
CM_z = input_amr.GetFieldData().GetArray("CenterOfMass").GetValue(2)

Pos_CM = np.array([CM_x, CM_y, CM_z])

# star data is available as 7 numbers (Mass, Positions[3], and Velocity[3]
# use indexing carefully

m_BH = input_amr.GetFieldData().GetArray("SimpleAccretor").GetValue(0)  * m_sun
m_DO = input_amr.GetFieldData().GetArray("IsotrInfWind").GetValue(0)    * m_sun

p_DO = np.array([input_amr.GetFieldData().GetArray("IsotrInfWind").GetValue(1),
                      input_amr.GetFieldData().GetArray("IsotrInfWind").GetValue(2),
                      input_amr.GetFieldData().GetArray("IsotrInfWind").GetValue(3)])
                      
p_BH = np.array([input_amr.GetFieldData().GetArray("SimpleAccretor").GetValue(1),
                      input_amr.GetFieldData().GetArray("SimpleAccretor").GetValue(2),
                      input_amr.GetFieldData().GetArray("SimpleAccretor").GetValue(3)])
                      
d_BD = np.linalg.norm(p_DO - p_BH)  # BH--DO separation, in coordinate units

output_amr.CopyStructure(input_amr)
output_amr.ShallowCopy(input_amr)

# angular velocity Omega of the rotating system - get it from DATA*amr5 files
Omega = input_amr.GetFieldData().GetArray("Omega").GetValue(0)

omega_vec = np.array([0.0, 0.0, Omega], dtype=np.float64)

# compute position angle of BH in order to establish a spherical coordinate frame
# independent of it. phi_rot is the 0-direction of the spherical system
# pointing always in direction of the BH motion

# alpha: angle between x-axis of simulation frame (SF) and 
# line (vector) connecting from donor to BH (line goes also through CM)
# line does not go throuh (0/0/0)!!!
# 0.5 pi accounts for the fact that BH moves 90 degrees to that line
alpha = np.arctan2(p_BH[1] - p_DO[1], p_BH[0] - p_DO[0])
if (alpha<0.): alpha= alpha+2. *np.pi
phi_rot = alpha + 0.5 * np.pi
cos_rot, sin_rot = np.cos(phi_rot), np.sin(phi_rot)

def add_vector(point_data, name, values):
    arr = numpy_to_vtk(np.ascontiguousarray(values), deep=1)
    arr.SetName(OUTPUT_PREFIX + name)
    arr.SetNumberOfComponents(3)
    point_data.AddArray(arr)


for level in range(input_amr.GetNumberOfLevels()):
    for index in range(input_amr.GetNumberOfBlocks(level)):
        src = input_amr.GetDataSetAsImageData(level, index)
        if src is None:
            continue

        # Empty AMR slots may have no data set.
        if src.GetNumberOfPoints() == 0:
            continue

        velocity_vtk = src.GetPointData().GetArray(VECTOR_ARRAY)
        if velocity_vtk is None:
            raise RuntimeError(
                "Block at level %d, index %d has no POINTS array named \'%s\'"
                % (level, index, VECTOR_ARRAY)
            )
        if velocity_vtk.GetNumberOfComponents() != 3:
            raise RuntimeError("Input velocity array must have 3 components")

        # Make a copy of the block
        output_grid = src.NewInstance()
        output_grid.ShallowCopy(src)

        vel_SF = np.asarray(vtk_to_numpy(velocity_vtk), dtype=np.float64)
        npoints = src.GetNumberOfPoints()

        # Image point coordinates, including non-zero extents.
        extent = src.GetExtent()
        origin = np.asarray(src.GetOrigin(), dtype=np.float64)
        spacing = np.asarray(src.GetSpacing(), dtype=np.float64)
        ix = np.arange(extent[0], extent[1] + 1, dtype=np.float64)
        iy = np.arange(extent[2], extent[3] + 1, dtype=np.float64)
        iz = np.arange(extent[4], extent[5] + 1, dtype=np.float64)
        # Fortran flattening makes x vary fastest, matching VTK point ordering.
        X, Y, Z = np.meshgrid(ix, iy, iz, indexing="ij")
        coords = np.column_stack((
            (origin[0] + X.ravel(order="F") * spacing[0]),
            (origin[1] + Y.ravel(order="F") * spacing[1]),
            (origin[2] + Z.ravel(order="F") * spacing[2]),
        ))
        if coords.shape[0] != npoints:
            raise RuntimeError("Internal coordinate/point count mismatch")

        CM2GP = coords - Pos_CM
        BH2GP = coords - p_BH
        DO2GP = coords - p_DO

        # Cartesian velocities in frames translating with BH and CM.
        vel_BH = vel_SF - np.cross(
            np.broadcast_to(omega_vec, DO2GP.shape), CM2GP - BH2GP
        )

        # Spherical coordinates centered at BH. Preserve the convention from
        # CoordFramesDemo.py: phi=0 is the BH orbital-motion direction and
        # e_theta = -e_r x e_phi.
        r = np.linalg.norm(BH2GP, axis=1)
        rho = np.hypot(BH2GP[:, 0], BH2GP[:, 1])
        theta = np.zeros(npoints, dtype=np.float64)
        valid_r = r > 0.0
        theta[valid_r] = np.arccos(
            np.clip(BH2GP[valid_r, 2] / r[valid_r], -1.0, 1.0)
        )
        phi_sf = np.arccos( BH2GP[:, 0]/ rho)
        for ip in range(len(phi_sf)):
            if (BH2GP[ip, 1] < 0.): phi_sf[ip]= 2.*np.pi- phi_sf[ip]

        phi = phi_sf - alpha - 0.5 * np.pi   # phi = phi_sf - phi_rot
        phi = np.mod(phi, 2.0 * np.pi)

        vx, vy, vz = vel_BH[:, 0], vel_BH[:, 1], vel_BH[:, 2]
        vx_rot = vx * cos_rot + vy * sin_rot
        vy_rot = -vx * sin_rot + vy * cos_rot
        vz_rot = vz

        sin_theta, cos_theta = np.sin(theta), np.cos(theta)
        sin_phi, cos_phi = np.sin(phi), np.cos(phi)
        vr = sin_theta * cos_phi * vx_rot + sin_theta * sin_phi * vy_rot + cos_theta * vz_rot
        vphi = -sin_phi * vx_rot + cos_phi * vy_rot
        vtheta = cos_theta * cos_phi * vx_rot + cos_theta * sin_phi * vy_rot - sin_theta * vz_rot

        vel_BH_spherical = np.column_stack((vr, vphi, vtheta))

        # Corotating-system velocity fields, following the signs in the demo.
        vel_SYS_COROT = np.column_stack((
             Omega * CM2GP[:, 1],
            -Omega * CM2GP[:, 0],
            np.zeros(npoints),
         ))
        
        # this is the co-rotating velocity in the simulation frame
        #
        vel_CoRot = np.column_stack((
            vel_SF[:, 0] + Omega * CM2GP[:, 1],
            vel_SF[:, 1] - Omega * CM2GP[:, 0],
            vel_SF[:, 2],
        ))

        pd = output_grid.GetPointData()
        add_vector(pd, "Velocity_BH_frame", vel_BH)
        add_vector(pd, "Velocity_BH_spherical", vel_BH_spherical)
        # add_vector(pd, "Velocity_corotating_frame", vel_CoRot)

        # Put block into AMR output
        output_amr.SetDataSet(level, index, output_grid)

# Output remains vtkOverlappingAMR with the same hierarchy and per-block geometry.""",
)


# create a new 'Extract Block'
extractBlock1 = ExtractBlock(registrationName='ExtractBlock1', Input=coordinatesFrames)
extractBlock1.Set(
    Assembly='Hierarchy',
    Selectors=['/Root/Level5', '/Root/Level6', '/Root/Level7', '/Root/Level8', '/Root/Level9', '/Root/Level10', '/Root/Level11'],
)

# create a new 'Resample To Image'
resampleToImage1 = ResampleToImage(registrationName='ResampleToImage1', Input=extractBlock1)

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
velocity = Glyph(registrationName='Velocity', Input=maskPoints1,
    GlyphType='Arrow')
velocity.Set(
    OrientationArray=['POINTS', 'Velocity [cm/s]'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=2.5e11,
    GlyphMode='All Points',
)

# create a new 'Glyph'
velocity_BH_spherical = Glyph(registrationName='Velocity_BH_spherical', Input=maskPoints1,
    GlyphType='Arrow')
velocity_BH_spherical.Set(
    OrientationArray=['POINTS', 'Velocity_BH_spherical'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=2.5e11,
    GlyphMode='All Points',
)


# create a new 'Glyph'
velocity_BH_frame = Glyph(registrationName='Velocity_BH_frame', Input=maskPoints1,
    GlyphType='Arrow')
velocity_BH_frame.Set(
    OrientationArray=['POINTS', 'Velocity_BH_frame'],
    ScaleArray=['POINTS', 'No scale array'],
    ScaleFactor=2.5e11,
    GlyphMode='All Points',
)


# ----------------------------------------------------------------
# setup the visualization in view 'renderView2'
# ----------------------------------------------------------------

# show data from reader_1
reader_1Display = Show(OutputPort(reader, 1), renderView2, 'GeometryRepresentation')

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

# init the 'Grid Axes Representation' selected for 'DataAxesGrid'
reader_1Display.DataAxesGrid.Set(
    XTitleBold=0,
    XTitleFontSize=12,
    YTitleBold=0,
    YTitleFontSize=12,
    ZTitleBold=0,
    ZTitleFontSize=12,
    XLabelBold=0,
    XLabelFontSize=12,
    YLabelBold=0,
    YLabelFontSize=12,
    ZLabelBold=0,
    ZLabelFontSize=12,
)

# init the 'Polar Axes Representation' selected for 'PolarAxes'
reader_1Display.PolarAxes.Set(
    PolarAxisTitleBold=0,
    PolarAxisTitleFontSize=12,
    PolarAxisLabelBold=0,
    PolarAxisLabelFontSize=12,
    LastRadialAxisTextBold=0,
    LastRadialAxisTextFontSize=12,
    SecondaryRadialAxesTextBold=0,
    SecondaryRadialAxesTextFontSize=12,
)

# show data from extractBlock1
extractBlock1Display = Show(extractBlock1, renderView2, 'UniformGridRepresentation')

densitygrcm3LUT = GetColorTransferFunction('Densitygrcm3')
densitygrcm3LUT.Set(
    RGBPoints=[
        # scalar, red, green, blue
        3.083465133715197e-15, 0.0564, 0.0564, 0.47,
        5.273912119552955e-14, 0.243, 0.46035, 0.81,
        4.305801958306158e-13, 0.356814, 0.745025, 0.954368,
        3.930094489031055e-12, 0.6882, 0.93, 0.91791,
        1.2081958547194573e-11, 0.899496, 0.944646, 0.768657,
        5.201717021545207e-11, 0.957108, 0.833819, 0.508916,
        3.6602191964106567e-10, 0.927521, 0.621439, 0.315357,
        3.8049304864845845e-09, 0.8, 0.352, 0.16,
        4.7340805232406797e-08, 0.59, 0.0767, 0.119475,
    ],
    UseLogScale=1,
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
extractBlock1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Density [gr/cm^3]'],
    LookupTable=densitygrcm3LUT,
    Assembly='Assembly',
)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView4'
# ----------------------------------------------------------------

# show data from reader_1
reader_1Display_2 = Show(OutputPort(reader, 1), renderView4, 'GeometryRepresentation')

# trace defaults for the display properties.
reader_1Display_2.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
reader_1Display_2.ScaleTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 19.2, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
reader_1Display_2.OpacityTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 19.2, 1.0, 0.5, 0.0]

# show data from extractBlock1
extractBlock1Display_2 = Show(extractBlock1, renderView4, 'UniformGridRepresentation')

# get color transfer function/color map for 'Velocitycms'
velocitycmsLUT = GetColorTransferFunction('Velocitycms')
velocitycmsLUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=GenerateRGBPoints(
        range_min=4036860.458995151,
        range_max=200000000.0,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
extractBlock1Display_2.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity [cm/s]'],
    LookupTable=velocitycmsLUT,
    Assembly='Assembly',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
extractBlock1Display_2.ScaleTransferFunction.Points = [3.0834651337151977e-15, 0.0, 0.5, 0.0, 7.896898444052566e-09, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
extractBlock1Display_2.OpacityTransferFunction.Points = [3.0834651337151977e-15, 0.0, 0.5, 0.0, 7.896898444052566e-09, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
extractBlock1Display_2.SliceFunction.Origin = [49781250000000.0, 50000000000000.0, 50000000000000.0]

# show data from velocity
velocityDisplay = Show(velocity, renderView4, 'GeometryRepresentation')

# trace defaults for the display properties.
velocityDisplay.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
velocityDisplay.ScaleTransferFunction.Points = [3.567656136532961e-15, 0.0, 0.5, 0.0, 8.083883587512268e-12, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
velocityDisplay.OpacityTransferFunction.Points = [3.567656136532961e-15, 0.0, 0.5, 0.0, 8.083883587512268e-12, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for velocitycmsLUT in view renderView4
velocitycmsLUTColorBar = GetScalarBar(velocitycmsLUT, renderView4)
velocitycmsLUTColorBar.Set(
    Title='Velocity [cm/s]',
    ComponentTitle='Magnitude',
)

# set color bar visibility
velocitycmsLUTColorBar.Visibility = 1

# show color legend
extractBlock1Display_2.SetScalarBarVisibility(renderView4, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView6'
# ----------------------------------------------------------------

# show data from reader_1
reader_1Display_4 = Show(OutputPort(reader, 1), renderView6, 'GeometryRepresentation')

# trace defaults for the display properties.
reader_1Display_4.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
reader_1Display_4.ScaleTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 19.2, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
reader_1Display_4.OpacityTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 19.2, 1.0, 0.5, 0.0]

# show data from extractBlock1
extractBlock1Display_4 = Show(extractBlock1, renderView6, 'UniformGridRepresentation')

# get color transfer function/color map for 'Velocity_BH_frame'
velocity_BH_frameLUT = GetColorTransferFunction('Velocity_BH_frame')
velocity_BH_frameLUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=GenerateRGBPoints(
        range_min=2558701.2061252985,
        range_max=100000000.0,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
extractBlock1Display_4.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity_BH_frame'],
    LookupTable=velocity_BH_frameLUT,
    Assembly='Assembly',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
extractBlock1Display_4.ScaleTransferFunction.Points = [3.0834651337151977e-15, 0.0, 0.5, 0.0, 4.734080523240675e-08, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
extractBlock1Display_4.OpacityTransferFunction.Points = [3.0834651337151977e-15, 0.0, 0.5, 0.0, 4.734080523240675e-08, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
extractBlock1Display_4.SliceFunction.Origin = [49781250000000.0, 50000000000000.0, 50000000000000.0]

# show data from velocity_BH_frame
velocity_BH_frameDisplay = Show(velocity_BH_frame, renderView6, 'GeometryRepresentation')

# trace defaults for the display properties.
velocity_BH_frameDisplay.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
velocity_BH_frameDisplay.ScaleTransferFunction.Points = [3.567656136532961e-15, 0.0, 0.5, 0.0, 8.083883587512268e-12, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
velocity_BH_frameDisplay.OpacityTransferFunction.Points = [3.567656136532961e-15, 0.0, 0.5, 0.0, 8.083883587512268e-12, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for velocity_BH_frameLUT in view renderView6
velocity_BH_frameLUTColorBar = GetScalarBar(velocity_BH_frameLUT, renderView6)
velocity_BH_frameLUTColorBar.Set(
    Title='Velocity_BH_frame',
    ComponentTitle='Magnitude',
)

# set color bar visibility
velocity_BH_frameLUTColorBar.Visibility = 1

# show color legend
extractBlock1Display_4.SetScalarBarVisibility(renderView6, True)

# ----------------------------------------------------------------
# setup the visualization in view 'renderView7'
# ----------------------------------------------------------------

# show data from reader_1
reader_1Display_5 = Show(OutputPort(reader, 1), renderView7, 'GeometryRepresentation')

# trace defaults for the display properties.
reader_1Display_5.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
    SelectNormalArray='Normals',
    Assembly='Hierarchy',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
reader_1Display_5.ScaleTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 19.2, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
reader_1Display_5.OpacityTransferFunction.Points = [14.8, 0.0, 0.5, 0.0, 19.2, 1.0, 0.5, 0.0]

# show data from extractBlock1
extractBlock1Display_5 = Show(extractBlock1, renderView7, 'UniformGridRepresentation')

# get color transfer function/color map for 'Velocity_BH_spherical'
velocity_BH_sphericalLUT = GetColorTransferFunction('Velocity_BH_spherical')
velocity_BH_sphericalLUT.Set(
    AutomaticRescaleRangeMode='Never',
    RGBPoints=GenerateRGBPoints(
        range_min=4036860.458995151,
        range_max=100000000.0,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
extractBlock1Display_5.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', 'Velocity_BH_spherical'],
    LookupTable=velocity_BH_sphericalLUT,
    Assembly='Assembly',
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
extractBlock1Display_5.ScaleTransferFunction.Points = [3.0834651337151977e-15, 0.0, 0.5, 0.0, 4.734080523240675e-08, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
extractBlock1Display_5.OpacityTransferFunction.Points = [3.0834651337151977e-15, 0.0, 0.5, 0.0, 4.734080523240675e-08, 1.0, 0.5, 0.0]

# init the 'Plane' selected for 'SliceFunction'
extractBlock1Display_5.SliceFunction.Origin = [49781250000000.0, 50000000000000.0, 50000000000000.0]

# show data from velocity_BH_spherical
velocity_BH_sphericalDisplay = Show(velocity_BH_spherical, renderView7, 'GeometryRepresentation')

# trace defaults for the display properties.
velocity_BH_sphericalDisplay.Set(
    Representation='Surface',
    ColorArrayName=[None, ''],
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
velocity_BH_sphericalDisplay.ScaleTransferFunction.Points = [3.567656136532961e-15, 0.0, 0.5, 0.0, 8.083883587512268e-12, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
velocity_BH_sphericalDisplay.OpacityTransferFunction.Points = [3.567656136532961e-15, 0.0, 0.5, 0.0, 8.083883587512268e-12, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for velocity_BH_sphericalLUT in view renderView7
velocity_BH_sphericalLUTColorBar = GetScalarBar(velocity_BH_sphericalLUT, renderView7)
velocity_BH_sphericalLUTColorBar.Set(
    Title='Velocity_BH_spherical',
    ComponentTitle='Magnitude',
)

# set color bar visibility
velocity_BH_sphericalLUTColorBar.Visibility = 1

# show color legend
extractBlock1Display_5.SetScalarBarVisibility(renderView7, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'Velocity_corotation_only'
velocity_corotation_onlyPWF = GetOpacityTransferFunction('Velocity_corotation_only')
velocity_corotation_onlyPWF.Set(
    Points=[0.0, 0.0, 0.5, 0.0, 40000000.0, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'Velocity_BH_frame'
velocity_BH_framePWF = GetOpacityTransferFunction('Velocity_BH_frame')
velocity_BH_framePWF.Set(
    Points=[2558701.2061252985, 0.0, 0.5, 0.0, 100000000.0, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'Velocitycms'
velocitycmsPWF = GetOpacityTransferFunction('Velocitycms')
velocitycmsPWF.Set(
    Points=[4036860.458995151, 0.0, 0.5, 0.0, 200000000.0, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'Velocity_CM_frame'
velocity_CM_framePWF = GetOpacityTransferFunction('Velocity_CM_frame')
velocity_CM_framePWF.Set(
    Points=[1000.0, 0.0, 0.5, 0.0, 100000000.0, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'Velocity_corotating_frame'
velocity_corotating_framePWF = GetOpacityTransferFunction('Velocity_corotating_frame')
velocity_corotating_framePWF.Set(
    Points=[7608457.51153899, 0.0, 0.5, 0.0, 200000000.0, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# get opacity transfer function/opacity map for 'Velocity_BH_spherical'
velocity_BH_sphericalPWF = GetOpacityTransferFunction('Velocity_BH_spherical')
velocity_BH_sphericalPWF.Set(
    Points=[4036860.458995151, 0.0, 0.5, 0.0, 100000000.0, 1.0, 0.5, 0.0],
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
timeKeeper1.SuppressedTimeSources = reader

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.Set(
    ViewModules=[renderView2, renderView4, renderView6, renderView7],
    Cues=timeAnimationCue1,
    AnimationTime=682111.484339536,
    StartTime=682111.484339536,
    EndTime=682111.484339536,
    PlayMode='Snap To TimeSteps',
)

# ----------------------------------------------------------------
# restore active source
SetActiveSource(extractBlock1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup camera links
AddCameraLink(renderView2, renderView4, 'CameraLink0')

AddCameraLink(renderView2, renderView6, 'CameraLink3')

AddCameraLink(renderView7, renderView2, 'CameraLink4')
# ----------------------------------------------------------------


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/nightly/python/
##--------------------------------------------
