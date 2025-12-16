# tested Wed May  6 16:15:16 CEST 2015
# with visit2.9.1 on cygnus

# simply execute:  visit -nowin -ni -cli -s jet3d.py


DefineScalarExpression("procid", "procid(mesh)")

opts = GetDefaultFileOpenOptions("AMAZE")
## 'Apply Length Scale Factor': 1
## 'Compute Log10': 1
## 'Use Shifted Grid Info': 0
## 'BoxLength': 3
opts['Compute Log10'] = 1
opts['BoxLength'] = 3
SetDefaultFileOpenOptions("AMAZE", opts)

OpenDatabase("localhost:/local/data/Walder/jet3d.amr5", 0, "AMAZE_1.0")
AddPlot("Pseudocolor", "Density", 1, 1)
AddOperator("Isosurface", 1)
AddOperator("Clip", 1)
SetActivePlots(0)

IsosurfaceAtts = IsosurfaceAttributes()
IsosurfaceAtts.contourNLevels = 10
IsosurfaceAtts.contourValue = (4)
IsosurfaceAtts.contourPercent = ()
IsosurfaceAtts.contourMethod = IsosurfaceAtts.Level  # Level, Value, Percent
IsosurfaceAtts.scaling = IsosurfaceAtts.Linear  # Linear, Log
IsosurfaceAtts.variable = "Density"
SetOperatorOptions(IsosurfaceAtts, 1)

ClipAtts = ClipAttributes()
ClipAtts.quality = ClipAtts.Fast  # Fast, Accurate
ClipAtts.funcType = ClipAtts.Plane  # Plane, Sphere
ClipAtts.plane1Status = 1
ClipAtts.plane2Status = 0
ClipAtts.plane3Status = 0
ClipAtts.plane1Origin = (0.25, 0, 0)
ClipAtts.plane2Origin = (0, 0, 0)
ClipAtts.plane3Origin = (0, 0, 0)
ClipAtts.plane1Normal = (-1, 0, 0)
ClipAtts.plane2Normal = (0, 1, 0)
ClipAtts.plane3Normal = (0, 0, 1)
ClipAtts.planeInverse = 0
ClipAtts.planeToolControlledClipPlane = ClipAtts.Plane1  # None, Plane1, Plane2, Plane3
ClipAtts.center = (0, 0, 0)
ClipAtts.radius = 1
ClipAtts.sphereInverse = 0
SetOperatorOptions(ClipAtts, 1)

PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.colorTableName = "hot_desaturated"
SetPlotOptions(PseudocolorAtts)

annot = AnnotationAttributes()
annot.axes3D.visible = 0
annot.axes3D.triadFlag = 1
annot.axes3D.bboxFlag = 0
annot.userInfoFlag = 1
annot.databaseInfoFlag = 1
SetAnnotationAttributes(annot)

DrawPlots()

# Begin spontaneous state
View3DAtts = View3DAttributes()
View3DAtts.viewNormal = (-0.833305, 0.262716, -0.486399)
View3DAtts.focus = (0.13334, 0.291087, 0.477535)
View3DAtts.viewUp = (0.266919, 0.961712, 0.0621539)
View3DAtts.viewAngle = 30
View3DAtts.parallelScale = 0.612372
View3DAtts.nearPlane = -1.22474
View3DAtts.farPlane = 1.22474
View3DAtts.imagePan = (0.320388, 0.0294118)
View3DAtts.imageZoom = 6.77344
View3DAtts.perspective = 1
View3DAtts.eyeAngle = 2
View3DAtts.centerOfRotationSet = 1
View3DAtts.centerOfRotation = (0.2609, 0.238377, 0.0704325)
View3DAtts.axis3DScaleFlag = 0
View3DAtts.axis3DScales = (1, 1, 1)
View3DAtts.shear = (0, 0, 1)
SetView3D(View3DAtts)
# End spontaneous state

SaveWindowAtts = SaveWindowAttributes()
SaveWindowAtts.outputToCurrentDirectory = 1
SaveWindowAtts.outputDirectory = "."
SaveWindowAtts.fileName = "jet3d"
SaveWindowAtts.family = 1
SaveWindowAtts.format = SaveWindowAtts.PNG  # BMP, CURVE, JPEG, OBJ, PNG, POSTSCRIPT, POVRAY, PPM, RGB, STL, TIFF, ULTRA, VTK, PLY
SaveWindowAtts.width = 1024
SaveWindowAtts.height = 1024
SaveWindowAtts.screenCapture = 0
SetSaveWindowAttributes(SaveWindowAtts)
SaveWindow()

ChangeActivePlotsVar("procid")
SaveWindow()



