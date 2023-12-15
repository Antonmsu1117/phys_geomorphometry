# Description: the script is aimed to compute indices for physical geomorphometry analysis.
# It applies DEM to compute hierarchical land surface variables up to third order and shapefile
# of segmentation to compute the general indices of physical geomorphometry.
# Author: Anton Popov
# The script was tested in IDLE environment for ArcGIS pro, ver. 3.9.18

import arcpy
import numpy as np
import math
from arcpy.sa import *

arcpy.env.overwriteOutput = True
arcpy.env.nodata = "NONE"
arcpy.env.workspace=r"D:\Dissertation\Toolbox_Physical_Geomorphometry\outputs\output_variables"

DEM = r"D:\Dissertation\Toolbox_Physical_Geomorphometry\inputs\DEM.tif"
segment_layer = r"D:\Dissertation\Toolbox_Physical_Geomorphometry\inputs\Segmentation_18_2slope_2ktn_ss20.shp"

DEM_array = arcpy.RasterToNumPyArray(DEM)


# Get the coordinate system of the DEM
desc = arcpy.Describe(DEM)
coordinate_system = desc.spatialReference

def set_coordinate_system(output_raster):
    arcpy.DefineProjection_management(output_raster, coordinate_system)


# Get raster metrics for coordinate info
imgRaster = arcpy.sa.Raster(DEM)


# Lower left coordinate of block (in map units)
mx = imgRaster.extent.XMin
my = imgRaster.extent.YMin


#Computation of partial derivatives
zy = np.gradient(DEM_array, axis=0)
zx = np.gradient(DEM_array, axis=1)
zxx = np.gradient(zx, axis=1)
zyy = np.gradient(zy, axis=0)
zxy = np.gradient(zx, axis=0)
zxxx = np.gradient(zxx, axis=1)
zxxy = np.gradient(zxy, axis=1)
zxxy = np.gradient(zxxy, axis=0)
zxyy = np.gradient(zx, axis=1)
zxyy = np.gradient(zxyy, axis=0)
zyyy = np.gradient(zyy, axis=0)


#Computation of coefficients for simplification of the formula
A = zxy**2-(zxx*zyy)
B = 2*zx*(A+zxxy*zy)-zy**2*zxxx-zx**2*zxyy
C = 2*zy*(A+zxyy*zx)-zy**2.*zxxy-zx**2*zyyy
D = zxy*(zxx+zyy)
E = -2*( zx*(zxxy*zy+zxy**2+zxx**2)+D*zy ) - zxxx*zx**2-zxyy*zy**2
M = zy**2+zx**2
L = -2*zxy*zx*zy-zyy*zy**2-zxx*zx**2
F = -2*( zy*(zxyy*zx+zxy**2+zyy**2)+D*zx ) - zyyy*zy**2-zxxy*zx**2
G = zy**2-zx**2
N = zxx*zx+zxy*zy
H = zxx-zyy
K = 2*zxy*zx*zy-zxx*zy**2-zyy*zx**2
O = zyy*zy+zxy*zx
I = zxxy*G+2*zxy*(zxy*zy-zxx*zx)+H*(zxx*zy+zxy*zx)+zx*zy*(zxxx-zxyy)
P = 1+M
R = zx*zy*H+zxy*G
J = zxyy*G+2*zxy*(zyy*zy-zxy*zx)+H*(zxy*zy+zyy*zx)+zx*zy*(zxxy-zyyy)
U = zxx*(1+zy**2)
V = zyy*(1+zx**2)
W = 2*zxy*zx*zy
Y = 1 + zx**2 + zy**2


##Computation of land surface variables


####SLOPE
slope = np.arctan(np.sqrt(zx**2 + zy**2)) *(180 / np.pi)
sin_slope = np.sin(np.radians(slope))
cos_slope = np.cos(np.radians(slope))
ctg_slope = cos_slope / sin_slope
##slope_deg=np.degrees(slope)

arcpy.NumPyArrayToRaster(slope, arcpy.Point(mx, my), 1, 1).save("slope_deg.tif")
set_coordinate_system("slope_deg.tif")

arcpy.NumPyArrayToRaster(sin_slope, arcpy.Point(mx, my), 1, 1).save("sin_slope.tif")
set_coordinate_system("sin_slope.tif")

arcpy.NumPyArrayToRaster(cos_slope, arcpy.Point(mx, my), 1, 1).save("cos_slope.tif")
set_coordinate_system("cos_slope.tif")

arcpy.NumPyArrayToRaster(ctg_slope, arcpy.Point(mx, my), 1, 1).save("ctg_slope.tif")
set_coordinate_system("ctg_slope.tif")
##arcpy.RasterToASCII_conversion(slope_raster,
##                               r"D:\Dissertation\Toolbox_Physical_Geomorphometry\Variables\Slope_deg.asc")


####ASPECT
Aspect = -90*(1-np.sign(zy))*(1-np.abs(np.sign(zx))) + 180*(1+np.sign(zx)) - 180*((np.sign(zx))*np.arccos((-zy)/(np.sqrt(zx**2+zy**2))))/np.pi
arcpy.NumPyArrayToRaster(Aspect, arcpy.Point(mx, my), 1, 1).save("asp_rad.tif")
set_coordinate_system("asp_rad.tif")

sin_asp = np.sin(Aspect*np.pi/180)
cos_asp = np.cos(Aspect*np.pi/180)
arcpy.NumPyArrayToRaster(sin_asp, arcpy.Point(mx, my), 1, 1).save("sin_asp.tif")
set_coordinate_system("sin_asp.tif")
arcpy.NumPyArrayToRaster(cos_asp, arcpy.Point(mx, my), 1, 1).save("cos_asp.tif")
set_coordinate_system("cos_asp.tif")


##Kn
Kn = -1*(zxx*zx**2+2*zxy*zx*zy+zyy*zy**2)/(
    (zx**2+zy**2)*np.sqrt((1+zx**2+zy**2)**3))

arcpy.NumPyArrayToRaster(Kn, arcpy.Point(mx,my), 1, 1).save("Kn.tif")
set_coordinate_system("Kn.tif")


##At
At = -(zxx*zy**2-2*zxy*zx*zy+zyy*zx**2)/np.sqrt((zx**2+zy**2)**3)
At_raster=arcpy.NumPyArrayToRaster(At, arcpy.Point(mx,my), 1, 1).save("At.tif")
set_coordinate_system("At.tif")


##Kt
Kt = -(zxx*zy**2-2*zxy*zx*zy+zyy*zx**2)/((zx**2+zy**2)*np.sqrt(1+zx**2+zy**2))
Kt_raster=arcpy.NumPyArrayToRaster(Kt, arcpy.Point(mx,my), 1, 1).save("Kt.tif")
set_coordinate_system("Kt.tif")


##St
St = (zx*zy*(zxx-zyy)-zxy*(zx**2-zy**2))/((zx**2+zy**2)*(1+zx**2+zy**2))
St_raster=arcpy.NumPyArrayToRaster(St, arcpy.Point(mx,my), 1, 1).save("St.tif")
set_coordinate_system("St.tif")

##Kd
Kd= (Kn-Kt)/2
Kd_raster=arcpy.NumPyArrayToRaster(Kd, arcpy.Point(mx,my), 1, 1).save("Kd.tif")
set_coordinate_system("Kd.tif")

##Kmean
Kmean = -(((1 + zy**2) * zxx - 2 * zxy * zx * zy + (1 + zx**2) * zyy) / (2 * np.sqrt((1 + zx**2 + zy**2)**3)))
Kmean_raster=arcpy.NumPyArrayToRaster(Kmean, arcpy.Point(mx,my), 1, 1).save("Kmean.tif")
set_coordinate_system("Kmean.tif")

##Ku
Ku = np.sqrt((((((1 + zy**2) * zxx - 2 * zxy * zx * zy + (1 + zx**2) * zyy) / (2 * np.sqrt((1 + zx**2 + zy**2)**3)))**2) - ((zxx * zyy - zxy**2) / (Y**2))))
Ku_raster=arcpy.NumPyArrayToRaster(Ku, arcpy.Point(mx,my), 1, 1).save("Ku.tif")
set_coordinate_system("Ku.tif")

##Kmax
Kmax = Kmean + Ku
Kmax_raster=arcpy.NumPyArrayToRaster(Kmax, arcpy.Point(mx,my), 1, 1).save("Kmax.tif")
set_coordinate_system("Kmax.tif")

##Kmin
Kmin = Kmean - Ku
Kmin_raster=arcpy.NumPyArrayToRaster(Kmin, arcpy.Point(mx,my), 1, 1).save("Kmin.tif")
set_coordinate_system("Kmin.tif")

##Kc
Kc = np.sqrt((Kmin**2 + Kmax**2) / 2)
Kc_raster=arcpy.NumPyArrayToRaster(Kc, arcpy.Point(mx,my), 1, 1).save("Kc.tif")
set_coordinate_system("Kc.tif")

##Knn
Knn = ( P*(zx*(2*L*N-E*M)-zy*(F*M-2*L*O))+3*L*M*(O*zy+N*zx))/((P*M)**2.5)
Knn_raster=arcpy.NumPyArrayToRaster(Knn, arcpy.Point(mx,my), 1, 1).save("Knn.tif")
set_coordinate_system("Knn.tif")

####Knt
##Knt = ( P*(zy*(E*M-L*N)-zx*(F*M-2*L*O))+3*L*M*(O*zx-N*zy))/((P*M)**2.5)
##Knt_raster=arcpy.NumPyArrayToRaster(Knt, arcpy.Point(mx,my), 1, 1).save("Knt.tif")
##set_coordinate_system("Knt.tif")

##Ktn
Ktn = (zx*(K*(2*N*P+M*N)-B*M*P)-zy*(C*M*P-K*(2*O*P+M*O)))/((P**1.5)*(M**2.5))
Ktn_raster=arcpy.NumPyArrayToRaster(Ktn, arcpy.Point(mx,my), 1, 1).save("Ktn.tif")
set_coordinate_system("Ktn.tif")

##Ktt
Ktt = (zy*(B*M*P- K*(2*N*P+M*N))-zx*(C*M*P-K*(2*O*P+M*O)))/((P**1.5)*(M**2.5))
Ktt_raster=arcpy.NumPyArrayToRaster(Ktt, arcpy.Point(mx,my), 1, 1).save("Ktt.tif")
set_coordinate_system("Ktt.tif")

####Stn
##Stn = (zx*(2*N*R*(1+2*M)-I*M*P)+zy*(2*O*R*(1+2*M)-J*M*P))/((P**2)*(M**2.5))
##Stn_raster=arcpy.NumPyArrayToRaster(Stn, arcpy.Point(mx,my), 1, 1).save("Stn.tif")
##set_coordinate_system("Stn.tif")


##Physical_geomorphometry indices

##ISED
ISED = abs(50 * (Kd / sin_slope))
ISED_raster=arcpy.NumPyArrayToRaster(ISED, arcpy.Point(mx,my), 1, 1).save("ISED.tif")
set_coordinate_system("ISED.tif")

##IGDED
IGDED = 100 * (abs(St) / ctg_slope)
IGDED_raster=arcpy.NumPyArrayToRaster(IGDED, arcpy.Point(mx,my), 1, 1).save("IGDED.tif")
set_coordinate_system("IGDED.tif")

##PESe
PESe = abs(2 * Kd)
PESe_raster=arcpy.NumPyArrayToRaster(PESe, arcpy.Point(mx,my), 1, 1).save("PESe.tif")
set_coordinate_system("PESe.tif")

##PESst
PESst = abs(St) * ctg_slope
PESst_raster=arcpy.NumPyArrayToRaster(PESst, arcpy.Point(mx,my), 1, 1).save("PESst.tif")
set_coordinate_system("PESst.tif")

##PESD
PESD = abs(Kmean * 2)
PESD_raster=arcpy.NumPyArrayToRaster(PESD, arcpy.Point(mx,my), 1, 1).save("PESD.tif")
set_coordinate_system("PESD.tif")

##IPESC
IPESC = 2 * Kc
IPESC_raster=arcpy.NumPyArrayToRaster(IPESC, arcpy.Point(mx,my), 1, 1).save("IPESC.tif")
set_coordinate_system("IPESC.tif")


print("Computation of variables successful")
#-------------------------------------------------------------------------------------

##Get zonal statistic for segmented territory 

sin_slope = "sin_slope.tif"
sin_asp = "sin_asp.tif"
cos_asp = "cos_asp.tif"
Kn = "Kn.tif"
Knn = "Knn.tif"
Kt = "Kt.tif"
Ktn = "Ktn.tif"
Ktt = "Ktt.tif"
St = "St.tif"
ISED = "ISED.tif"
IGDED = "IGDED.tif"
PESe = "PESe.tif"
PESst = "PESst.tif"
PESD = "PESD.tif"
IPESC = "IPESC.tif"

#################DEM##################
# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", DEM,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", DEM,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", DEM,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", DEM,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_DEM', 'Std_DEM', 'Range_DEM', 'Median_DEM']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_DEM'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_DEM'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_DEM'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Median_DEM'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_DEM', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "(3.46 * !STD_DEM!) / (math.sqrt(!AREA_DEM!) * math.tan(0.2))"
fieldname = "CV_DEM"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)

#################SIN_SLOPE##################
# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", sin_slope,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", sin_slope,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", sin_slope,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", sin_slope,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_s_sl', 'Std_s_sl', 'Range_s_sl', 'Median_s_s']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_s_sl'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_s_sl'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_s_sl'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Median_s_s'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_sin_sl', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_s_sl!/abs(!Mean_s_sl!)"
fieldname = "CV_sin_sl"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)

#################SIN_ASP##################
# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", sin_asp,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", sin_asp,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", sin_asp,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", sin_asp,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_s_a', 'Std_s_a', 'Range_s_a', 'Median_s_a']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_s_a'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_s_a'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_s_a'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Median_s_a'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_sin_asp', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_s_a!/abs(!Mean_s_sl!)"
fieldname = "CV_sin_asp"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)

#################COS_ASP##################
# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", cos_asp,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", cos_asp,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", cos_asp,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", cos_asp,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_c_a', 'Std_c_a', 'Range_c_a', 'Median_c_a']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_c_a'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_c_a'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_c_a'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Median_c_a'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_cos_asp', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_c_a!/abs(!Mean_c_a!)"
fieldname = "CV_cos_asp"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)


#################KN##################
# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Kn,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Kn,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Kn,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Kn,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_Kn', 'Std_Kn', 'Range_Kn', 'Median_Kn']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_Kn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_Kn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_Kn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Median_Kn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_Kn', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_Kn!/abs(!Mean_Kn!)"
fieldname = "CV_Kn"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)

#################KNN##################
# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Knn,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Knn,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Knn,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Knn,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_Knn', 'Std_Knn', 'Range_Knn', 'Median_Knn']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_Knn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_Knn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_Knn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Median_Knn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_Knn', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_Knn!/abs(!Mean_Knn!)"
fieldname = "CV_Knn"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)

#################KT##################
# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Kt,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Kt,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Kt,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Kt,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_Kt', 'Std_Kt', 'Range_Kt', 'Median_Kt']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_Kt'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_Kt'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_Kt'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Median_Kt'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_Kt', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_Kt!/abs(!Mean_Kt!)"
fieldname = "CV_Kt"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)


#################KTN##################
# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Ktn,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Ktn,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Ktn,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Ktn,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_Ktn', 'Std_Ktn', 'Range_Ktn', 'Median_Ktn']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_Ktn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_Ktn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_Ktn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Median_Ktn'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_Ktn', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_Ktn!/abs(!Mean_Ktn!)"
fieldname = "CV_Ktn"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)


#################KTT##################
# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Ktt,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Ktt,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Ktt,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", Ktt,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_Ktt', 'Std_Ktt', 'Range_Ktt', 'Median_Ktt']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_Ktt'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_Ktt'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_Ktt'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Median_Ktt'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_Ktt', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_Ktt!/abs(!Mean_Ktt!)"
fieldname = "CV_Ktt"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)

#################ST##################
# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", St,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", St,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", St,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", St,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_St', 'Std_St', 'Range_St', 'Median_St']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_St'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_St'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_St'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Median_St'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_St', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_St!/abs(!Mean_St!)"
fieldname = "CV_St"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)


#################ISED##################
# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", ISED,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", ISED,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", ISED,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", ISED,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile
fields_to_add = ['Mean_ISED', 'Std_ISED', 'Range_ISED', 'Med_ISED']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_ISED'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_ISED'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_ISED'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Med_ISED'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_ISED', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_ISED!/abs(!Mean_ISED!)"
fieldname = "CV_ISED"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)


#################IGDED##################


# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", IGDED,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", IGDED,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", IGDED,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", IGDED,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_IGDED', 'Std_IGDED', 'Rng_IGDED', 'Med_IGDED']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_IGDED'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_IGDED'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Rng_IGDED'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Med_IGDED'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_IGDED', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_IGDED!/abs(!Mean_IGDED!)"
fieldname = "CV_IGDED"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)


#################PESe##################


# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESe,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESe,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESe,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESe,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_PESe', 'Std_PESe', 'Range_PESe', 'Med_PESe']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_PESe'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_PESe'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_PESe'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Med_PESe'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_PESe', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_PESe!/abs(!Mean_PESe!)"
fieldname = "CV_PESe"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)


#################PESst##################


# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESst,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESst,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESst,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESst,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile for each statistic
fields_to_add = ['Mean_PESst', 'Std_PESst', 'Rng_PESst', 'Med_PESst']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_PESst'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_PESst'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Rng_PESst'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Med_PESst'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_PESst', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_PESst!/abs(!Mean_PESst!)"
fieldname = "CV_PESst"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)


#################PESD##################


# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESD,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESD,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESD,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", PESD,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile
fields_to_add = ['Mean_PESD', 'Std_PESD', 'Range_PESD', 'Med_PESD']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_PESD'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_PESD'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Range_PESD'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Med_PESD'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_PESD', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_PESD!/abs(!Mean_PESD!)"
fieldname = "CV_PESD"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)


#################IPESC##################


# ZonalStatisticsAsTable for Mean, STD, Range and Median
table_mean = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", IPESC,
                                             "table_mean",
                                             "DATA", "MEAN")

table_std = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", IPESC,
                                            "table_std",
                                            "DATA", "STD")

table_range = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", IPESC,
                                              "table_range",
                                              "DATA", "RANGE")

table_median = arcpy.ia.ZonalStatisticsAsTable(segment_layer, "FID", IPESC,
                                               "table_median",
                                               "DATA", "MEDIAN")

# Add new fields to the segments shapefile
fields_to_add = ['Mean_IPESC', 'Std_IPESC', 'Rng_IPESC', 'Med_IPESC']

for field in fields_to_add:
    arcpy.AddField_management(segment_layer, field, 'DOUBLE')

# Create dictionary to store FID and corresponding values
stat_dict = {}

# Populate the dictionary with values from the output tables
stat_dict['Mean_IPESC'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_mean, ['FID', 'MEAN'])}
stat_dict['Std_IPESC'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_std, ['FID', 'STD'])}
stat_dict['Rng_IPESC'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_range, ['FID', 'RANGE'])}
stat_dict['Med_IPESC'] = {row[0]: row[1] for row in arcpy.da.SearchCursor(table_median, ['FID', 'MEDIAN'])}

# Update the fields in the segments shapefile
with arcpy.da.UpdateCursor(segment_layer, ['FID'] + fields_to_add) as cursor:
    for row in cursor:
        fid = row[0]
        for field in fields_to_add:
            if fid in stat_dict[field]:
                row[fields_to_add.index(field) + 1] = stat_dict[field][fid]
        cursor.updateRow(row)

arcpy.AddField_management(segment_layer, 'CV_IPESC', 'DOUBLE')

# Calculate Coefficient of Variation
expression = "!Std_IPESC!/abs(!Mean_IPESC!)"
fieldname = "CV_IPESC"
arcpy.CalculateField_management(segment_layer, fieldname, expression, "PYTHON")

cv_names = ["expression","fieldname"]
intermediate_tables = ["table_mean", "table_std", "table_range", "table_median"]

# Delete intermediate tables and variables for CV
for table in intermediate_tables:
    if arcpy.Exists(table):
        arcpy.management.Delete(table)

for var in cv_names:
    if arcpy.Exists(var):
        arcpy.management.Delete(var)         
#----------------------------------------------------------------------------------------------------------------------------------------------

print("Fields added to the shape file")


####ELEMENTARY FORM TYPE
arcpy.AddField_management(segment_layer, 'Min_CV', 'DOUBLE')
arcpy.AddField_management(segment_layer, 'Min_CV2', 'DOUBLE')
arcpy.AddField_management(segment_layer, 'Min_CV3', 'DOUBLE')
arcpy.AddField_management(segment_layer, 'EF_type', 'TEXT', field_length=3)
arcpy.AddField_management(segment_layer, 'EF_type2', 'TEXT', field_length=3)
arcpy.AddField_management(segment_layer, 'EF_type3', 'TEXT', field_length=3)

# Calculate the minimum CV values
expression_min_cv = "min(!CV_DEM!, !CV_Kn!, !CV_Knn!, !CV_Kt!, !CV_Ktn!, !CV_Ktt!, !CV_sin_asp!, !CV_sin_sl!, !CV_St!)"
fieldname_min_cv = "Min_CV"
arcpy.CalculateField_management(segment_layer, fieldname_min_cv, expression_min_cv, "PYTHON")

expression_min_cv2 = "sorted([!CV_DEM!, !CV_Kn!, !CV_Knn!, !CV_Kt!, !CV_Ktn!, !CV_Ktt!, !CV_sin_asp!, !CV_sin_sl!, !CV_St!])[1]"
fieldname_min_cv2 = "Min_CV2"
arcpy.CalculateField_management(segment_layer, fieldname_min_cv2, expression_min_cv2, "PYTHON")

expression_min_cv3 = "sorted([!CV_DEM!, !CV_Kn!, !CV_Knn!, !CV_Kt!, !CV_Ktn!, !CV_Ktt!, !CV_sin_asp!, !CV_sin_sl!, !CV_St!])[2]"
fieldname_min_cv3 = "Min_CV3"
arcpy.CalculateField_management(segment_layer, fieldname_min_cv3, expression_min_cv3, "PYTHON")

# Create a dictionary to assign names of variable with minimal value
expression_ef_type = "get_ef_type(!Min_CV!, !CV_DEM!, !CV_Kn!, !CV_Knn!, !CV_Kt!, !CV_Ktn!, !CV_Ktt!, !CV_sin_asp!, !CV_sin_sl!, !CV_St!)"
code_block_ef_type = """
def get_ef_type(min_cv, cv_dem, cv_kn, cv_knn, cv_kt, cv_ktn, cv_ktt, cv_sin_asp, cv_sin_sl, cv_st):
    ef_type_mapping = {
        cv_dem: 'Z',
        cv_kn: 'Kn',
        cv_knn: 'Knn',
        cv_kt: 'Kt',
        cv_ktn: 'Ktn',
        cv_ktt: 'Ktt',
        cv_sin_asp: 'A',
        cv_sin_sl: 'S',
        cv_st: 'St'
    }
    return ef_type_mapping[min_cv]
"""

fieldname_ef_type = "EF_type"
arcpy.CalculateField_management(segment_layer, fieldname_ef_type, expression_ef_type, "PYTHON", code_block_ef_type)

expression_ef_type2 = "get_ef_type(!Min_CV2!, !CV_DEM!, !CV_Kn!, !CV_Knn!, !CV_Kt!, !CV_Ktn!, !CV_Ktt!, !CV_sin_asp!, !CV_sin_sl!, !CV_St!)"
fieldname_ef_type2 = "EF_type2"  
arcpy.CalculateField_management(segment_layer, fieldname_ef_type2, expression_ef_type2, "PYTHON", code_block_ef_type)

expression_ef_type3 = "get_ef_type(!Min_CV3!, !CV_DEM!, !CV_Kn!, !CV_Knn!, !CV_Kt!, !CV_Ktn!, !CV_Ktt!, !CV_sin_asp!, !CV_sin_sl!, !CV_St!)"
fieldname_ef_type3 = "EF_type3"
arcpy.CalculateField_management(segment_layer, fieldname_ef_type3, expression_ef_type3, "PYTHON", code_block_ef_type)


##AFFINITY OF LAND SURFACE (LS) SEGMENTS TO VARIABLES

# Computation of ideal elementary forms (IEF) with constant values
# of following variables

arcpy.AddField_management(segment_layer, 'Z', 'DOUBLE')
expression_Z = "!CV_DEM!"
fieldname_Z = "Z"  
arcpy.CalculateField_management(segment_layer, fieldname_Z, expression_Z, "PYTHON")

arcpy.AddField_management(segment_layer, 'S_A', 'DOUBLE')
expression_S_A = "(!CV_sin_sl! + !CV_sin_asp!) / 2"
fieldname_S_A = "S_A"  
arcpy.CalculateField_management(segment_layer, fieldname_S_A, expression_S_A, "PYTHON")

arcpy.AddField_management(segment_layer, 'A_Kn', 'DOUBLE')
expression_A_Kn = "(!CV_sin_asp! + !CV_Kn!) / 2"
fieldname_A_Kn = "A_Kn"  
arcpy.CalculateField_management(segment_layer, fieldname_A_Kn, expression_A_Kn, "PYTHON")

arcpy.AddField_management(segment_layer, 'A_Knn', 'DOUBLE')
expression_A_Knn = "(!CV_sin_asp! + !CV_Knn!) / 2"
fieldname_A_Knn = "A_Knn"  
arcpy.CalculateField_management(segment_layer, fieldname_A_Knn, expression_A_Knn, "PYTHON")

arcpy.AddField_management(segment_layer, 'S_Ktt', 'DOUBLE')
expression_S_Ktt = "(!CV_sin_sl! + !CV_Ktt!) / 2"
fieldname_S_Ktt = "S_Ktt"  
arcpy.CalculateField_management(segment_layer, fieldname_S_Ktt, expression_S_Ktt, "PYTHON")

arcpy.AddField_management(segment_layer, 'S_Ktn', 'DOUBLE')
expression_S_Ktn = "(!CV_sin_sl! + !CV_Ktn!) / 2"
fieldname_S_Ktn = "S_Ktn"  
arcpy.CalculateField_management(segment_layer, fieldname_S_Ktn, expression_S_Ktn, "PYTHON")

arcpy.AddField_management(segment_layer, 'Kn_Ktn', 'DOUBLE')
expression_Kn_Ktn = "(!CV_Kn! + !CV_Ktn!) / 2"
fieldname_Kn_Ktn = "Kn_Ktn"  
arcpy.CalculateField_management(segment_layer, fieldname_Kn_Ktn, expression_Kn_Ktn, "PYTHON")

arcpy.AddField_management(segment_layer, 'Knn_Ktn', 'DOUBLE')
expression_Knn_Ktn = "(!CV_Knn! + !Kn_Ktn!) / 2"
fieldname_Knn_Ktn = "Knn_Ktn"  
arcpy.CalculateField_management(segment_layer, fieldname_Knn_Ktn, expression_Knn_Ktn, "PYTHON")

arcpy.AddField_management(segment_layer, 'St_Kn', 'DOUBLE')
expression_St_Kn = "(!CV_St! + !CV_Kn!) / 2"
fieldname_St_Kn = "St_Kn"  
arcpy.CalculateField_management(segment_layer, fieldname_St_Kn, expression_St_Kn, "PYTHON")

arcpy.AddField_management(segment_layer, 'St_Knn', 'DOUBLE')
expression_St_Knn = "(!CV_St! + !CV_Knn!) / 2"
fieldname_St_Knn = "St_Knn"  
arcpy.CalculateField_management(segment_layer, fieldname_St_Knn, expression_St_Knn, "PYTHON")


#Computation of the highest affinity of LS segments
#to ideal elementary form
arcpy.AddField_management(segment_layer, 'Min_IEF', 'DOUBLE')
arcpy.AddField_management(segment_layer, 'Min_IEF2', 'DOUBLE')
arcpy.AddField_management(segment_layer, 'EF', 'TEXT', field_length=10)
arcpy.AddField_management(segment_layer, 'EF2', 'TEXT', field_length=10)

expression_min_IEF = "min(!Z!, !S_A!, !A_Kn!, !A_Knn!, !S_Ktt!, !S_Ktn!, !Kn_Ktn!, !Knn_Ktn!, !St_Kn!, !St_Knn!)"
fieldname_min_IEF = "Min_IEF"
arcpy.CalculateField_management(segment_layer, fieldname_min_IEF, expression_min_IEF, "PYTHON")

expression_min_IEF2 = "sorted([!Z!, !S_A!, !A_Kn!, !A_Knn!, !S_Ktt!, !S_Ktn!, !Kn_Ktn!, !Knn_Ktn!, !St_Kn!, !St_Knn!])[1]"
fieldname_min_IEF2 = "Min_IEF2"
arcpy.CalculateField_management(segment_layer, fieldname_min_IEF2, expression_min_IEF2, "PYTHON")

expression_ef = "get_ef(!Min_IEF!,!Z!, !S_A!, !A_Kn!, !A_Knn!, !S_Ktt!, !S_Ktn!, !Kn_Ktn!, !Knn_Ktn!, !St_Kn!, !St_Knn!)"
code_block_ef = """
def get_ef(min_ief,z, s_a, a_kn, a_knn, s_ktt, s_ktn, kn_ktn, knn_ktn, st_kn, st_knn):
    ef_mapping = {
        z: 'Z',
        s_a: 'S_A',
        a_kn: 'A_Kn',
        a_knn: 'A_Knn',
        s_ktt: 'S_Ktt',
        s_ktn: 'S_Ktn',
        kn_ktn: 'Kn_Ktn',
        knn_ktn: 'Knn_Ktn',
        st_kn: 'St_Kn',
        st_knn: 'St_Knn'
    }
    return ef_mapping[min_ief]
"""

fieldname_ef = "EF"
arcpy.CalculateField_management(segment_layer, fieldname_ef, expression_ef, "PYTHON", code_block_ef)

#Find affinity of IEF to each segment by the formula:
# 1 - (MinS1 / 2) and assign the value to each segment
#in column "Affinity"
arcpy.AddField_management(segment_layer, 'Affinity', 'DOUBLE')
expression_affinity = "1 - (!Min_IEF! / 2)"
fieldname_affinity = "Affinity"
arcpy.CalculateField_management(segment_layer, fieldname_affinity, expression_affinity, "PYTHON")

expression_ef2 = "get_ef(!Min_IEF2!,!Z!, !S_A!, !A_Kn!, !A_Knn!, !S_Ktt!, !S_Ktn!, !Kn_Ktn!, !Knn_Ktn!, !St_Kn!, !St_Knn!)"
fieldname_ef2 = "EF2"  
arcpy.CalculateField_management(segment_layer, fieldname_ef2, expression_ef2, "PYTHON", code_block_ef)

arcpy.AddField_management(segment_layer, 'Affinity2', 'DOUBLE')
expression_affinity2 = "1 - (!Min_IEF2! / 2)"
fieldname_affinity2 = "Affinity2"
arcpy.CalculateField_management(segment_layer, fieldname_affinity2, expression_affinity2, "PYTHON")


output_excel = r"D:\Dissertation\Toolbox_Physical_Geomorphometry\outputs\output_table.xlsx"

# Use TableToExcel to export the attribute table to Excel
arcpy.conversion.TableToExcel(segment_layer, output_excel)
