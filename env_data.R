##Loading environmental data from Copernicus Marine##

library(ncdf4)
library(lubridate)

#Loading in dissolved oxygen monthly mean product for October 2022

oxygen <- nc_open("mercatorbiomer4v2r1_global_mean_bio_202210.nc")
print(oxygen) #gives metadata

#See the structure of the file
names(oxygen)
names(oxygen$var) #o2 is what we want
names(oxygen$dim) #time, longitude, latitude, depth

#Extract the coordinates
dim_lon <- ncvar_get(oxygen, "longitude") #lon and lat are by 0.25 of a degree
dim_lat <- ncvar_get(oxygen, "latitude")
dim_depth <- ncvar_get(oxygen, "depth") #50 depth options
dim_time <- ncvar_get(oxygen, "time") #only one time: 638076

#Make coordinates into data  frames for easier viewing
dim_lon <- as.data.frame(dim_lon)
dim_lat <- as.data.frame(dim_lat)
dim_depth <- as.data.frame(dim_depth)

#Time conversion
t_units <- ncatt_get(oxygen, "time", "units")
t_ustr <- strsplit(t_units$value, " ")
t_dstr <- strsplit(unlist(t_ustr)[3], "-")
date <- ymd(t_dstr) + dseconds(dim_time)
date

#Extract the oxygen variable
o2 <- ncvar_get(oxygen, "o2", collapse_degen=FALSE)

#Extract the o2 value for each sampling site
#Format float o2[longitude,latitude,depth,time]
ga1_o2 <- o2[396, 446, 6, 1] #1st GA site ind 1-26 lon=-81.25, lat=31.25, depth=6.44
ga2_o2 <- o2[396, 445, 7, 1] #2nd GA site ind 27-35 lon=-81.25, lat=31, depth=7.93

ob1_o2 <- o2[410, 457, 8, 1] #1st OB site ind 36-50 lon=-77.75, lat=34, depth=9.57
ob2_o2 <- o2[410, 458, 8, 1] #2nd OB site ind 51-60 lon=-77.75, lat=34.25, depth=9.57
ob3_o2 <- o2[413, 459, 7, 1] #3rd OB site ind 61-70 lon=-77, lat=34.5, depth=7.93

lb1_o2 <- o2[408, 456, 8, 1] #1st and 2nd LB site ind 71-102 lon=-78.25, lat=33.75, depth=9.57

fl_o2 <- o2[396, 442, 8, 1] #FL site ind 103-137 lon=-81.25, lat=30.25, depth=9.57

sc_o2 <- o2[403, 452, 7, 1] #SC site ind 138-172 lon=-79.5, lat=32.75, depth=7.93

de1_o2 <- o2[425, 481, 11, 1] #DE site 1 ind 173 lon=-74, lat=40, depth=15.81
de2_o2 <- o2[424, 479, 11, 1] #DE site 2 ind 185-188 lon=-74.25, lat=39.5, depth=15.81
de3_o2 <- o2[422, 475, 11, 1] #DE site 3 ind 219-241 lon=-74.75, lat=38.5, depth=15.81

md1_o2 <- o2[420, 471, 12, 1] #MD site 1 and ind 189-211 and 218 lon=-75.25, lat=37.5, depth=18.5
md2_o2 <- o2[422, 474, 12, 1] #MD site 2 ind 242-261 lon=-74.75, lat=38.25, depth=18.5
md3_o2 <- o2[421, 473, 12, 1] #MD site 3 ind 262-284 lon=-75, lat=38, depth=18.5

va1_o2 <- o2[418, 467, 11, 1] #VA site 1 ind 174-178 and 180-184 lon=-75.75, lat=36.75, depth=15.81
va2_o2 <- o2[418, 467, 8, 1] #VA site 2 ind. 179 lon=-75.75, lat=36.5, depth=9.57
va3_o2 <- o2[419, 470, 12, 1] #VA site 3  ind. 285-311 lon=-75.5, lat=37.25, depth=18.5
va4_o2 <- o2[417, 469, 9, 1] #VA site 4  ind. 312-351 lon=-76, lat=37, depth=11.41

nc1_o2 <- o2[418, 467, 9, 1] #NC site 1 ind. 212-217 lon=-75.75, lat=36.5, depth=11.41
nc2_o2 <- o2[419, 465, 11, 1] #NC site 1 ind. 351-377 lon=-75.5, lat=36, depth=15.81
nc3_o2 <- o2[419, 462, 9, 1] #NC site 1 ind. 378-400 lon=-75.5, lat=35.25, depth=15.81

#Add the extracted o2 values to the env_oct data frame
env_oct$oxygen <- "NA"

env_oct[1:26, 10] <- ga1_o2
env_oct[27:35, 10] <- ga2_o2
env_oct[36:50, 10] <- ob1_o2
env_oct[51:60, 10] <- ob2_o2
env_oct[61:70, 10] <- ob3_o2
env_oct[71:102, 10] <- lb1_o2
env_oct[103:137, 10] <- fl_o2
env_oct[138:172, 10] <- sc_o2
env_oct[173, 10] <- de1_o2
env_oct[185:188, 10] <- de2_o2
env_oct[219:241, 10] <- de3_o2
env_oct[189:211, 10] <- md1_o2
env_oct[218, 10] <- md1_o2
env_oct[242:261, 10] <- md2_o2
env_oct[262:284, 10] <- md3_o2
env_oct[174:178, 10] <- va1_o2
env_oct[180:184, 10] <- va1_o2
env_oct[179, 10] <- va2_o2
env_oct[285:311, 10] <- va3_o2
env_oct[312:351, 10] <- va4_o2
env_oct[212:217, 10] <- nc1_o2
env_oct[352:377, 10] <- nc2_o2
env_oct[378:400, 10] <- nc3_o2

#Loading in temperature and salinity monthly mean product for October 2022

temp_sal <- nc_open("glo12_rg_1m-m_202210-202210_2D_hcst.nc")
print(temp_sal) #gives metadata

#See the structure of the file
names(temp_sal)
names(temp_sal$var) #sob is sea bottom salinity in PSU and tob is sea bottom temp
#[longitude, latitude, time]
names(temp_sal$dim) #time, longitude, latitude, depth (no depth for sob or tob)

#Extract the coordinates
dim_lon_ts <- ncvar_get(temp_sal, "longitude") #lon and lat are by 1/12 of a degree
dim_lat_ts <- ncvar_get(temp_sal, "latitude")
dim_depth_ts <- ncvar_get(temp_sal, "depth")
dim_time_ts <- ncvar_get(temp_sal, "time") #only one time: 638076

#Make coordinates into data  frames for easier viewing
dim_lon_ts <- as.data.frame(dim_lon_ts)
dim_lat_ts <- as.data.frame(dim_lat_ts)

#Extract the temp and salinity variables
temp <- ncvar_get(temp_sal, "tob", collapse_degen=FALSE)
sal <- ncvar_get(temp_sal, "sob", collapse_degen=FALSE)

#Extract the temp value for each sampling site
#Format float tob[longitude,latitude,time]
ga1_temp <- temp[1186, 1335, 1] #1st GA site ind 1-26 lon=-81.25, lat=31.167
ga2_temp <- temp[1185, 1333, 1] #2nd GA site ind 27-35 lon=-81.33, lat=31

ob1_temp <- temp[1227, 1370, 1] #1st OB site ind 36-50 lon=-77.83, lat=34.08
ob2_temp <- temp[1229, 1373, 1] #2nd OB site ind 51-60 lon=-77.66, lat=34.33
ob3_temp <- temp[1237, 1377, 1] #3rd OB site ind 61-70 lon=-77, lat=34.66

lb1_temp <- temp[1222, 1368, 1] #1st and 2nd LB site ind 71-102 lon=-78.25, lat=33.92

fl_temp <- temp[1185, 1324, 1] #FL site ind 103-137 lon=-81.33, lat=30.25

sc_temp <- temp[1207, 1356, 1] #SC site ind 138-172 lon=-79.5, lat=32.92

de1_temp <- temp[1273, 1442, 1] #DE site 1 ind 173 lon=-74.0, lat=40.08
de2_temp <- temp[1270, 1435, 1] #DE site 2 ind 185-188 lon=-74.25, lat=39.5
de3_temp <- temp[1263, 1423, 1] #DE site 3 ind 219-241 lon=-74.83, lat=38.5

md1_temp <- temp[1256, 1411, 1] #MD site 1 and ind 189-211 and 218 lon=-75.41, lat=37.5
md2_temp <- temp[1262, 1422, 1] #MD site 2 ind 242-261 lon=-74.91, lat=38.41
md3_temp <- temp[1260, 1417, 1] #MD site 3 ind 262-284 lon=-75.08, lat=38

va1_temp <- temp[1251, 1402, 1] #VA site 1 ind 174-178 and 180-184 lon=-75.83, lat=36.75
va2_temp <- temp[1251, 1400, 1] #VA site 2 ind. 179 lon=-75.83, lat=36.58
va3_temp <- temp[1254, 1411, 1] #VA site 3  ind. 285-311 lon=-75.58, lat=37.5
va4_temp <- temp[1251, 1404, 1] #VA site 4  ind. 312-351 lon=-75.83, lat=36.91

nc1_temp <- temp[1253, 1395, 1] #NC site 1 ind. 212-217 lon=-75.67, lat=36.17
nc2_temp <- temp[1255, 1393, 1] #NC site 1 ind. 352-377 lon=-75.5, lat=36
nc3_temp <- temp[1255, 1384, 1] #NC site 1 ind. 378-400 lon=-75.5, lat=35.25

#Add the extracted o2 values to the env_oct data frame
env_oct$temp <- "NA"

env_oct[1:26, 11] <- ga1_temp
env_oct[27:35, 11] <- ga2_temp
env_oct[36:50, 11] <- ob1_temp
env_oct[51:60, 11] <- ob2_temp
env_oct[61:70, 11] <- ob3_temp
env_oct[71:102, 11] <- lb1_temp
env_oct[103:137, 11] <- fl_temp
env_oct[138:172, 11] <- sc_temp
env_oct[173, 11] <- de1_temp
env_oct[185:188, 11] <- de2_temp
env_oct[219:241, 11] <- de3_temp
env_oct[189:211, 11] <- md1_temp
env_oct[218, 11] <- md1_temp
env_oct[242:261, 11] <- md2_temp
env_oct[262:284, 11] <- md3_temp
env_oct[174:178, 11] <- va1_temp
env_oct[180:184, 11] <- va1_temp
env_oct[179, 11] <- va2_temp
env_oct[285:311, 11] <- va3_temp
env_oct[312:351, 11] <- va4_temp
env_oct[212:217, 11] <- nc1_temp
env_oct[352:377, 11] <- nc2_temp
env_oct[378:400, 11] <- nc3_temp

#Extract the salinity value for each sampling site
#Format float sob[longitude,latitude,time]
ga1_sal <- sal[1186, 1335, 1] #1st GA site ind 1-26 lon=-81.25, lat=31.167
ga2_sal <- sal[1185, 1333, 1] #2nd GA site ind 27-35 lon=-81.33, lat=31

ob1_sal <- sal[1227, 1370, 1] #1st OB site ind 36-50 lon=-77.83, lat=34.08
ob2_sal <- sal[1229, 1373, 1] #2nd OB site ind 51-60 lon=-77.66, lat=34.33
ob3_sal <- sal[1237, 1377, 1] #3rd OB site ind 61-70 lon=-77, lat=34.66

lb1_sal <- sal[1222, 1368, 1] #1st and 2nd LB site ind 71-102 lon=-78.25, lat=33.92

fl_sal <- sal[1185, 1324, 1] #FL site ind 103-137 lon=-81.33, lat=30.25

sc_sal <- sal[1207, 1356, 1] #SC site ind 138-172 lon=-79.5, lat=32.92

de1_sal <- sal[1273, 1442, 1] #DE site 1 ind 173 lon=-74.0, lat=40.08
de2_sal <- sal[1270, 1435, 1] #DE site 2 ind 185-188 lon=-74.25, lat=39.5
de3_sal <- sal[1263, 1423, 1] #DE site 3 ind 219-241 lon=-74.83, lat=38.5

md1_sal <- sal[1256, 1411, 1] #MD site 1 and ind 189-211 and 218 lon=-75.41, lat=37.5
md2_sal <- sal[1262, 1422, 1] #MD site 2 ind 242-261 lon=-74.91, lat=38.41
md3_sal <- sal[1260, 1417, 1] #MD site 3 ind 262-284 lon=-75.08, lat=38

va1_sal <- sal[1251, 1402, 1] #VA site 1 ind 174-178 and 180-184 lon=-75.83, lat=36.75
va2_sal <- sal[1251, 1400, 1] #VA site 2 ind. 179 lon=-75.83, lat=36.58
va3_sal <- sal[1254, 1411, 1] #VA site 3  ind. 285-311 lon=-75.58, lat=37.5
va4_sal <- sal[1251, 1404, 1] #VA site 4  ind. 312-351 lon=-75.83, lat=36.91

nc1_sal <- sal[1253, 1395, 1] #NC site 1 ind. 212-217 lon=-75.67, lat=36.17
nc2_sal <- sal[1255, 1393, 1] #NC site 1 ind. 352-377 lon=-75.5, lat=36
nc3_sal <- sal[1255, 1384, 1] #NC site 1 ind. 378-400 lon=-75.5, lat=35.25

#Add the extracted o2 values to the env_oct data frame
env_oct$salinity <- "NA"

env_oct[1:26, 12] <- ga1_sal
env_oct[27:35, 12] <- ga2_sal
env_oct[36:50, 12] <- ob1_sal
env_oct[51:60, 12] <- ob2_sal
env_oct[61:70, 12] <- ob3_sal
env_oct[71:102, 12] <- lb1_sal
env_oct[103:137, 12] <- fl_sal
env_oct[138:172, 12] <- sc_sal
env_oct[173, 12] <- de1_sal
env_oct[185:188, 12] <- de2_sal
env_oct[219:241, 12] <- de3_sal
env_oct[189:211, 12] <- md1_sal
env_oct[218, 12] <- md1_sal
env_oct[242:261, 12] <- md2_sal
env_oct[262:284, 12] <- md3_sal
env_oct[174:178, 12] <- va1_sal
env_oct[180:184, 12] <- va1_sal
env_oct[179, 12] <- va2_sal
env_oct[285:311, 12] <- va3_sal
env_oct[312:351, 12] <- va4_sal
env_oct[212:217, 12] <- nc1_sal
env_oct[352:377, 12] <- nc2_sal
env_oct[378:400, 12] <- nc3_sal
