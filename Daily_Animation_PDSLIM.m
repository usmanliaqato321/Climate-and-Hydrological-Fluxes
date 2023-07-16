clear all
tic
 addpath('F:\NALTAR\PDSLIM_Karakoram_2_15\NETCDF_create') %function folder
% cd 'F:/NALTAR/PDSLIM_Karakoram_2_13' %data folder

% import mask
maschera = importdata('Mask_New.asc');
maschera = maschera';
maschera(maschera == 0) = NaN;

% prepare the HCONVXY 3d matrix
file = 'SWE.dat'
columns = 260;
binformat = 'real*4' %check the fread function 
tsteps = 182; %number of timesteps
SWE = bin2arr(file, binformat,columns,tsteps,maschera);
%HCONVXY = int16(HCONVXY);
% 40269 for 01-04-2010
% 41729 for 01-04-2014
% 54878 for 01-04-2050
% 69488 for 01-04-2090
%CREATE TIME VECTOR
startdate = 54878; %check the date using datetime(1900,1,1,startdate,0,0) [hours] datetime(1900,1,changehere)
datetime(1900,1,startdate)
dt = 1; %in horurs
timetot = 182; %24*ndays

for i = 1:timetot
    time(i) = startdate+i-1;
end

longitude = importdata('LONGITUDE.txt');
latitude = importdata('LATITUDE.txt');

longitude = single(longitude);
latitude = single(latitude);
time = int32(time');

%creation of netcdf file
nccreate('Karakoram_test.nc','SWE',"Dimensions",{'longitude',size(longitude,1),'latitude',size(latitude,1),'time',size(time,1)},"Datatype","double","Format", "64bit");
nccreate('Karakoram_test.nc','longitude',"Dimensions",{'longitude',size(longitude,1)},"Datatype","single","Format", "64bit");
nccreate('Karakoram_test.nc','latitude',"Dimensions",{'latitude',size(latitude,1)},"Datatype","single","Format", "64bit");
nccreate('Karakoram_test.nc','time',"Dimensions",{'time',size(time,1)},"Datatype","int32","Format", "64bit");

ncwrite("Karakoram_test.nc","SWE",SWE);
ncwrite("Karakoram_test.nc","time",time);
ncwrite("Karakoram_test.nc","longitude",longitude);
ncwrite("Karakoram_test.nc","latitude",latitude);

toc
%% WRITE ATTRIBUTES
ncwriteatt("Karakoram_test.nc","time","units","days since 1900-01-01")