function [alon,alat,adata]=cropregion(lon,lat,data,bbox)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%maskregion is a part of Atmospheric Science ToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by Ankur Kumar                        Wednesday; June 09, 2019
%                                              Version: 1.0
%
% PURDUE UNIVERISTY
% West Lafayette, IN 47907
% Department of Earth, Atmospheric, and Planetary Sciences
% Email: ankurk017@gmail.com
%        kumar409@purdue.edu
%
% 
% Function:
%           cropregion crops the data in the nc file as per the bounding
%           box provided as the fourth argument
% Syntax:
%           [lon_cropped, lat_cropped, data_cropped]=cropregion(lon,lat,data,[lon1 lon2 lat1 lat2]);
% 
% Inputs:
%           It takes the first three arguments as the latitude, longitude, and
%           data respectively. The fourth argument should be the bounding box of the region(s) you wish to crop.
% 
% 
% Example:
% 
%            W=shaperead('us_states.shp');
%           [xx,yy,zz]=cropregion(lon,lat,data,[-79.5683  -74.7735 33.5789   36.7302])
%
% Please send your suggestions to the email id: ankurk017@gmail.com or
%                                               416AS2025@nitrkl.ac.in
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ncdispread is a part of Atmospheric Science ToolBox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1=find(lon>bbox(1) & lon <bbox(2));
id2=find(lat>bbox(3) & lat <bbox(4));
alon=lon(id1);
alat=lat(id2);
if length(lon)~=size(data,1)
    adata=data(id2,id1,:);
else
    adata=data(id1,id2,:);
end
end
