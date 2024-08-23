% The input, R, is a map.rasterref.MapCellsReference object indicating that you are working in a
% projected coordinate system. If so, then specify a projected coordinate system by setting the
% appropriate values for the 'CoordRefSysCode' or 'GeoKeyDirectoryTag' optional parameters.

info = geotiffinfo(([root_path,'/BedTopo/BedmachineV4_BedTopography.tif']));
geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
[A,R] = readgeoraster([root_path,'/BedTopo/BedmachineV4_BedTopography.tif']);
geotiffwrite('bed_error.tif', bed_err, R,'GeoKeyDirectoryTag',geoTags);
