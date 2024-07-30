function [DEMstruct]= ThicknessEstimates(old_thickness,new_thickness,aero_err,arctic_err,FrontalAB, Poly)
% Function determines Surface elevation change per year based on two
% DEMS, and calculates ice thickness using surface elvation and bedrock 
% topography (automatically cropped to DEM extent)for each available 
% date of terminus positions. 
%
% arguments: (input)
% PathToLatestDEM - Path to directory of latest DEM. DEM has to be in Polar 
%           Stereographic projection (EPSG:3413).
% PathToOldestDEM - Path to directory of oldest available DEM.DEM has to be 
%           in Polar Stereographic projection (EPSG:3413).
% PathToBedrockTopo - Path the Bedmachine v4 data set.
% FrontalAB - Matrix containing area of ROIs (not used here) and dates of 
%           digitised termini
%
% Input raster data is automatically downsampled to the resolution of the
% bedrock topography data (150 m)
%
% arguments: (output)
% DEM.X - X coordinates of corresponding to caluclated ice thickness values
% DEM.Y - Y coordinates of corresponding to caluclated ice thickness values     
% DEM.H - Ice thickness values calculated assuming linear change in surface
%         elevation. Annual change of surface elevation is determined by
%         differencing the two DEMs and dividing values by number of years
%         between both. Surface elevation is heightend by adding the annual 
%         elevation change multiplied by the elapsed time between terminus 
%         observations for years earlier than the creation date of the
%         oldest DEM. Surface elevation is lowered by subtracting the annual 
%         elevation change multiplied by the elapsed time between terminus 
%         observations for years later than the creation date of the
%         oldest DEM. 
%
%
% NOTE: These calculations assume linear changes in surface elevation and
% therefore do not account for annual dynamic variability
%         

% PathToLatestDEM = 'Data/DEMs/Narsap_Sermia/ArcticDEM_2017_Narsap.tif';
% PathToOldestDEM = 'Data/DEMs/Narsap_Sermia/1985_DEM_NS_Clipped.tif';
% PathToBedrockTopo = 'Data/BedTopo/BedmachineV4_BedTopography.tif';

    [Old_thick,Rold] = geotiffread(old_thickness);
    [Latest_thick,Rlatest] =geotiffread(new_thickness);
    [aero_err,~] = geotiffread(aero_err);
    [arctic_err,~] = geotiffread(arctic_err);

    %% Find acquisition year of earliest and latest DEM and calculate yearly change
  
   
    k = strfind(old_thickness,'19');
    EarliestDEMYear = str2num(old_thickness(k:k+3));
    l = strfind(new_thickness,'20');
    LatestDEMYear = str2num(new_thickness(l:l+3));
    DEMYearDifference = LatestDEMYear-EarliestDEMYear;
    [x,y] = pixcenters(Rold,size(Old_thick));
    [xDEMfinal,yDEMfinal] = meshgrid(x,y);
    DEM{1}.X = xDEMfinal;
    DEM{1}.Y = yDEMfinal;
    l1 = length(FrontalAB);
    ROImask{1} = inpolygon(DEM{1}.X,DEM{1}.Y,Poly{1}.Vertices(:,1),Poly{1}.Vertices(:,2));
    ROImask{2} = inpolygon(DEM{1}.X,DEM{1}.Y,Poly{l1}.Vertices(:,1),Poly{l1}.Vertices(:,2));
    Old_thick(Old_thick>1e+38) = NaN;
    Old_thick(Old_thick>1e+38) = NaN;
    Latest_thick(Latest_thick<-1e+38) = NaN;
    Latest_thick(Latest_thick<-1e+38) = NaN;
    aero_err(aero_err>1e+38 | aero_err<-1e+38) = NaN;
    arctic_err(arctic_err>1e+38 | arctic_err<-1e+38) = NaN;
    DiffDEM = nanmean(Old_thick(ROImask{2}),'all') - nanmean(Latest_thick(ROImask{2}),"all");

    DiffDEM(DiffDEM<-1e+38)=NaN;
    YearlyDiff = DiffDEM/DEMYearDifference;
    yearlyChange = YearlyDiff;
    dailyChange = yearlyChange./365;
    Change_new = DEMYearDifference*abs(yearlyChange);
    old_resampled = Latest_thick + Change_new;

    %% Get Dates of observations of terminus positions
    Tdates = datetime(single(FrontalAB(:,1)),'ConvertFrom','datenum');

    TdateLatestDEM = datetime(new_thickness(l:l+3),'Format','yyyy');

    DifferenceToLatestDEM.Hours = TdateLatestDEM-Tdates;


    %%
    % Subtract/add surface elevation change form 1985 DEM based ontimedelta
    % between Terminus position observations and get Ice Thickness (H) by
    % subtracting surface elvation from bedrock elevation
    for i=1:length(FrontalAB)
         [x,y] = pixcenters(Rlatest,size(Latest_thick));
         [xDEMfinal,yDEMfinal] = meshgrid(x,y);
            if year(Tdates(i)) <= EarliestDEMYear
                DEM{i}.X = xDEMfinal;
                DEM{i}.Y = yDEMfinal;
                fprintf('Date is the same as Earliest DEM \n')
                fprintf('No change to surface height \n')
%                 Change{i} = DEMYearDifference*abs(yearlyChange);
%                 NewSurfaceHeight{i} = Latest_thick + Change{i}; 
                NewSurfaceHeight{i} = Old_thick;
                DEM{i}.H = NewSurfaceHeight{i};
                DEM{i}.U1 = aero_err;DEM{i}.U2 = aero_err;
                DEM{i}.L1 = aero_err;DEM{i}.L2 = aero_err;
            elseif year(Tdates(i)) > EarliestDEMYear && year(Tdates(i)) < LatestDEMYear 
                DEM{i}.X = xDEMfinal;
                DEM{i}.Y = yDEMfinal;
                timeDelta(i) = LatestDEMYear-year(Tdates(i)); 
                timeDelta_upp(i) = year(Tdates(i)) - EarliestDEMYear;
                fprintf('Date before Latest DEM \n')
                fprintf('Adding surface height \n') 
                Change{i} = timeDelta(i)*abs(yearlyChange);
                NewSurfaceHeight{i} = Latest_thick + Change{i}; %gave much more realistic values
                DEM{i}.H = NewSurfaceHeight{i}; 
                DEM{i}.L1 = Old_thick - aero_err; DEM{i}.L2 = Latest_thick - arctic_err; 
                DEM{i}.U1 = Old_thick + aero_err; DEM{i}.U2 = Latest_thick + arctic_err;  
            elseif year(Tdates(i)) >= LatestDEMYear 
                DEM{i}.X = xDEMfinal;
                DEM{i}.Y = yDEMfinal;
                timeDelta(i) = 1;
                timeDelta_upp(i) = 1;
                fprintf('Date is the same as Latest DEM \n')
                fprintf('No change to surface height \n')
                Change{i} = DEMYearDifference*abs(yearlyChange);
                NewSurfaceHeight{i} = Latest_thick; %
                DEM{i}.H = NewSurfaceHeight{i}; 
                DEM{i}.U1 = arctic_err; DEM{i}.U2 = arctic_err;
                DEM{i}.L1 = arctic_err; DEM{i}.L2 = arctic_err;
            end
     DEMstruct = DEM;
%      ChangeStruct = Change;

    end
