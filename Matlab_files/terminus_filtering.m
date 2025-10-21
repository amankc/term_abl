clear
close all;

%% 1. Reads the terminus shapefiles and Centerline data
root_path = '/Users/amankc/Terminus_Ablation/';
cd(root_path)
glacier_name = 'Rimfaxe';
region = 'NE';
glacierID = 31; vel = 21.92; % glacier speed in m/day
cenline_shp = append(glacier_name,'.shp');
glacier_shp = append('Combined_',glacier_name,'.shp');
centerline_path = fullfile(root_path,glacier_name,'Centerline',['Cenline_',cenline_shp]);
cline = shaperead(centerline_path);

% Process Centerline
x = rmmissing(cline.X');
y = rmmissing(cline.Y');
full_length = sqrt((x(end)-x(1)).^2 + (y(end)-y(1)).^2);
dx = 2; % ideal spacing in meters
inds = round(full_length/dx);
[points] = interparc(inds,x,y,'linear');
points = rmmissing(points);

% Calculate Centerline Distance from Origin (upstream end)
centerline_dist(1) = 0;
for k = 2:length(points)
    centerline_dist(k) = centerline_dist(k-1)+ sqrt((points(k,1)-points(k-1,1)).^2 + (points(k,2)-points(k-1,2)).^2);
end

% Read Terminus Positions
S = shaperead(fullfile(root_path,glacier_name,'Terminus_Positions',glacier_shp));

% Assign DecDate and Sort 
for i = 1:length(S)
    if isnan(S(i).DecDate)
        S(i).DecDate = decyear(datetime(S(i).Date,'InputFormat','yyyy-MM-dd'));
    end
end
[~,index] = sortrows([S.DecDate].'); 
S = S(index); 

%% 2. Quality control for multiple traces on a same day
Save = {};
mm=1;
i=1;
while i<=length(S)
    n = 2;
    temp = S(i).DecDate;
    skip = 1;
    index = zeros(0);
    for j = i+1:length(S)
        if temp == S(j).DecDate
            index(1) = S(i).flag;
            index(n) = S(j).flag;
            n = n+1;
        end
    end
    if length(index)>1
        [~,id] = min(index);
        temp2 = S(i+id);
        skip = length(index);

    else
        temp2 = S(i);
    end
    Save {mm} = temp2;
    mm=mm+1;
    i = i + skip;
end

S_filt = cell2mat(Save);
indss = [];
con = 1;
for i=1:length(S_filt)-1
    if S_filt(i).DecDate == S_filt(i+1).DecDate
        indss(con) = i;
        con = con+1;
    end
end
%% Getting rid of the repeated data
% RUN ONLY ONCE!
S_filt(indss) = [];

%%
% Convert S_filt array to cell array
S_fil = cell(1, length(S_filt));
for i = 1:length(S_filt)
    S_fil{i} = S_filt(i);
end

%% 9. Centerline Intersections and Distances (SETUP for velocity filter)
centerdist = zeros(1, length(S_fil));
termdate = zeros(1, length(S_fil));
for k = 1:length(S_fil)
    decdate(k) = S_fil{1,k}.DecDate;
    termdate(k) = datenum(S_fil{1,k}.Date);
%     S(i).X = rmmissing(S(i).X); S(i).Y = rmmissing(S(i).Y); 
    [xi,yi,ii] = polyxpoly(points(:,1),points(:,2),S_fil{1,k}.X',S_fil{1,k}.Y'); %xi, yi are intersection coordinates
    if ~isempty(ii) 
        centerdist(k) = term(ii).center_dist;
    else
        centerdist(k) = NaN;
    end
    clear ii;
end
termdate(isnan(centerdist)) = []; centerdist(isnan(centerdist)) = []; 

% Remove NaNs (where terminus didn't intersect centerline)
non_nan_idx = ~isnan(centerdist);
termdate = termdate(non_nan_idx);
centerdist = centerdist(non_nan_idx);
S_fil = S_fil(non_nan_idx);


%% 10. Filtering using ITS_LIVE Velocity (CORRECTED)

% Anomalous: rate > rate_thresh (fast retreat/spike)
% Anomalous: rate < -rate_thresh (fast advance/dip)

threshold = 2;
rate_thresh = threshold * vel; % 2 * 21.92 m/day
Shp_filtered = S_fil;
ind = []; % indices to remove
counter = 1;

for k = 1:length(centerdist) - 2
    % Calculate the rate (change in distance / change in time) from k to k+1
    % Positive rate = Retreat; Negative rate = Advance
    rate = (centerdist(k + 1) - centerdist(k)) / (termdate(k + 1) - termdate(k));
    
    % Anomalous Retreat
    if rate > rate_thresh
        ind(counter) = k + 1; % Flag the outlier position
        counter = counter + 1;
    end
    
    % Anomalous Advance
    % Directly check for an extremely fast advance (highly negative rate)
    if rate < -rate_thresh 
        ind(counter) = k + 1;
        counter = counter + 1;
    % Check for the classic dip-and-spike pattern: small advance followed by massive retreat
    elseif rate < 0 % If the terminus has advanced (or small retreat)
       
        rate_next = (centerdist(k + 2) - centerdist(k + 1)) / (termdate(k + 2) - termdate(k + 1));
        
        % If the rate from k+1 to k+2 is an anomalously fast retreat (the spike after the dip)
        if rate_next > rate_thresh
            ind(counter) = k + 1; % Flag k+1, which is the dip
            counter = counter + 1;
        end
    end
end

% Remove unique flagged indices
Shp_filtered(:,unique(ind)) = [];

for i=1:length(Shp_filtered)
    plot(Shp_filtered{1,i}.X,Shp_filtered{1,i}.Y, 'DisplayName', datestr(Shp_filtered{1,i}.Date))
    hold on;
end
title('Filtered Terminus Positions')
xlabel('Easting (m)')
ylabel('Northing (m)')
axis equal
hold off;

%% 13. Clearing the length (OPTIONAL filter)
% Note: Skipping all steps for Fjord Boundaries, BLN, and Midpoint velocity as requested.
% This length-based filter is a common final quality check.
try
    px = zeros(1, length(Shp_filtered));
    for i=1:length(Shp_filtered)
        x_coord = rmmissing(Shp_filtered{1,i}.X);
        y_coord = rmmissing(Shp_filtered{1,i}.Y);
        d = diff([x_coord(:) y_coord(:)]);
        total_length = sum(sqrt(sum(d.*d,2)));
        px(i) = total_length/1000;
    end
    
    ml= median(px); ctr=1;
    iidx_len = []; 
    for i=1:length(px)
        if px(i) > 1.5 * ml || px(i) < 0.5 * ml
            iidx_len(ctr) = i;
            ctr = ctr+1;
        end
    end
    
    if isempty(iidx_len)
        disp('Length filter: No points found that require removal.')
    else
        Shp_filtered(:,unique(iidx_len)) = [];
        disp(['Length filter: Removed ', num2str(length(unique(iidx_len))), ' traces based on length anomaly.'])
    end

catch 
    % Catch if Shp_filtered is empty after velocity filtering
    disp('Skipping length filter: No data to process or an error occurred.')
end

%% Final plotting

for i=1:length(Shp_filtered)
    plot(Shp_filtered{1,i}.X,Shp_filtered{1,i}.Y)
    hold on;
end
title('Final Filtered Terminus Positions (After Velocity and Length Checks)')
xlabel('Easting (m)')
ylabel('Northing (m)')
axis equal