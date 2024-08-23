clear
close all;
%% 1. Reads the terminus shapefiles from the directory
root_path = '/Users/amankc/Terminus_Ablation/';
cd(root_path)
glacier_name = 'Rimfaxe';
region = 'NE';
glacierID = 31; vel = 21.92; %in m/day
cenline_shp = append(glacier_name,'.shp');
glacier_shp = append('Combined_',glacier_name,'.shp');
centerline_path = fullfile(root_path,glacier_name,'Centerline',['Cenline_',cenline_shp]);
cline = shaperead(centerline_path);

% Run this if you are getting input as a line feature
x = cline.X';
y = cline.Y';
x = rmmissing(x);%removes NaN values
y = rmmissing(y);%removes NaN values

%determine ideal number of points for semi-standard spacing across all glaciers
full_length = sqrt((x(end)-x(1)).^2 + (y(end)-y(1)).^2);
dx = 2; %ideal spacing in meters
inds = round(full_length/dx);

%Apply the function to get evenly sapced points
[pt,dudt,fofthandle] = interparc(inds,x,y,'linear');
points = rmmissing(pt);

%convert the centerline coordinates to along-profile distance from the origin
term(1).center_dist = 0;
for k = 2:length(points)
    term(k).center_dist = term(k-1).center_dist+ sqrt((points(k,1)-points(k-1,1)).^2 + (points(k,2)-points(k-1,2)).^2);
end

%instead of loading centerline intersections, load full terminus traces
terminuspos_path = fullfile(root_path,glacier_name,'Terminus_Positions',glacier_shp);
S = shaperead(terminuspos_path);

%Assigning quality flag to the authors
names = {'Brough';'Bjork';'Black';'Black_Taryn';'Bunce';'Catania';'Fahrner';'Hill';'Korsgaard';'Moon_Twila';'Sole';'Termpicks';'Wood';'Zhang';
    'Carr';'Murray';'Bevan';'Cheng';'Cheng_D';'ESA';'PROMICE'};
qual_flag = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,3,3,3,3,3];

for j = 1:length(S)
    match = strcmp(S(j).Author,names);
    if ~isempty(qual_flag(match==1))
        S(j).flag = qual_flag(match==1);
    else
        S(j).flag = 1;
        S(j).Author = 'KC';
    end
end
%% Adding Decdate to the Nan values and arrange in ascending order
for i = 1:length(S)
    if isnan(S(i).DecDate)
        S(i).DecDate = decyear(datetime(S(i).Date,'InputFormat','yyyy-MM-dd'));
    end
end
[~,index] = sortrows([S.DecDate].'); 
S = S(index); 
clear index

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
% x = 0;
% for i = 1:length(indss)
%     del = indss(i);
%     S_filt(del-x) = [];
%     x=x+1;
% end


%% Plotting all the terminus
for i=1:length(S_filt)
    plot(S_filt(i).X,S_filt(i).Y)
    hold on;
end

%% 3. Initialize


%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Terminus Positions
AllTermini = S_filt;

% Satellite Image
DEM_path = fullfile(root_path,glacier_name,glacier_name);
% DEMs and Bedrock Topography
cd(root_path)
old_dir = dir([DEM_path '/*_19*' '.tif']);
OlderDEM = [old_dir.folder,'/',old_dir.name];
new_dir = dir([DEM_path '/*_20*' '.tif']);
ArcticDEM = [new_dir.folder,'/',new_dir.name];
BedDEM = [root_path,'BedTopo/BedmachineV4_BedTopography.tif'];
Landsat_Im = [root_path,glacier_name,'/',glacier_name,'.tif'];
thickness_path = [root_path, 'Thickness_Geotiffs','/'];
% glacier_name_2 = 'Sermeq_Kujalleq(NW)'; only use this if there's a conflict in name and you wanna re-initialize
%if ~exist('glacier_name_2', 'var')
%    thickness_aero_path = dir([thickness_path, glacier_name, '*' '.tif']);
%else   
%    thickness_aero_path = dir([thickness_path, glacier_name_2, '*' '.tif']);
%end

thickness_aero = [thickness_path, thickness_aero_path(3).name];
thickness_arctic = [thickness_path, thickness_aero_path(4).name];
err_aero = [thickness_path, thickness_aero_path(1).name];
err_arctic = [thickness_path, thickness_aero_path(2).name];

% [L8Im,RL8] = readgeoraster(Landsat_Im); % Input Landsat Image

% Get Glacier Name
GlacierName = glacier_name;
outFolderData = [root_path,'Results/Data/'];
output_path = [root_path,'Results/Output/'];
outFolderImages = [root_path,'Results/Images/'];


% Create Output Folders
if ~exist(outFolderData, 'dir')
    if ~exist(outFolderData, 'dir')
        mkdir(outFolderData)
        mkdir(outFolderImages)
    end
end

% Specify Colors and loop length
Ncolor= [99 172 190]/255;
N= length(AllTermini);

%% 4. Draw Fjord Walls (Upper and Lower)
if exist(strcat(outFolderData,'FjordBoundaries/','UpperFW_',GlacierName,'.csv'))
    disp('Boundary files found - Loading files')
    FW = readmatrix(strcat(outFolderData,'FjordBoundaries/','UpperFW_',GlacierName,'.csv'));
    FWl = readmatrix(strcat(outFolderData,'FjordBoundaries/','LowerFW_',GlacierName,'.csv'));
    [UpperPoly,LowerPoly] = CreatePolygon(FW,FWl);
else
    [FW,FWl,UpperPoly,LowerPoly] = drawFjordWalls(L8Im,RL8,GlacierName,outFolderData);
end

%% 5. Get closest point on Fjordwall to end/startpoint of Terminus Position
% Checks which way the terminus is drawn (e.g. left to right/right to left 
% or up-down/down/up)and flips them to a consistent start location if
% necessary
f=waitbar(0,'Processing');
for i=1:N
%     tePos{i} = shaperead(fullfile(SplitTerminiFolder,TP_Names{i}));
    tePos{i} = AllTermini(i);
    if isnan(tePos{i}.X(end))
        tePos{i}.X = tePos{i}.X(:,any(~isnan(tePos{i}.X),1));
        tePos{i}.Y = tePos{i}.Y(:,any(~isnan(tePos{i}.Y),1));
    else
        tePos{i}.X = tePos{i}.X;
        tePos{i}.Y  =tePos{i}.Y;     
    end
    if ~isnan(tePos{i}.X(end))
        TerminusEndPoint = tePos{i}.X(end);
        TerminusEndPoint(:,2) =tePos{i}.Y(end);
    else
        TerminusEndPoint = tePos{i}.X(end-1);
        TerminusEndPoint(:,2) =tePos{i}.Y(end-1);
    end
    TerminusStartPoint = tePos{i}.X(1,1);
    TerminusStartPoint(:,2) = tePos{i}.Y(1,1);
    %Check distance of Terminus EndPoints to either fjord wall to see if
    % they are consistently draw N/S, S/N, or E/W, W/E
    [~,DEndToUpper,~] = distance2curve(FW,TerminusEndPoint);
    [~,DEndToLower,~] = distance2curve(FWl,TerminusEndPoint);

    if min(DEndToLower)<min(DEndToUpper) 
        tePos{i}.X = fliplr([tePos{i}.X]);
        tePos{i}.Y = fliplr([tePos{i}.Y]);
    else
        tePos{i}.X =tePos{i}.X ;
        tePos{i}.Y =tePos{i}.Y;
    end

    NewTerminusEndPoint = tePos{i}.X(end);
    NewTerminusEndPoint(:,2) =tePos{i}.Y(end);
    NewTerminusStartPoint = tePos{i}.X(1);
    NewTerminusStartPoint(:,2) =tePos{i}.Y(1);

waitbar(i/N, f, sprintf('Processing drawing direction: %d%%',floor(i/N*100)));
end
close(f)


%% 6. Check where Margin Endposition is located in relation to upper fjord wall boundary 
% Extrapolates terminus position endpoint to upper fjord wall boundary if it 
% doesn't intersect with the fjord wall or crops it to the point of
% intersection
disp('Adjusting terminus end position relative to upper fjord wall...');
f= waitbar(0,'Processing');

for i=1:N
    % Check if terminus position intersects fjord walls more than once
    [xi,yi,] = polyxpoly(tePos{i}.X,tePos{i}.Y,FW(:,1),FW(:,2));
    if size(xi)>1
        EndCoords(:,1) = xi(1);
        EndCoords(:,2) = yi(1);
    else
        EndCoords(:,1) = tePos{i}.X(end);
        EndCoords(:,2) = tePos{i}.Y(end);
    end

    [InUpper,OnUpper,Number] = InPolygonCheck(tePos{i},UpperPoly,0,1);
PointsAbove = Number;
% Crop to intersection with fjord wall point if terminus endpoint is
% located above fjord wall
if InUpper==1 
        disp('Above Upper Fjord Wall');
        PointsAbove = Number;
        if PointsAbove>0
            [InterX,InterY]= polyxpoly(tePos{i}.X,tePos{i}.Y,FW(:,1),FW(:,2));
            if isempty(InterY)
                [~,InterIdx] = pdist2(FW, EndCoords,'euclidean','Smallest',1);
                NewFWEndPointUp = [FW(InterIdx,1),FW(InterIdx,2)]; 
                tePos{i}.X(end) = NewFWEndPointUp(:,1);
                tePos{i}.Y(end) = NewFWEndPointUp(:,2);
            else
                StorageEndPoint = [InterX(1),InterY(1)];
                [~,IntIdxUp] = pdist2(FW, StorageEndPoint,'euclidean','Smallest',1);
                NewFWEndPointUp= [FW(IntIdxUp,1),FW(IntIdxUp,2)];
                tePos{i}.X = tePos{i}.X(1:end-PointsAbove);
                tePos{i}.Y = tePos{i}.Y(1:end-PointsAbove);
                tePos{i}.X(end) = NewFWEndPointUp(:,1);
                tePos{i}.Y(end) = NewFWEndPointUp(:,2);
            end
        end
        if InUpper==1 && OnUpper==1 && PointsAbove>0 %if the terminus intersects fjord walls twice
            [InterXDouble,InterYDouble]= polyxpoly(tePos{i}.X,tePos{i}.Y,FW(:,1),FW(:,2));
            StorageEndPoint = [InterXDouble(1),InterYDouble(1)];
            [~,IntIdxUp] = pdist2(FW, StorageEndPoint,'euclidean','Smallest',1);
            NewFWEndPointUp= [FW(IntIdxUp,1),FW(IntIdxUp,2)];
            tePos{i}.X = tePos{i}.X(1:end-PointsAbove);
            tePos{i}.Y = tePos{i}.Y(1:end-PointsAbove);
            tePos{i}.X(end) = NewFWEndPointUp(:,1);
            tePos{i}.Y(end) = NewFWEndPointUp(:,2);
        end
% Extrapolate to intersection with fjord walls if terminus endpoint is 
% located below fjord wall
    elseif InUpper==0 && PointsAbove == 0
        disp('Below Upper Fjord Wall');
        [~,InterIdx] = pdist2(FW, EndCoords,'euclidean','Smallest',1);
        NewFWEndPointUp = [FW(InterIdx,1),FW(InterIdx,2)];
        tePos{i}.X(end) = NewFWEndPointUp(:,1);
        tePos{i}.Y(end) = NewFWEndPointUp(:,2);
    elseif OnUpper==1
        % Do nothing if terminus endpoint is already on fjord wall
        disp('On Upper Fjord Wall');
        tePos{i}.X(end) = tePos{i}.X(end);
        tePos{i}.Y(end) = tePos{i}.Y(end);
end
    waitbar(i/N, f, sprintf('Processing: %d%%',floor(i/N*100)));
end
close(f)



%%
f= waitbar(0,'Processing');
for i=1:N
    % Check if terminus position intersects fjord walls more than once
    [xi,yi,] = polyxpoly(tePos{i}.X,tePos{i}.Y,FW(:,1),FW(:,2));
    if size(xi)>1
        EndCoords(:,1) = xi(1);
        EndCoords(:,2) = yi(1);
    else
        EndCoords(:,1) = tePos{i}.X(end);
        EndCoords(:,2) = tePos{i}.Y(end);
    end

    [InUpper,OnUpper,Number] = InPolygonCheck(tePos{i},UpperPoly,0,1);
% Crop to intersection with fjord wall point if terminus endpoint is
% located above fjord wall
    if InUpper==1 
        disp('Above Upper Fjord Wall');
        PointsAbove = Number;
%         if PointsAbove>0
%             [InterX,InterY]= polyxpoly(tePos{i}.X,tePos{i}.Y,FW(:,1),FW(:,2));
%             StorageEndPoint = [InterX(1),InterY(1)];
%             [~,IntIdxUp] = pdist2(FW, StorageEndPoint,'euclidean','Smallest',1);
%             NewFWEndPointUp= [FW(IntIdxUp,1),FW(IntIdxUp,2)];
%             tePos{i}.X = tePos{i}.X(1:end-PointsAbove);
%             tePos{i}.Y = tePos{i}.Y(1:end-PointsAbove);
%             tePos{i}.X(end) = NewFWEndPointUp(:,1);
%             tePos{i}.Y(end) = NewFWEndPointUp(:,2);
%         end
        if InUpper==1 && OnUpper==1 && PointsAbove>0
            [InterXDouble,InterYDouble]= polyxpoly(tePos{i}.X,tePos{i}.Y,FW(:,1),FW(:,2));
            StorageEndPoint = [InterXDouble(1),InterYDouble(1)];
            [~,IntIdxUp] = pdist2(FW, StorageEndPoint,'euclidean','Smallest',1);
            NewFWEndPointUp= [FW(IntIdxUp,1),FW(IntIdxUp,2)];

            tePos{i}.X = tePos{i}.X(1:end-PointsAbove);
            tePos{i}.Y = tePos{i}.Y(1:end-PointsAbove);
            tePos{i}.X(end) = NewFWEndPointUp(:,1);
            tePos{i}.Y(end) = NewFWEndPointUp(:,2);
        end
% Extrapolate to intersection with fjord walls if terminus endpoint is 
% located below fjord wall
    elseif InUpper==0 
        disp('Below Upper Fjord Wall');
        [~,InterIdx] = pdist2(FW, EndCoords,'euclidean','Smallest',1);
        NewFWEndPointUp = [FW(InterIdx,1),FW(InterIdx,2)];
        tePos{i}.X(end) = NewFWEndPointUp(:,1);
        tePos{i}.Y(end) = NewFWEndPointUp(:,2);
    elseif OnUpper==1
% Do nothing if terminus endpoint is already on fjord wall
        disp('On Upper Fjord Wall');
        tePos{i}.X(end) = tePos{i}.X(end);
        tePos{i}.Y(end) = tePos{i}.Y(end);
    end
   
clear  StorageEndPoint InterX InterY InterXDouble InterYDouble InterIdx InterSdist 
clear xi yi NewFWEndPointUp InUpper OnUpper PointsAbove EndCoords IntIdxUp InterS Int 
clear TempEndP ClosestIntIdx EndP

    % Repeat above steps for lower fjord wall and terminus start points
    StartCoords(:,1) = tePos{i}.X(1);
    StartCoords(:,2) = tePos{i}.Y(1);
    [xiD,~,LSeg] = polyxpoly(tePos{i}.X,tePos{i}.Y,FWl(:,1),FWl(:,2));
    if length(xiD)>1
        tePos{i}.X = tePos{i}.X(max(LSeg)-1:end);
        tePos{i}.Y = tePos{i}.Y(max(LSeg)-1:end);
    end
    [InLower,OnLower,PointsBelow] = InPolygonCheck(tePos{i},LowerPoly,1,0);
    if InLower==0 
        [~,InterIdx] = pdist2(FWl, StartCoords,'euclidean','Smallest',1);
        NewXYDown = zeros(length(tePos{i}.X)+1,2);
        NewXYDown(1,1) = FWl(InterIdx,1);
        NewXYDown(1,2) = FWl(InterIdx,2);
        NewXYDown(2:end,1) = tePos{i}.X;
        NewXYDown(2:end,2) = tePos{i}.Y;
        tePos{i}.X = NewXYDown(:,1);
        tePos{i}.Y = NewXYDown(:,2);
    elseif InLower==1 
        % Check if terminus position intersects fjord walls more than once and if so, crop to first intersection point
        if PointsBelow>0
           [InterXLower,InterYLower]= polyxpoly(tePos{i}.X,tePos{i}.Y,FWl(:,1),FWl(:,2)); 
            tePos{i}.X = tePos{i}.X(1+PointsBelow:end);
            tePos{i}.Y = tePos{i}.Y(1+PointsBelow:end);
            NewXYDown = zeros(length(tePos{i}.X)+1,2);
            NewXYDown(1,1) = InterXLower(1);
            NewXYDown(1,2) = InterYLower(1);
            NewXYDown(2:end,1) = tePos{i}.X;
            NewXYDown(2:end,2) = tePos{i}.Y;
            tePos{i}.X = NewXYDown(:,1); 
            tePos{i}.Y = NewXYDown(:,2);
        end
    elseif InLower==1 && OnLower==1 && PointsBelow>0
            [InterXLower,InterYLower]= polyxpoly(tePos{i}.X,tePos{i}.Y,FWl(:,1),FWl(:,2)); 
            tePos{i}.X = tePos{i}.X(1+PointsBelow:end);
            tePos{i}.Y = tePos{i}.Y(1+PointsBelow:end);
            NewXYDown = zeros(length(tePos{i}.X)+1,2);
            NewXYDown(1,1) = InterXLower(1);
            NewXYDown(1,2) = InterYLower(1);
            NewXYDown(2:end,1) = tePos{i}.X;
            NewXYDown(2:end,2) = tePos{i}.Y;
            tePos{i}.X = NewXYDown(:,1);
            tePos{i}.Y = NewXYDown(:,2);
    elseif OnLower==1
        tePos{i}.X(1) = tePos{i}.X(1);
        tePos{i}.Y(1) = tePos{i}.Y(1);
    end
       waitbar(i/N, f, sprintf('Processing terminus positions: %d%%',floor(i/N*100)));
end
close(f)
clear NewXYDown InterXLower InterYLower InterIdx  
clear OnLower InLower PointsBelow StartCoords xiD LSeg

%% 7. Check Extrapolation
% % Check if length of extrapolation exceeds actual terminus delineation
for i=1:length(tePos)
    ExCoordsStart{i} =[tePos{i}.X(1:2)';tePos{i}.Y(1:2)'];
    ExCoordsEnd{i} = [tePos{i}.X(end-1:end)';tePos{i}.Y(end-1:end)'];
    DelinCoords{i} = [tePos{i}.X';tePos{i}.Y'];
    dDelin{i} = diff(DelinCoords{i}); 
    
    DelinLen(i,1) = sum(sqrt(sum(dDelin{i}.*dDelin{i},2)));
    ExTotalLen(i,1)= norm(ExCoordsStart{i}(1:2,1)-ExCoordsStart{i}(1:2,2)) + norm(ExCoordsEnd{i}(1:2,1)-ExCoordsEnd{i}(1:2,2));
    if (ExTotalLen(i,1)<=DelinLen(i,1))==1 
        tePosTemp{i}=tePos{i};
    else
        tePosTemp{i}=[];
        disp('Extrapolation > Delineation')
    end
end
tePos=tePosTemp(~cellfun('isempty',tePosTemp));
N=length(tePos);
clear tePosTemp dTerminus TerminusLength dExtrapolation ExtrapolationLength

%% 8. Define Boundary upstream of last Terminus Position
% as the most retreated position of the terminus might not be the most
% retreated across the terminus, an arbitrary inland boundary should be
% drawn upstream of the most retreated terminus position. The discrepancy
% will be subtracted later on. 
if exist(strcat(outFolderData,'UpstreamBoundary/',GlacierName,'.csv'))
    disp('Boundary files found - Loading files')
    BLN = readmatrix(strcat(outFolderData,'UpstreamBoundary/',GlacierName,'.csv'));
else
    BLN =LowerBoxBoundary(L8Im,RL8, FW,FWl,tePos);
end
%% 9. Getting the center velocity

S_fil = tePos;
for k = 1:length(tePos)
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

%calculate terminus velocity from distance change over time P
%then filter according to if dTermdt > 2*speed

%% 10. Filtering using ITS_LIVE Velocity
%% Run only once
%rate = glacier_name.velocity;
rate_thresh = 2*vel;
Shp_filtered = S_fil;
ind= []; counter = 1;
for k = 1:length(centerdist) - 2
    rate = (centerdist(k + 1) - centerdist(k)) / (termdate(k + 1) - termdate(k));
    
    if rate > rate_thresh
        ind(counter) = k + 1;
        counter = counter + 1;
    end
    
    if rate < 0
        date_next = termdate(k + 2);
        dist_next = centerdist(k + 2);
        rate_next = dist_next / date_next;
        
        if rate_next > rate_thresh
            ind(counter) = k + 1;
            counter = counter + 1;
        end
    end
end

Shp_filtered(:,ind) = [];  % removing the cells with dips and spikes
for i=1:length(Shp_filtered)
    plot(Shp_filtered{1,i}.X,Shp_filtered{1,i}.Y)
    hold on;
end

%% Doing with the midpoints velocity; also only run once

S_data = Shp_filtered;
ind= []; counter = 1;
first_x = (BLN(1,1)+BLN(2,1))/2; first_y = (BLN(1,2)+BLN(2,2))/2;
for i = 1:length(S_data)
    x1 = Shp_filtered{1, i}.X(1); x2 = Shp_filtered{1, i}.X(end);
    y1 = Shp_filtered{1, i}.Y(1); y2 = Shp_filtered{1, i}.Y(end);
    date_shp(i) = datenum(Shp_filtered{1, i}.Date);
    x_mid(i) = (x1+x2)/2; y_mid(i) = (y1+y2)/2;
end
for k = 1:length(S_data)-2
    dist(k) = sqrt((x_mid(k) - first_x)^2 + (y_mid(k)-first_y)^2);
end
for k = 1:length(S_data)-4
    rate_vel(k,1) = (dist(k+1)-dist(k))/(date_shp(k+1)-date_shp(k));
    rate_vel(k,2) = date_shp(k+1);
end
for k = 1:length(S_data)-4
    if rate_vel(k)>rate_thresh*10
            ind(counter) = k+1;
            counter = counter+1;
    elseif rate_vel(k)<0
        date_diff = date_shp(k+2)-date_shp(k+1); dist_next = dist(k+2)-dist(k+1);
        rate_next = dist_next/date_diff;
        if rate_next > rate_thresh*10
            ind(counter) = k+1;
            counter = counter+1;
        end
    end
end
S_test = Shp_filtered;
S_test(:,ind) = [];
Shp_filtered(:,ind) = [];
%% Run this if you wanna plot
for i=1:length(S_test)
    plot(S_test{1,i}.X,S_test{1,i}.Y);
    hold on;
end

%%11. Asking if the user wants to change the upstream boundary; have to have
% the landsat image as well; make sure to load that in step 3
plot([BLN(1,1),BLN(2,1)],[BLN(1,2),BLN(2,2)],'r','linewidth',1.5)
hold on;
for i = 1:length(Shp_filtered)
    plot(Shp_filtered{1,i}.X,Shp_filtered{1,i}.Y)
    hold on;
end
prompt = "Press 1 if you want to redraw the upstream boundary ";
inp = input(prompt);
while inp == 1
    x_trans = input('Enter the value for translation in x direction');
    y_trans = input(['Enter the value for translation in y direction']);
    BLN(1,1) = BLN(1,1) + x_trans;
    BLN(1,2) = BLN(1,2) + y_trans;
    BLN(2,1) = BLN(2,1) + x_trans;
    BLN(2,2) = BLN(2,2) + y_trans;
    plot([BLN(1,1),BLN(2,1)],[BLN(1,2),BLN(2,2)],'r','linewidth',1.5)
    hold on;
    for i = 1:length(Shp_filtered)
        plot(Shp_filtered{1,i}.X,Shp_filtered{1,i}.Y)
        hold on;
    end
    prompt = "Press 1 if you want to redraw the upstream boundary ";
    inp = input(prompt);
end

close;
%% 12. Crop the fjord wall boundaries to the specific terminus boundaries
% Needed to create a polygon in the next step. 

for i=1:N
    FWUpperStart = [tePos{i}.X(end),tePos{i}.Y(end)];
    FWUpperEnd =[BLN(1,1),BLN(1,2)];
    FWLowerStart=[tePos{i}.X(1),tePos{i}.Y(1)];
    FWLowerEnd= [BLN(end,1),BLN(end,2)];

    UpperStartIdx = find(FW==round(FWUpperStart,25));
    if isempty(UpperStartIdx)==1
        UpperStartIdx = dsearchn(FW, FWUpperStart);
    end
    [EndXUp,EndYUp] = polyxpoly(FW(:,1),FW(:,2),FWUpperEnd(:,1),FWUpperEnd(:,2));
    StorageEndPointUp = [EndXUp,EndYUp];
    [~,UpperEndIdx] = pdist2(FW, StorageEndPointUp,'euclidean','Smallest',1);

    IndFW{i}.X = FW(UpperStartIdx:UpperEndIdx,1);
    IndFW{i}.Y = FW(UpperStartIdx:UpperEndIdx,2);

    %%%%%% Same for lower boundary %%%%%
    LowerStartIdx = find(FWl==FWLowerStart);
    if isempty(LowerStartIdx)==1 
        LowerStartIdx = dsearchn(FWl, FWLowerStart);
        StartIdxLoKeep=1;
    end
    [EndXLo,EndYLo] = polyxpoly(FWl(:,1),FWl(:,2),FWLowerEnd(:,1),FWLowerEnd(:,2));

    StorageStartPointLo = [EndXLo,EndYLo];
    [~,LowerEndIdx] = pdist2(FWl, StorageStartPointLo,'euclidean','Smallest',1);

    IndFWl{i}.X = FWl(LowerStartIdx:LowerEndIdx,1);
    IndFWl{i}.Y = FWl(LowerStartIdx:LowerEndIdx,2);
    IndFWl{i}.X(1)=FWLowerStart(:,1);
    IndFWl{i}.Y(1)=FWLowerStart(:,2);

end
S_temp = Shp_filtered;
disp('Fjord wall boundaries cropped');

%% 13.Clearing the length that are longer than 1.5 or 0.5 times the mean of all the termini
% ONLY RUN ONCE!
for i=1:length(Shp_filtered)
    x_coord = rmmissing(Shp_filtered{1,i}.X); %every last elemnt is a NaN
    y_coord = rmmissing(Shp_filtered{1,i}.Y);
    d = diff([x_coord(:) y_coord(:)]);
    total_length = sum(sqrt(sum(d.*d,2)));
    px(i) = total_length/1000;
end
% Try to do with median fitering
% other idea would be with Orientation of terminus; tangent to every point 
ml= median(px); ctr=1;
for i=1:length(px)
    if px(i) > 1.5 * ml || px(i) < 0.5 * ml
        iidx(ctr) = i;
        ctr = ctr+1;
    end
end
if ~exist("iidx")
    disp('No points found: so does not require step 12(b)')
else
    disp('Run step 12(b)')
end
% Deleting the flagged termini
iidx = iidx';
Shp_filtered(:,iidx) = [];
disp('Flagged termini deleted')

%% if you wanna plot, run this section.
for i=1:length(Shp_filtered)
    plot(Shp_filtered{1,i}.X,Shp_filtered{1,i}.Y)
    hold on;
end
hold on
plot(FW(:,1),FW(:,2));
hold on

plot(FWl(:,1),FWl(:,2));%%
% for i=1:length(Shp_filtered)
%     plot(Shp_filtered{1,i}.X,Shp_filtered{1,i}.Y)
%     hold on;
% end
% hold on
% plot(FW(:,1),FW(:,2));
% plot(FWl(:,1),FWl(:,2));

%% 14. Saving cells into a structure
fieldNames = fieldnames(Shp_filtered{1});
fn = {};

for i = 1:length(fieldNames)
    fn{i} = fieldNames{i};
end
for i = 1:length(Shp_filtered)
    combinedStruct(i) = Shp_filtered{i};
end
%% 15. To save the filtered termini
% Can only run once
for i=1:length(combinedStruct)
    if isnan(combinedStruct(i).Center_X)
        combinedStruct(i).Center_X = combinedStruct(i).pointLon;
        combinedStruct(i).Center_Y = combinedStruct(i).pointLat;
    end
end
combinedStruct = rmfield(combinedStruct, 'Time');
for i=1:length(combinedStruct)
    if isempty(combinedStruct(i).flag)
        combinedStruct(i).flag = 1;
    end
end

output_file = fullfile(root_path, glacier_name, 'Terminus_Positions', append('Filtered_Termini_', glacier_name, '.shp'));
shapewrite(combinedStruct, output_file);
%% Not plotting greater than 2023-03-31; only run once
latest_date = datenum('2023-03-31');
count = 1;
for i =1:length(Shp_filtered)
    if  datenum(Shp_filtered{1,i}.Date)>latest_date
        inds(count) = i;
        count = count+1;
    end
end
if ~exist("iidx2")
    disp('No points found: so does not require step below')
else
    disp('Run step below')
end
iidx2 = iidx2';
Shp_filtered(:,iidx2) = [];
%% 14. Concatenate Boundaries, create polygon and calculate area
% After Filtering terminus velocity
%create polygon for each terminus delineation

tePos = Shp_filtered;
for i=1:length(tePos)
    if size(tePos{i}.X,1) > size(tePos{i}.X,2)
        TerminusPosition = tePos{i}.X; TerminusPosition(:,2)=tePos{i}.Y;
    else
        TerminusPosition = tePos{i}.X'; TerminusPosition(:,2)=tePos{i}.Y';
    end
        UpperBoundary=IndFW{i}.X; UpperBoundary(:,2)=IndFW{i}.Y;
        LowerBoundary=IndFWl{i}.X; LowerBoundary(:,2)=IndFWl{i}.Y;
        Poly{i}= polyshape([TerminusPosition; UpperBoundary;BLN;flipud(LowerBoundary)]);
        VolumeTS(i,1) = datenum(tePos{i}.Date);
        VolumeTS(i,2) = area(Poly{i})/1e6;
end

%identify polygon bounding box
for j=1:length(Poly)
    minVals(j,1) = min(Poly{j}.Vertices(:,1));
    minVals(j,2) = min(Poly{j}.Vertices(:,2));
    maxVals(j,1) = max(Poly{j}.Vertices(:,1));
    maxVals(j,2) = max(Poly{j}.Vertices(:,2));
end
%extract data and add to a structure
for i = 1:length(tePos)
    AreaTS(i,1) = datenum(tePos{i}.Date);
    AreaTS(i,2) = area(Poly{i})/1e6;
end
%plot a terminus area time series
if length(tePos)>100 %plot only a subset since there are a ton
    plot_inc = round(length(tePos)/100);
    for i=1:100
        figure(5)
        pause(0.01);
        set(gcf,'WindowState','maximized')
        [n,~]=numSubplots(100);
        subplot(n(1),n(2),i)
        N=100;
        cmap=parula(N);
        plot(Poly{i*plot_inc+1},'FaceColor',cmap(i,:),'EdgeColor','k');
        set(gca,'XTick',[], 'YTick', [])
        hold on
        xlim([min(minVals(:,1)) max(maxVals(:,1))]);
        ylim([min(minVals(:,2)) max(maxVals(:,2))]);
        title(tePos{i*plot_inc+1}.Date)
        sgtitle('Determined areas for available terminus positions');
    end
    saveas(gcf,strcat(outFolderImages,'Area_',GlacierName,'_2.jpg'));
    close()
else
    for i=1:length(tePos)
%         TerminusPosition = tePos{i}.X'; TerminusPosition(:,2)=tePos{i}.Y';
%         UpperBoundary=IndFW{i}.X; UpperBoundary(:,2)=IndFW{i}.Y;
%         LowerBoundary=IndFWl{i}.X; LowerBoundary(:,2)=IndFWl{i}.Y;
%         Poly{i}= polyshape([TerminusPosition; UpperBoundary;BLN;flipud(LowerBoundary)]);

        figure(5)
        pause(0.00001);
        set(gcf,'WindowState','maximized')
        [n,~]=numSubplots(length(tePos));
        subplot(n(1),n(2),i)
        N=length(tePos);
        cmap=parula(N);
        plot(Poly{i},'FaceColor',cmap(i,:),'EdgeColor','k');
        set(gca,'XTick',[], 'YTick', [])
        hold on
        xlim([min(minVals(:,1)) max(maxVals(:,1))]);
        ylim([min(minVals(:,2)) max(maxVals(:,2))]);
%         title(tePos{i}.Date)
        sgtitle('Determined areas for available terminus positions');
    end
    saveas(gcf,strcat(outFolderImages,'Area_',GlacierName,'.jpg'));
    close()
end

%save datasets
save([outFolderData,GlacierName,'_terminus-polygons.mat'],'tePos','-v7.3');
save([outFolderData,GlacierName,'_area-timeseries.mat'],'AreaTS','-v7.3');

%% 15. Calculate surface elevation change based on old and new DEM 

%create a new matrix for frontal ablation calculations
VolumeTS = AreaTS;
% Calculate the volume using thickness from BedMachine thickness estimates
k = strfind(thickness_aero,'19');
EarliestDEMYear = str2num(thickness_aero(k:k+3));
l = strfind(thickness_arctic,'20');
LatestDEMYear = str2num(thickness_arctic(l:l+3));
Tdates = datetime(single(VolumeTS(:,1)),'ConvertFrom','datenum');
[DEM_list] = ThicknessEstimates(thickness_aero,thickness_arctic,err_aero,err_arctic,VolumeTS(:,1:2),Poly);
for i = 1:length(Poly)
    clear IceThickness
    ROImask = inpolygon(DEM_list{i}.X,DEM_list{i}.Y,Poly{i}.Vertices(:,1),Poly{i}.Vertices(:,2));
    IceThickness.X = DEM_list{i}.X(ROImask);
    IceThickness.Y = DEM_list{i}.Y(ROImask);
    IceThickness.H = DEM_list{i}.H(ROImask);
    IceThickness.L1 = DEM_list{i}.L1(ROImask);IceThickness.L2 = DEM_list{i}.L2(ROImask);
    IceThickness.U1 = DEM_list{i}.U1(ROImask);IceThickness.U2 = DEM_list{i}.U2(ROImask);
    thickness = (nanmean(IceThickness.H,'all'))/1000;%use the MEAN, not the median          
    if year(Tdates(i)) > EarliestDEMYear && year(Tdates(i)) < LatestDEMYear 
        if nanmean(IceThickness.U1,'all')>nanmean(IceThickness.U2,'all')
            upp_lim = abs(IceThickness.U1 -  IceThickness.H)/1000;% subtracting from the thickness to get the range
            upper_thick_err(i) = sqrt (nansum (upp_lim.^2))/length(upp_lim);
        else
            upp_lim = abs(IceThickness.U2 -  IceThickness.H)/1000;
            upper_thick_err(i) = sqrt (nansum (upp_lim.^2))/length(upp_lim);
        end
        if nanmean(IceThickness.L1,'all')<nanmean(IceThickness.L2,'all')
            low_lim = abs(IceThickness.H -  IceThickness.L1)/1000;
            lower_thick_err(i) = sqrt (nansum (low_lim.^2))/length(low_lim);
        else
            low_lim = abs(IceThickness.H -  IceThickness.L2)/1000;
            lower_thick_err(i) = sqrt (nansum (low_lim.^2))/length(low_lim);
         end
    else 
        upper_thick_err(i) = sqrt (nansum ((IceThickness.U1/1000).^2))/length(IceThickness.U1);
        lower_thick_err(i) = sqrt (nansum ((IceThickness.L1/1000).^2))/length(IceThickness.U1);
    end
    VolumeTS(i,3) = thickness;
    VolumeTS(i,4) = VolumeTS(i,2)*VolumeTS(i,3);%volume in km^3
    clear upp_lim low_lim;
end
err_bounds(:,1) = upper_thick_err';
err_bounds(:,2) = lower_thick_err';
thick_err = err_bounds;
dates_num = datetime(VolumeTS(:,1),'ConvertFrom','datenum');
dates = [dates_num;flipud(dates_num)];
err_upp = [VolumeTS(:,3)+err_bounds(:,1) ; flipud(VolumeTS(:,3))-flipud(err_bounds(:,2))];
fill(dates(1:end),err_upp(1:end),[.7 .7 .7])
hold on;
plot(dates_num(1:end),VolumeTS(1:end,3),'red')

%%16. Saving the dataset
Frontal_Mass = VolumeTS;
Frontal_Mass(:,5) = Frontal_Mass(:,4)*0.917;
date = datetime(Frontal_Mass(:,1),'convertfrom','datenum');
area = Frontal_Mass(:,2); thickness = Frontal_Mass(:,3); volume = Frontal_Mass(:,4); terminus_mass = Frontal_Mass(:,5); 
term_mass = table(date,area,thickness,volume,terminus_mass);
output_path_mass = [output_path,'Terminus_Mass/'];
writetable(term_mass, [output_path_mass,strcat('Term_Mass_updated_',GlacierName,'_',region,'.csv')], 'Delimiter', ',');
%% if you wanna plot
two_thick = readtable('/Users/amankc/Terminus_Ablation/Results/Output/Terminus_Mass/Term_Mass_Sermeq_Avanarleq_NW.csv');
bedmachine = readtable('/Users/amankc/Terminus_Ablation/Results/Output/Terminus_Mass/Term_Mass_testSermeq_Avanarleq_NW.csv');
plot(two_thick.date,two_thick.thickness)
hold on
plot(bedmachine.date,bedmachine.thickness)
legend('old thickness','bedmachine thickness')

