%%% Calculate frontal ablation and terminus volume change over time.

%% 3. Initialize

%root_path = '/Users/amankc/Downloads/Frontal_Ablation/Data/';
% AllTerminiPath = '/Users/amankc/Downloads/Frontal_Ablation/Data/TP/Helheim_Gletsjer.shp';
% AllTermini = shaperead(AllTerminiPath);

%%%%%%%%%%%%%%%%%%%%%% Load Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Terminus Positions
%AllTerminiPath = [root_path,'TermPicks/Kangerdlussuaq152.shp'];
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
% [L8Im,RL8] = readgeoraster(Landsat_Im); % Input Landsat Image
[L8Im,RL8] = readgeoraster(ArcticDEM); % Input Landsat Image
% dhdt = readmatrix('Data/SurfaceChange/dHdt_Mankoff.csv');
% Split Terminus Positions
% TerminusPositionFolder = [root_path,'TP/'];
% TP_files = dir(fullfile(TerminusPositionFolder,'*.shp'));
% TP_Names = {TP_files.name};

% % Velocities
% VelFolder = 'Velocity_Helheim/';
% Vel_files = dir(fullfile(VelFolder,'*.tif'));
% Vel_Names = {Vel_files.name};
% % Create list of available velocity dates
% for i=1:length(Vel_Names)
%     DateList(i,1)=datetime(Vel_Names{i}(1:end-4),'InputFormat','yyyy-MM');
%     DateNumList(i,1) = single(datenum(DateList(i,1)));
% end

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
    tePos{i} = AllTermini(i);
    if isnan(tePos{i}.X(end))
        tePos{i}.X = tePos{i}.X(:,any(~isnan(tePos{i}.X),1));
        tePos{i}.Y = tePos{i}.Y(:,any(~isnan(tePos{i}.Y),1));
    else
        tePos{i}.X = tePos{i}.X;
        tePos{i}.Y  =tePos{i}.Y;     
    end
    TerminusEndPoint = tePos{i}.X(end);
    TerminusEndPoint(:,2) =tePos{i}.Y(end);
    TerminusStartPoint = tePos{i}.X(1,1);
    TerminusStartPoint(:,2) = tePos{i}.Y(1,1);
    %Check distance of Terminus EndPoints to either fjord wall to see if
    % they are consistently draw N/S, S/N, or E/W, W/E
    [~,DEndToUpper,~] = distance2curve(FW,TerminusEndPoint);
    [~,DEndToLower,~] = distance2curve(FWl,TerminusEndPoint);
    [~,DStartToUpper,~] = distance2curve(FW,TerminusStartPoint);
    [~,DStartToLower,~] = distance2curve(FWl,TerminusStartPoint);

    if min(DEndToLower)<min(DEndToUpper) 
        tePos{i}.X = fliplr([tePos{i}.X]);
        tePos{i}.Y = fliplr([tePos{i}.Y]);
%         disp('Flipping the Table')
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
close()

%% 7. Check where Margin Endposition is located in relation to lower fjord wall boundary 
% Repeats the steps above for the lower boundary
disp('Adjusting terminus end position relative to lower fjord wall...');
f= waitbar(0,'Processing');

for i=1:N
    clear InLower OnLower Number
    % Repeat above steps for lower fjord wall and terminus start points
    StartCoords(:,1) = tePos{i}.X(1);
    StartCoords(:,2) = tePos{i}.Y(1);
    [xyEndToLower,distL,~]=distance2curve(FWl,StartCoords);

    IdxLower = dsearchn(FWl,xyEndToLower);
    if (IdxLower<size(FWl,1))==1 && IdxLower<998
        LineLower(1,:) = FWl(IdxLower(1),:);
        LineLower(2,:) = FWl(IdxLower(1),:);
    else
        LineLower(1,:) = FWl(IdxLower(1)-2,:);
        LineLower(2,:) = FWl(IdxLower(1),:);
    end
    PointXLower = StartCoords(:,1);
    PointYLower = StartCoords(:,2);
    givenXLower = LineLower(:,1);
    givenYLower = LineLower(:,2);

    [InLower,OnLower,Number] = InPolygonCheck(tePos{i},LowerPoly,1,0);
    if InLower==0 
%         disp('Above Lower Fjord Wall');
        [~,InterIdx] = pdist2(FWl, StartCoords,'euclidean','Smallest',1);
        NewStartPoint=[FWl(InterIdx,1), FWl(InterIdx,2)];
        tePos{i}.X(1) = NewStartPoint(:,1);
        tePos{i}.Y(1) = NewStartPoint(:,2);
    elseif InLower==1 
%         disp('Below Lower Fjord Wall');
        PointsBelow = Number;
        if PointsBelow>0
           [InterXLower,InterYLower]= polyxpoly(tePos{i}.X,tePos{i}.Y,FWl(:,1),FWl(:,2)); 
            NewStartPoint = [InterXLower(1),InterYLower(1)];
            tePos{i}.X = tePos{i}.X(1+PointsBelow:end);
            tePos{i}.Y = tePos{i}.Y(1+PointsBelow:end);
            tePos{i}.X(1) = NewStartPoint(:,1);
            tePos{i}.Y(1) = NewStartPoint(:,2);
        end
        if InLower==1 && OnLower==1 && PointsBelow>0
            [InterXLower,InterYLower]= polyxpoly(tePos{i}.X,tePos{i}.Y,FWl(:,1),FWl(:,2)); 
            NewStartPoint = [InterXLower(1),InterYLower(1)];
            tePos{i}.X = tePos{i}.X(1+PointsBelow:end);
            tePos{i}.Y = tePos{i}.Y(1+PointsBelow:end);
            tePos{i}.X(1) = NewStartPoint(:,1);
            tePos{i}.Y(1) = NewStartPoint(:,2);
        end
    elseif OnLower==1
%         disp('On Lower Fjord Wall');
        tePos{i}.X(1) = tePos{i}.X(1);
        tePos{i}.Y(1) = tePos{i}.Y(1);
    end

%plot new terminus positions

%     figure(4)
%     plot(FW(:,1),FW(:,2),'k');
%     hold on
%     plot(FWl(:,1),FWl(:,2),'k');
%     hold on
%     plot(tePos{i}.X(1),tePos{i}.Y(1),'or');
%     title(num2str(i));
    waitbar(i/N, f, sprintf('Processing: %d%%',floor(i/N*100)));
    disp(i)
%     pause(0.1)
end

close()
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
%% Asking if the user wants to change the upstream boundary; have to have
% the landsat image as well; make sure to load that in step 3
plot([BLN(1,1),BLN(2,1)],[BLN(1,2),BLN(2,2)],'r','linewidth',1.5)
hold on;
for i = 1:length(Shp_filtered)
    plot(Shp_filtered{1,i}.X,Shp_filtered{1,i}.Y)
    hold on;
end
prompt = "Press 1 if you want to redraw the upstream boundary ";
inp = input(prompt);
% inp = 1;
% if inp == 1
%     BLN_redrawn =LowerBoxBoundary(L8Im,RL8, FW,FWl,tePos);
% else
%     close
% end
while inp == 1
    x_trans = input(['En' ...
        'ter the value for translation in x direction']);
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
%% 9. Crop the fjord wall boundaries to the specific terminus boundaries
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
disp('Fjord wall boundaries cropped');
%% 12.aClearing the length that are longer than 1.2 or 0.8 times the mean of all the termini

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
%% 12.b Deleting the flagged termini (DO NOT RUN MORE THAN ONCE!)
iidx = iidx';
Shp_filtered(:,iidx) = [];
%%
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
%% 13.a Saving cells into a structure
fieldNames = fieldnames(Shp_filtered{1});
fn = {};

for i = 1:length(fieldNames)
    fn{i} = fieldNames{i};
end
for i = 1:length(Shp_filtered)
    combinedStruct(i) = Shp_filtered{i};
end
%% 13.b To save the filtered termini
% Can only run once
for i=1:length(combinedStruct)
    if isnan(combinedStruct(i).Center_X)
        combinedStruct(i).Center_X = combinedStruct(i).pointLon;
        combinedStruct(i).Center_Y = combinedStruct(i).pointLat;
    end
end
% combinedStruct = rmfield(combinedStruct, 'Time');
for i=1:length(combinedStruct)
    if isempty(combinedStruct(i).flag)
        combinedStruct(i).flag = 1;
    end
end

output_file = fullfile(root_path, glacier_name, 'Terminus_Positions', append('Filtered_Termini_', glacier_name, '.shp'));
shapewrite(combinedStruct, output_file);
%% Not plotting greater than 2023-03-31
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
%%
iidx2 = iidx2';
Shp_filtered(:,iidx2) = [];
%% 14. Concatenate Boundaries, create polygon and calculate area
% After Filtering terminus velocity
%create polygon for each terminus delineation
tePos = Shp_filtered;
for i=1:length(tePos)
        if  length(tePos{i}.X) ~= 1
            x_coord = tePos{i}.X; y_coord = tePos{i}.Y;
            TerminusPosition = x_coord; 
            TerminusPosition(:,2)=y_coord;
        else
            TerminusPosition = tePos{i}.X; 
            TerminusPosition(:,2)=tePos{i}.Y;
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
VolumeTS = AreaTS; %dummy matrix to fill in
%VolumeTS(1,:) = [AreaTS(1,1) round(mean(AreaTS(1:2,1))) AreaTS(1,2)];
%duplicate each observation so that the matrix contains the average date
%with the preceeding observation AND the average with the subsequent
%observation... needed when differencing volumes so that the same elevation
%map is used for consecutive dates!
date = datetime(AreaTS(:,1),'convertfrom','datenum');

%% 16. calculate the volume using the AVERAGE dates (truncate VolumeTS matrix as input to crop off the actual dates)
[DEM_list,SurfaceChange,DEM_error2] = SurfaceElevationChange(ArcticDEM,OlderDEM,BedDEM,VolumeTS(:,1:2),Poly);
for i=1:length(DEM_list)
    %index into Poly is (round(i/2)) since VolumeTS is roughly 2 times the
    %length of Poly and Poly(1) corresponds to VolumeTS(1,:), Poly(2)
    %corresponds to VolumeTS(2,:) AND VolumeTS(3,:), and so on
    %ROImask{i} = inpolygon(DEM_list{i}.X,DEM_list{i}.Y,Poly{round(i/2)}.Vertices(:,1),Poly{round(i/2)}.Vertices(:,2));
    ROImask{i} = inpolygon(DEM_list{i}.X,DEM_list{i}.Y,Poly{i}.Vertices(:,1),Poly{i}.Vertices(:,2));
    IceThickness{i}.X = DEM_list{i}.X(ROImask{i});
    IceThickness{i}.Y = DEM_list{i}.Y(ROImask{i});
    IceThickness{i}.H = DEM_list{i}.H(ROImask{i});
    IceThickness{i}.H2 = DEM_list{i}.H;
    VolumeTS(i,3) = (nanmean(IceThickness{i}.H,'all'))/1000;%use the MEAN, not the median
    VolumeTS(i,4) = VolumeTS(i,2)*VolumeTS(i,3);%volume in km^3
end
dates_num = datetime(VolumeTS(:,1),'ConvertFrom','datenum');
dates = [dates_num;flipud(dates_num)];
DEM_err = DEM_error2./1000;
err_upp = [VolumeTS(:,3)+DEM_err(:,1) ; flipud(VolumeTS(:,3))-flipud(DEM_err(:,2))];
fill(dates(10:end),err_upp(10:end),[.7 .7 .7])
hold on;
plot(dates_num(1:end),VolumeTS(1:end,3),'red')
%%
Frontal_Mass = VolumeTS;
Frontal_Mass(:,5) = Frontal_Mass(:,4)*0.917;
date = datetime(Frontal_Mass(:,1),'convertfrom','datenum');
area = Frontal_Mass(:,2); thickness = Frontal_Mass(:,3); volume = Frontal_Mass(:,4); terminus_mass = Frontal_Mass(:,5); 
term_mass = table(date,area,thickness,volume,terminus_mass);
output_path_mass = [output_path,'Terminus_Mass/'];
writetable(term_mass, [output_path_mass,strcat('Term_Mass_',GlacierName,'_',region,'.csv')], 'Delimiter', ',');
%% save datasets
save([root_path,GlacierName,'_volume-timeseries.mat'],'VolumeTS','-v7.3');
save([root_path,GlacierName,'_FA-timeseries.mat'],'FrontalAB','-v7.3');

