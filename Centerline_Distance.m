clear
close all;
%% 1.
root_path = '/Users/amankc/Terminus_Ablation/';
cd(root_path)
glacier_name = 'Zachariae_Isström';
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
%terminuspos_path = fullfile(root_path,glacier_name,'Merged_Termini',glacier_shp);
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

%% 10. Getting the center velocity
% 
% fieldNames = fieldnames(S_filt{1});
% fn = {};
% for i = 1:length(fieldNames)
%     fn{i} = fieldNames{i};
% end
% for i = 1:length(S_fil)
%     S_fil(i) = S_filt{i};
% end
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

%calculate terminus velocity from distance change over time PROBABLY DON'T
%NEED
%then filter according to if dTermdt > 3*speed using Jukes code

%% 11. Filtering using ITS_LIVE Velocity
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
%% Run only once
Shp_filtered(:,ind) = [];  % removing the cells with dips and spikes

%%
for i=1:length(Shp_filtered)
    plot(Shp_filtered{1,i}.X,Shp_filtered{1,i}.Y)
    hold on;
end
%% Doing with the midpoints velocity

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
%%
S_test = Shp_filtered;
S_test(:,ind) = [];
Shp_filtered(:,ind) = [];
%%
for i=1:length(S_test)
    plot(S_test{1,i}.X,S_test{1,i}.Y);
    hold on;
end