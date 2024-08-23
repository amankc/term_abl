%% 1. This step creates dummy arrays to store the data based on the date of the dataset
clear NearTermMassChange;
NearTerm = readtable([output_path_mass,'Term_Mass_updated_',glacier_name,'_',region,'.csv']);
NearTermMassChange (:,1) = table(Frontal_Mass(:,1));
NearTermMassChange(:,2) = table(Frontal_Mass(:,5));
NearTermMassChange.Properties.VariableNames{1} = 'Dates';
NearTermMassChange.Properties.VariableNames{2} = 'mass';
first_date = datetime(NearTermMassChange.Dates(1),'ConvertFrom','datenum'); 
len = length(NearTermMassChange.Dates);
last_date = datetime(NearTermMassChange.Dates(len),'ConvertFrom','datenum');
duration = calmonths(between(first_date,last_date));
tdt = datetime(year(first_date),month(first_date):duration+month(first_date)+1,1);%generating array of 
% first day of each month
t = datenum(tdt);
mass_rates = zeros(duration,4);%Storing dates, mass changes and number of termini position 
% used to build the average
varnames = {'dates','datenumb','mass'};
TPb = table(datetime(NearTermMassChange.Dates,'ConvertFrom','Datenum'),datenum(NearTermMassChange.Dates),NearTermMassChange.mass,'VariableNames',varnames);
[C,ia] = unique(TPb.datenumb);
TP = TPb(ia,:);
for j = 2:length(t)-2
    terms = find(TP.datenumb>=t(j) & TP.datenumb<t(j+1));
    num_term = length(terms);
    if isempty (terms) %if no data is found between the month, lineraly interpolate 
        % between two earliest and latest dates from the month
        mean_date = round(mean([t(j) t(j+1)])); %getting the middle date of the month
        mass_rates (j-1,1) = mean_date;
        %mass_rates(j,1) = datetime(mean_date,'ConvertFrom','datenum') (To
        %convert into regular dates) Make sure predefined arrays will
        %accept datetime values
        mass_rates (j-1,2) = linear_fun(j,TP,t);
        mass_rates (j-1,3) = 0;
    else
        mass_rates (j-1,1) = round(mean([t(j) t(j+1)]));
        month_dates = t(j+1)-t(j);
        mass_rates(j-1,2) = weighted_fun(terms,TP,month_dates,t);
        mass_rates(j-1,3) = num_term;
    end
end

%% 2. Calculation of Terminus Ablation
Run this, if a glacier is fed through two disharge gates
disch = readtable('/Users/amankc/Terminus_Ablation/gate_D.csv'); %reads the discharge dataset 
date_dis = datetime(table2array(disch(:,1)));
Val_dis = table2array(disch(:,2:end));
disch2 = zeros(size(Val_dis,1)-1,2);
% Manually input Glacier ID
count = 0;
for i = 1:size(Val_dis,2)
    count = count+1;
    if Val_dis(1,i) == glacierID
        disch2(:,1) = datenum(date_dis(2:end,1));
        disch2(:,2) = Val_dis(2:end,count)/365;
    end
end
td = datenum(datetime(year(first_date),month(first_date):duration+month(first_date)+1,16));%generating array of 
% middle day of each month
td = td(2:end-1);
disch2 = unique (disch2,'rows'); %Eliminating duplicated rows
C = td';
years = datetime(mass_rates(:,1),'convertfrom','datenum');
for i = 1:length(t)-2
    terms = find(disch2(:,1)>=t(i) & disch2(:,1)<t(i+1));
    if isempty(terms)
        x1 = find(disch2(:,1)<t(i),1,'last');
        x2 = find(disch2(:,1)>t(i),1,'first');
        if isempty(x1)
            discharge(i) = disch2(x2,2);
        else
            y1 = disch2(x1,2);
            y2 = disch2(x2,2);
            mid_month = (t(i) + t(i+1))/2; t1 = disch2(x1,1); t2 = disch2(x2,1);
            discharge(i) = interp1([t1,t2],[y2,y2],mid_month,'linear');
        end
    else
        sum = 0;
        obs = 0;
        for k = 1:length(terms)
            sum = sum + disch2(terms(k),2);
            obs = obs+1;
        end
        discharge(i) = sum/obs;
    end
end

% discharge = interp1(disch2(:,1),disch2(:,2),C,'linear','extrap');

%% if a glacier is fed through two disharge gates
clear td discharge disch2 disch3 disch4;
glacierID1 = 125; %input these numbers
glacierID2 = 126; %input these numbers
count = 0;
disch = readtable('/Users/amankc/Terminus_Ablation/gate_D.csv');
date_dis = datetime(table2array(disch(:,1)));
Val_dis = table2array(disch(:,2:end));
disch3 = zeros(size(Val_dis,1)-1,2);
disch4 = zeros(size(Val_dis,1)-1,2);
for i = 1:size(Val_dis,2)
    count = count+1;
    if Val_dis(1,i) == glacierID1
        disch3(:,1) = datenum(date_dis(2:end,1));
        disch3(:,2) = Val_dis(2:end,count)/365;
    end
end

count = 0;
for i = 1:size(Val_dis,2)
    count = count+1;
    if Val_dis(1,i) == glacierID2
        disch4(:,1) = datenum(date_dis(2:end,1));
        disch4(:,2) = Val_dis(2:end,count)/365;
    end
end
td = datenum(datetime(year(first_date),month(first_date):duration+month(first_date)+1,16));%generating array of 
% middle day of each month
td = td(2:end-1);
disch2(:,1) = disch3(:,1);
disch2(:,2) = disch3(:,2) + disch4(:,2);
disch2 = unique (disch2,'rows'); %Eliminating duplicated rows
C = td';
years = datetime(mass_rates(:,1),'convertfrom','datenum');
for i = 1:length(t)-1
    terms = find(disch2(:,1)>=t(i) & disch2(:,1)<t(i+1));
    if isempty(terms)
        x1 = find(disch2(:,1)<t(i),1,'last');
        x2 = find(disch2(:,1)>t(i),1,'first');
        if isempty(x1)
            discharge(i) = disch2(x2,2);
        else
            y1 = disch2(x1,2);
            y2 = disch2(x2,2);
            mid_month = (t(i) + t(i+1))/2; t1 = disch2(x1,1); t2 = disch2(x2,1);
            discharge(i) = interp1([t1,t2],[y2,y2],mid_month,'linear');
        end
    else
        sum = 0;
        obs = 0;
        for k = 1:length(terms)
            sum = sum + disch2(terms(k),2);
            obs = obs+1;
        end
        discharge(i) = sum/obs;
    end
end
indx = []; coun = 1;
for i = 1:length(mass_rates)
    mass_rates(i,5) = discharge(i);
end
%     if mass_rates(i,2)>discharge(i)
%          indx(coun) = i; % Getting the index/rows of values that are not possible
%          coun = coun+1;
%     end
% end
% 
% mass_rates(indx,:) = []; %Eliminating the rows

%%3. To get the positive terminus ablation
for i = 1:length(mass_rates)
    if mass_rates (i,2) == 0; % In case of stable terminus ablation
        mass_rates(i,4) = mass_rates(i,5);
    elseif mass_rates(i,2) < 0; % In case of retreating terminus ablation; denominator is always -ve, so it is reversed
        mass_rates(i,4) = (mass_rates(i,5)-mass_rates(i,2));
    else mass_rates(i,2) > 0; % In case of advancing terminus ablation; denominator is always -ve, so it is reversed
        mass_rates(i,4) = (mass_rates(i,5) - mass_rates(i,2));
    end
end
mass_rates(:,6) = mass_rates(:,4)*365; %converting to GT/year
columnHeaders = {'datenumb','nmc','obs','term_ablation','discharge','term_abl_year'};
tableWithHeaders = array2table(mass_rates, 'VariableNames', columnHeaders);
writetable(tableWithHeaders, [output_path,strcat('Term_Ablation_',GlacierName,'_',region,'.csv')], 'Delimiter', ',');

%%4. find the number of points which have negative terminus ablation; only RUN ONCE
count = 0;
inds = [];
for i = 1:length(mass_rates)
    if mass_rates(i,4)<0
        count = count+1;
        inds(count) = i;
    end
end
mass_rates(inds,3) = 0;
mass_rates(inds,4) = 0;mass_rates(inds,6) = 0;
% Plotting only the dates greater than 2013 and saving the figure
date_req = datenum('2013-01-01');
for ind_date = 1:length(mass_rates)   
    if mass_rates(ind_date,1)> date_req
        break;
    end
end
%%5. Plotting and exporting the data
bar(datetime(mass_rates(ind_date:end-1,1),'convertfrom','datenum'),mass_rates(ind_date:end-1,6),'FaceColor', [0.7 0.7 0.7]);
hold on
% low_lim = years(1);
% high_lim = years (length(years)-1);
% xticks(years(1):365*7:high_lim);
plot(datetime(mass_rates(ind_date:end-1,1),'convertfrom','datenum'),mass_rates(ind_date:end-1,5)*365,'LineWidth',1.5,'Color','red');
xlabel('Time'); ylabel('Ice Flux(Gt/year)');
grid on
ax = gca; 
xticks = get(gca,'ytick');
set(gca,'yticklabel',xticks,'TickDir','out');
% ax.YAxis.TickLabelFormat = '%d'; 
ax.FontSize = 15; 
legend('Terminus Ablation','Ice Discharge','Location','northeast','FontSize',10)
% add time 
exportgraphics(ax,[outFolderImages,'TA_updated2_',GlacierName,region,'.png'],'Resolution',300);
dates_latest = datetime(mass_rates(ind_date:end-1,1),'ConvertFrom','datenum');
columnHeaders = {'dates','discharge','term_abl_year'};
discharge_latest = mass_rates(ind_date:end-1,5)*365;
term_abl_latest = mass_rates(ind_date:end-1,6);
tableWithHeaders2 = table(dates_latest,discharge_latest,term_abl_latest, 'VariableNames', columnHeaders);
output_path_term = [output_path,'Term_Ablation_latest/'];
writetable(tableWithHeaders2, [output_path_term,strcat('Term_Ablation_latest_updated_',GlacierName,'_',region,'.csv')], 'Delimiter', ',');

%%6. Finding Error in Gate Meta
error_file = readtable([root_path,'/Mankoffs_Data/gate_err.csv']);
date_dis_err = datetime(table2array(error_file(:,1)));
Val_dis_err = table2array(error_file(:,2:end));
disch_err = zeros(size(Val_dis_err,1)-1,2); 

if ~exist('glacierID2', 'var')
    count = 0;
    for i = 1:size(Val_dis_err,2)
        count = count+1;
        if Val_dis_err(1,i) == glacierID
            disch_err(:,1) = datenum(date_dis_err(2:end,1));
            disch_err(:,2) = Val_dis_err(2:end,count); %In GT/year
        end
    end
    td = datenum(datetime(year(first_date),month(first_date):duration+month(first_date)+1,16));%generating array of 
    % middle day of each month
    td = td(2:end-1);
    disch2 = unique (disch_err,'rows'); %Eliminating duplicated rows
    C = td';
    years = datetime(mass_rates(:,1),'convertfrom','datenum');  
    for j = 1:length(years)
        err = find(disch2(:,1)>=(mass_rates(j,1)-1) & disch2(:,1)<(mass_rates(j,1)+1));
        if isempty (err)
            %linear error
            ind1 = find(disch2(:,1)<(mass_rates(j,1)-1),1);
            ind2 = find(disch2(:,1)>(mass_rates(j,1)+1),1);
            if isempty(ind1)
                ind1 = ind2;
            end
            %just do the average of two dates (initial and latest)
            discharg_err(j) = sqrt((disch2(ind1,2).^2)/4 + (disch2(ind2,2).^2)/4);
        else
            %weighted error
            ab = 0;
            for i = 1:length(err)
                ab =  ab + (disch2(err(i),2));
            end
            discharg_err(j) = sqrt(ab);
        end
    end
else
    glacierIDs = [glacierID1,glacierID2];
    disch_errs = zeros(size(Val_dis_err,1)-1,3);
    for j = 1:length(glacierIDs)
        count = 0;
        for i = 1:size(Val_dis_err,2)
            count = count+1;
            if Val_dis_err(1,i) == glacierIDs(j)
                disch_errs(:,1) = datenum(date_dis_err(2:end,1));
                disch_errs(:,j+1) = Val_dis_err(2:end,count); %In GT/year
            end
        end
    end
    disch_err = zeros(length(disch_errs),2);
    disch_err(:,1) = disch_errs(:,1);
    disch_sum = zeros(length(disch_errs),1);
    for i = 1:length(glacierIDs)
        disch_sum = disch_sum + disch_errs(:,i+1).^2;
    end
    disch_err(:,2) = sqrt(disch_sum);
    td = datenum(datetime(year(first_date),month(first_date):duration+month(first_date)+1,16));%generating array of 
    % middle day of each month
    td = td(2:end-1);
    disch2 = unique (disch_err,'rows'); %Eliminating duplicated rows
    C = td';
    years = datetime(mass_rates(:,1),'convertfrom','datenum');  
    for j = 1:length(years)
        err = find(disch2(:,1)>=(mass_rates(j,1)-1) & disch2(:,1)<(mass_rates(j,1)+1));
        if isempty (err)
            %linear error
            ind1 = find(disch2(:,1)<(mass_rates(j,1)-1),1);
            ind2 = find(disch2(:,1)>(mass_rates(j,1)+1),1);
            if isempty(ind1)
                ind1 = ind2;
            end
            %just do the average of two dates (initial and latest)
            discharg_err(j) = sqrt((disch2(ind1,2).^2)/4 + (disch2(ind2,2).^2)/4);
        else
            %weighted error
            ab = 0;
            for i = 1:length(err)
                ab =  ab + (disch2(err(i),2));
            end
            discharg_err(j) = sqrt(ab);
        end
    end
end
%% Error in bed data; Run this if you are propagating error from BedMachine
% This one takes time to run
bed_file = [root_path,'/BedTopo/BedMachineGreenland.nc'];
x_bed = double(ncread(bed_file,'x'));
y_bed = double(ncread(bed_file,'y'));
bed_err = ncread(bed_file,'errbed');
bed_err = bed_err'; 
[Rast,RA] = readgeoraster([root_path,'/BedTopo/bed_error.tif']);
% Get image info: 
R = geotiffinfo([root_path,'/BedTopo/bed_error.tif']); 
% Get x,y locations of pixels: 
[x_loc,y_loc] = pixcenters(R); 
% Convert x,y arrays to grid: 
[X,Y] = meshgrid(x_loc,y_loc);
for i = 1:length(Poly)
    polygon = Poly(i);
    x_req = Poly{1,i}.Vertices(:,1);
    y_req = Poly{1,i}.Vertices(:,2);
    x_req1 = x_req';
    y_req1 = y_req';
    mask = inpolygon(X,Y,x_req,y_req);
    %plot(X(mask),Y(mask),'r+')
    err_value = Rast(mask);
    bed_error(i) = mean(err_value);
end

%error in DEMs
arc_DEM_error = 2.28;
korsgaard_error = 5.4;
for i = 1:length(bed_error)
    thick_err(i,1) = (sqrt(bed_error(i)^2+DEM_err(i,1)^2));
    thick_err(i,2) = (sqrt(bed_error(i)^2+DEM_err(i,2)^2));
end

%%7. Error in Area
leng = length(Shp_filtered);
Satellites = strings(leng,1);
for i = 1:leng
    Satellites(i) = Shp_filtered{1,i}.Satellite;
end
Uniq_Sat = unique(Satellites);
for i=1:leng
    x_coord = Shp_filtered{1,i}.X;
    y_coord = Shp_filtered{1,i}.Y;
    total_length = 0;
    for j=1:length(x_coord)-1
        x_dist = (x_coord(j+1) - x_coord(j))^2;
        y_dist = (y_coord(j+1) - y_coord(j))^2;
        total_length = total_length + sqrt(x_dist+y_dist);
    end
    Satellites = Shp_filtered{1,i}.Satellite;
    if regexp(Satellites, regexptranslate('wildcard', 'L*')) == 1
        pix_size = 30;
    elseif regexp(Satellites, regexptranslate('wildcard', 'S*')) == 1
        pix_size = 10;
    elseif regexp(Satellites, regexptranslate('wildcard', 'A*')) == 1
        pix_size = 50;
    elseif regexp(Satellites, regexptranslate('wildcard', 'R*')) == 1;
        pix_size = 8;
    else
        pix_size = 30;
    end
    num_pixels = total_length/pix_size;
    Area_err(i) = num_pixels * pix_size^2;%in sq meters
    p(i,1) = total_length/1000;
    p(i,2) = pix_size;
end
for i = 1:leng %VolumeTS(i,3) is in kms
    Vol_err(i,1) = (abs(VolumeTS(i,4)*10^9))*sqrt((Area_err(i)./(AreaTS(i,2)*1e6)).^2+(thick_err(i,1)./(VolumeTS(i,3).*1000)).^2); % In m^3
    Vol_err(i,2) = (abs(VolumeTS(i,4)*10^9))*sqrt((Area_err(i)./(AreaTS(i,2)*1e6)).^2+(thick_err(i,2)./(VolumeTS(i,3).*1000)).^2);
end

NMC_err(:,1) = (abs(Vol_err(:,1)) * 917)/10^12; %In Gt
NMC_err(:,2) = (abs(Vol_err(:,2)) * 917)/10^12; %In gt
ds(1,1) = 1;

dates_mass = [term_mass.date;flipud(term_mass.date)];
m1 = (term_mass.terminus_mass+NMC_err(:,1));
m2 = flipud(term_mass.terminus_mass)-flipud(NMC_err(:,2));
TA_uppdown = [ m1; m2];
fill(dates_mass,TA_uppdown,[.7 .7 .7])
hold on;
plot(term_mass.date,term_mass.terminus_mass,'red')
xlabel('Time'); ylabel('Mass (Gt)');
grid on
ax = gca; 
legend('Error bound','Near Terminus Mass','Location','northeast','FontSize',10)
% add time 
exportgraphics(ax,[outFolderImages,'Near_Term_Mass_Updated_',GlacierName,region,'.png'],'Resolution',300);
upper_mass_lim = NMC_err(:,1); lower_mass_lim = NMC_err(:,2);dates = term_mass.date;terminus_mass = term_mass.terminus_mass;
term_mass_table = table(dates,terminus_mass,upper_mass_lim,lower_mass_lim);
writetable(term_mass_table, [output_path,strcat('Term_Mass_Error_Updated_',GlacierName,'_',region,'.csv')], 'Delimiter', ',');

%% Plotting a specific glacier (optional)
root_path = '/Users/amankc/Terminus_Ablation/';
outFolderImages = [root_path,'Results/Images/'];
GlacierName = glacier_name;
% region = 'NW';
term_path = '/Users/amankc/Terminus_Ablation/Results/Output/';
term_mass = readtable([term_path,strcat('Term_Mass_Error_Updated_',GlacierName,'_',region,'.csv')]);
date_req = datetime('2013-01-01','InputFormat','yyyy-MM-dd');
A = table2array(term_mass(:,1));
B = table2array(term_mass(:,2:4));
for ind_date = 1:length(A)   
    if A(ind_date) > date_req
        break;
    end
end
dates_mass = [term_mass.dates(ind_date:end);flipud(term_mass.dates(ind_date:end))];
m1 = (term_mass.terminus_mass(ind_date:end)+term_mass.upper_mass_lim((ind_date:end)));
m2 = flipud(term_mass.terminus_mass(ind_date:end))-flipud(term_mass.lower_mass_lim(ind_date:end));
TA_uppdown = [ m1; m2];
fill(dates_mass,TA_uppdown,[.7 .7 .7])
hold on;
plot(term_mass.dates(ind_date:end),term_mass.terminus_mass(ind_date:end),'black','linewidth',1)
xlabel('Time'); ylabel('Mass (Gt)');
grid on
ax = gca; 
ax.FontSize = 14;
legend('Error bound','Near Terminus Mass','Location','northeast','FontSize',10)
exportgraphics(ax,[outFolderImages,'Latest_NTM_Updated',GlacierName,region,'.png'],'Resolution',300);


%%8. Average_mass_rates
varnames2 = {'dates','datenumb','mass','error_upp','error_down'};
TPa = table(term_mass_table.dates,datenum(term_mass_table.dates),term_mass_table.terminus_mass,NMC_err(:,1),NMC_err(:,2),'VariableNames',varnames2);
[C,ia] = unique(TPa.datenumb);
TP1 = TPa(ia,:);
for i = 2:length(t)-2
    clear mass_rates err_upp err_down;
    terms = find(TP1.datenumb>=t(i) & TP1.datenumb<t(i+1));
    num_term = length(terms);
    if isempty (terms) %if no data is found between the month, lineraly interpolate 
        % between two earliest and latest dates from the month
        mean_date = round(mean([t(i) t(i+1)])); %getting the middle date of the month
        average_err(i,1) = mean_date;
        %mass_rates(j,1) = datetime(mean_date,'ConvertFrom','datenum') (To
        %convert into regular dates) Make sure predefined arrays will
        %accept datetime values
        [mass_rates,err_upp,err_down] = Linear_mass_rates(i,TP1,t);
        average_err(i,2) = mass_rates;
        average_err(i,3) = err_upp; average_err(i,4) = err_down;
    else
%         mass_rates (j-1,1) = round(mean([t(i) t(j+1)]));
%         month_dates = t(j+1)-t(j);
        mean_date = round(mean([t(i) t(i+1)])); %getting the middle date of the month
        average_err(i,1) = mean_date;
        [mass_rates,err_upp,err_down]  = Average_Mass_Rates(terms,TP1,t);
        average_err(i,2) = mass_rates;
        average_err(i,3) = err_upp; average_err(i,4) = err_down;
    end
end
date = datetime(average_err(2:end,1),'convertfrom','datenum');
avg_nmc_rate = average_err(2:end,2)*365;
max_rate = average_err(2:end,3)*365; 
min_rate = average_err(2:end,4)*365;
average_rate_table = table(date,avg_nmc_rate,max_rate,min_rate);
output_path_avg_rate = [output_path,'Average_Rate/'];
writetable(average_rate_table, [output_path_avg_rate,strcat('Average_Rate_Error_',GlacierName,'_',region,'.csv')], 'Delimiter', ',');

%%9. plotting the data
clear term_abl upper_bound lower_bound;
dates_avg = date;
dates_avg2 = [dates_avg;flipud(dates_avg)];
ma = max_rate;
mb = flipud(min_rate);
TA_uppdown2 = [ ma; mb];
fill(dates_avg2,TA_uppdown2,[.7 .7 .7])
hold on;
plot(dates_avg,avg_nmc_rate,'red')
xlabel('Time'); ylabel('Ice Flux(Gt/year)');
grid on
ax = gca; 
xticks = get(gca,'ytick');
set(gca,'yticklabel',xticks);
ax.YAxis.TickLabelFormat = '%d'; 
ax.FontSize = 15; 
legend('Error bound','Terminus Ablation','Location','northwest','FontSize',10)
exportgraphics(ax,[outFolderImages,'TA_Average_Error_',GlacierName,'.png'],'Resolution',300);
discharge2 = discharge(2:end)*365;
for i = 1:length(avg_nmc_rate)
    if avg_nmc_rate(i) == 0; % In case of stable terminus ablation
        term_abl(i) = discharge2(i);
    elseif avg_nmc_rate(i) < 0; % In case of retreating terminus ablation; denominator is always -ve, so it is reversed
        term_abl(i)  = (discharge2(i)-avg_nmc_rate(i));
    else avg_nmc_rate(i) > 0; % In case of advancing terminus ablation; denominator is always -ve, so it is reversed
        term_abl(i)  = (discharge2(i) - avg_nmc_rate(i));
    end
end

for i=1:length(max_rate)
    upper_bound(i) = sqrt(max_rate(i)^2+discharg_err(i+1)^2);
    lower_bound(i) = sqrt(min_rate(i)^2+discharg_err(i+1)^2);
end
term_abl = term_abl'; upper_bound = upper_bound'; lower_bound = lower_bound';
average_term_abl_table = table(dates_avg,term_abl,upper_bound,lower_bound);
output_path_avg_rate = [output_path,'Average_Rate/'];
writetable(average_term_abl_table, [output_path_avg_rate,strcat('Average_Term_Abl_Error_',GlacierName,'_',region,'.csv')], 'Delimiter', ',');

