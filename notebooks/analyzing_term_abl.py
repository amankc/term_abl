#importing the required libraries
import geopandas as gpd
import pandas as pd
import osgeo as osg
import numpy as np
import matplotlib.pyplot as plt
from osgeo import ogr
import os
import pyproj
import datetime
import cartopy.crs as ccrs
import scipy
import glob

#This section plots terminus ablation vales for all of the glaciers in the region with lighter values indicating individual glacier and darker color is the mean.
"""
CE = #00B1D8
CE_light = #B7FCFF
SE = #7900FF
SE_light = #E7D6FA
SW = #FF8403
SW_light = #F9DBBB
CW = #ff41d6
CW_light = #FFB8BF
NW = #014E05
NW_light = #A1FC77
NE = #848484
"""
root_path = '/Users/amankc/Terminus_Ablation/';
wild_card = '_NW'; #change with the regions you want to plot
region = 'NW';

if region == 'CE':
    color_code1 = '#B7FCFF'; # This is the light color
    color_code2  = '#00B1D8';
elif region == 'SE':
    color_code1 = '#E7D6FA'; # This is the light color
    color_code2  = '#7900FF';
elif region == 'SW':
    color_code1 = '#F9DBBB'; # This is the light color
    color_code2  = '#FF8403';
elif region == 'CW':
    color_code1 = '#FFB8BF'; # This is the light color
    color_code2  = '#ff41d6';
elif region == 'NW':
    color_code1 = '#A1FC77'; # This is the light color
    color_code2  = '#014E05';
else:
    color_code1 = '#848484'; # This is the light color
    color_code2  = '#848484';


directory_path = root_path + 'Results/Output/Term_Ablation_Latest/'
file_pattern = 'Term_Ablation_latest_updated_*' + region +'*.csv'  # Replace this with your desired file pattern

# Use glob to find files matching the pattern in the directory
matching_files = glob.glob(directory_path + file_pattern)

csv_file = [];
for file in matching_files:
    csv_file.append(pd.read_csv(file))

term_ablation2 = [];
discharge = [];
date_axis = [];
l = [];
for i in range(len(csv_file)):
    l.append(len(csv_file[i]));
len_max = np.min(l)
min_idx = pd.Series(l).idxmin()

#just to make the length of vectors equal   
for date in csv_file[min_idx].dates:
    date_axis.append(date)
for i in range(len(csv_file)):
    j = 0;
    term_ablation2.append([])# Initialize an empty list for each index in term_ablation2
    discharge.append([])
    for j in range(len_max):
        term_ablation2[i].append(csv_file[i].term_abl_year[j])
        discharge[i].append(csv_file[i].discharge[j])
        # date_axis[i].append(csv_file[i].dates[j]);
        # date_axis[i] = csv_file[i].dates[j]
        j+=1;
mean_term_abl = [];
mean_discharge = [];
for i in range(len(term_ablation2[0])):
    counter = 0;
    sum_term_abl = 0;
    sum_discharge  = 0;
    for j in range(len(term_ablation2)):
        sum_term_abl = sum_term_abl + term_ablation2[j][i]
        sum_discharge = sum_discharge + discharge[j][i]
        counter += 1;
    mean_term_abl.append(sum_term_abl/counter);
    mean_discharge.append(sum_discharge/counter);
plt.figure(figsize=(14,10)) 
ax = plt.plot(term_ablation2[0], color = color_code1)
for i in range(len(term_ablation2)):
    plt.plot(term_ablation2[i], color = color_code1)
plt.plot(date_axis,mean_term_abl, color = color_code2)
plt.xlabel('Years',fontsize = 14)
plt.ylabel('Terminus Ablation (Gt/year)',fontsize = 14)
plt.title('Time series for ' + region)
xz = 0;a=[];
interval = 12;
dates_ticks = [];
dates_label = [];
for date in date_axis:
    a.append(xz%interval)
    if xz%interval == 0:
        dates_ticks.append(date)
        dates_label.append(date[7:])     
    xz+=1
ax = plt.gca()
plt.xticks(dates_ticks,dates_label,fontsize = 13,weight='bold');
plt.yticks(fontsize = 13,weight='bold');
plt.tight_layout()
plt.grid('on',color = '#E4E1E1')
plt.savefig(root_path + 'Results/Images/Updated_Regional_Figures/'+'Time_series_for ' + region +'.png',dpi=300,bbox_inches='tight', pad_inches=0)

#This plots normailzed terminus ablation and discharge at the same time
root_path = '/Users/amankc/Terminus_Ablation/';
regions_name = ['NW','NE','CW','CE','SE','SW']
fig,ax = plt.subplots(3,2,sharex='all', sharey='all',figsize=(14, 12))
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.tight_layout()
for region in regions_name:
    if region == 'CE':
        color_code1 = '#B7FCFF'; # This is the light color
        color_code2  = '#00B1D8';x=1;y=1;
    elif region == 'SE':
        color_code1 = '#E7D6FA'; # This is the light color
        color_code2  = '#7900FF';x=2;y=1;
    elif region == 'SW':
        color_code1 = '#F9DBBB'; # This is the light color
        color_code2  = '#FF8403';x=2;y=0;
    elif region == 'CW':
        color_code1 = '#FFB8BF'; # This is the light color
        color_code2  = '#ff41d6';x=1;y=0;
    elif region == 'NW':
        color_code1 = '#A1FC77'; # This is the light color
        color_code2  = '#014E05'; x=0;y=0;
    else:
        color_code1 = '#848484'; # This is the light color
        color_code2  = '#848484';x=0;y=1;
    
    import glob
    directory_path = root_path + 'Results/Output/Term_Ablation_Latest/'
    file_pattern = 'Term_Ablation_latest_updated_*' + region +'*.csv'  # Replace this with your desired file pattern
    
    # Use glob to find files matching the pattern in the directory
    matching_files = glob.glob(directory_path + file_pattern)
    
    csv_file = [];
    for file in matching_files:
        csv_file.append(pd.read_csv(file))
    
    term_ablation2 = [];
    discharge = [];
    date_axis = [];
    l = [];
    for i in range(len(csv_file)):
        l.append(len(csv_file[i]));
    len_max = np.min(l)
    min_idx = pd.Series(l).idxmin()
    #just to make the length of vectors equal   
    for date in csv_file[min_idx].dates:
        date_axis.append(date)
    for i in range(len(csv_file)):
        j = 0;
        term_ablation2.append([])# Initialize an empty list for each index in term_ablation2
        discharge.append([])
        for j in range(len_max):
        # range(len(csv_file[i].term_abl_year)):
            term_ablation2[i].append(csv_file[i].term_abl_year[j])
            discharge[i].append(csv_file[i].discharge[j])
            # date_axis[i].append(csv_file[i].dates[j]);
            # date_axis[i] = csv_file[i].dates[j]
            j+=1;
    mean_term_abl = [];
    mean_discharge = [];
    for i in range(len(term_ablation2[0])):
        counter = 0;
        sum_term_abl = 0;
        sum_discharge  = 0;
        for j in range(len(term_ablation2)):
            sum_term_abl = sum_term_abl + term_ablation2[j][i]
            sum_discharge = sum_discharge + discharge[j][i]
            counter += 1;
        mean_term_abl.append(sum_term_abl/counter);
        mean_discharge.append(sum_discharge/counter);
    
    # new step
    csvs_df = pd.DataFrame()
    for i, file in enumerate(matching_files):
        df = pd.read_csv(file)
        df['site_number'] = i
        csvs_df = pd.concat([csvs_df, df])
    csvs_df.reset_index(drop=True, inplace=True)
    # convert dates to pandas.datetime
    csvs_df['dates'] = pd.to_datetime(csvs_df['dates'])
    
    # grab unique dates
    unique_dates = csvs_df['dates'].drop_duplicates()
    discharge_mean = [csvs_df.loc[csvs_df['dates']==date]['discharge'].mean() for date in unique_dates]
    term_ablation_mean = [csvs_df.loc[csvs_df['dates']==date]['term_abl_year'].mean() for date in unique_dates]
    
    #new step
    # Normailze the figure
    term_ablation2_normalized = [];
    discharge_normalized = [];
    for i in range(len(csv_file)):
        j = 0;
        term_ablation2_normalized.append([]) # Initialize an empty list for each index in term_ablation2
        discharge_normalized.append([]);
        mean_value = np.mean(csv_file[i].term_abl_year)
        max_value = np.max(csv_file[i].term_abl_year-mean_value)
    
        mean_value_dis = np.mean(csv_file[i].discharge)
        max_value_dis = np.max(csv_file[i].discharge-mean_value_dis)
        for j in range(len_max):
        # range(len(csv_file[i].term_abl_year)):
            term_ablation2_normalized[i].append((csv_file[i].term_abl_year[j]-mean_value)/max_value)
            discharge_normalized[i].append((csv_file[i].discharge[j]-mean_value_dis)/max_value_dis)
            # date_axis[i].append(csv_file[i].dates[j]);
            # date_axis[i] = csv_file[i].dates[j]
            j+=1;
    
    mean_term_abl_normalized = [];
    mean_discharge_normalized = [];
    for i in range(len(term_ablation2_normalized[0])):
        counter = 0;
        sum_term_abl = 0; discharge_sum = 0;
        for j in range(len(term_ablation2)):
            sum_term_abl = sum_term_abl + term_ablation2_normalized[j][i]
            discharge_sum = discharge_sum +  discharge_normalized[j][i]
            counter += 1;
        mean_term_abl_normalized.append(sum_term_abl/counter);
        mean_discharge_normalized.append(discharge_sum/counter);
    columns = {'date': date_axis,'normalized_term_abl': mean_term_abl_normalized,'normalized_discharge':mean_discharge_normalized}
    mean_df = pd.DataFrame(data=columns)
    mean_df.to_csv(root_path + 'Results/Output/Normalized_Terminus_Ablation/'+'Normalized_Time_series_for_' + region + '.csv')
    
    
    date_array = np.array(date_axis)
    mean_ta_array = np.array(mean_term_abl_normalized)
    term_abl_df = pd.DataFrame({'date': date_array,'mean_ta_normalized': mean_ta_array})
    # for i in range(len(term_ablation2_normalized)):
    #     term_abl_array[i] = np.array(term_ablation2_normalized[0])
    plt.figure(figsize=(14,8)) 
    term_abl_df['date'] = pd.to_datetime(term_abl_df['date'])
    date_array = np.array(date_axis)
    mean_dis_array = np.array(mean_discharge_normalized)
    dis_df = pd.DataFrame({'date': date_array,'mean_dis_normalized': mean_discharge_normalized})
    
    # Normailze the figure
    term_ablation2_normalized = [];
    discharge2_normalized = [];
    for i in range(len(csv_file)):
        j = 0;
        term_ablation2_normalized.append([]) # Initialize an empty list for each index in term_ablation2
        discharge2_normalized.append([]);
        mean_value = np.mean(csv_file[i].term_abl_year)
        max_value = np.max(csv_file[i].term_abl_year-mean_value)
    
        mean_value_dis = np.mean(csv_file[i].discharge)
        max_value_dis = np.max(csv_file[i].discharge-mean_value_dis)
        for j in range(len_max):
        # range(len(csv_file[i].term_abl_year)):
            term_ablation2_normalized[i].append((csv_file[i].term_abl_year[j]-mean_value)/max_value)
            discharge2_normalized[i].append((csv_file[i].discharge[j]-mean_value)/max_value)
            # date_axis[i].append(csv_file[i].dates[j]);
            # date_axis[i] = csv_file[i].dates[j]
            j+=1;
    
    mean_term_abl_normalized = [];
    mean_discharge_normalized2 = [];
    for i in range(len(term_ablation2_normalized[0])):
        counter = 0;
        sum_term_abl = 0; discharge_sum = 0;
        for j in range(len(term_ablation2)):
            sum_term_abl = sum_term_abl + term_ablation2_normalized[j][i]
            discharge_sum = discharge_sum +  discharge2_normalized[j][i]
            counter += 1;
        mean_term_abl_normalized.append(sum_term_abl/counter);
        mean_discharge_normalized2.append(discharge_sum/counter);
    term_abl_color = color_code2
    discharge_color = '#ff0000ff'
    ax[x][y].plot(term_abl_df['date'],mean_term_abl_normalized, color = term_abl_color)
    ax[x][y].plot(term_abl_df['date'],mean_discharge_normalized2, color = discharge_color)
    ax[x][y].set_xlim(np.datetime64('2013-01-01'), np.datetime64('2023-03-31'))
    ax[x][y].set_ylim(-0.45,1)
    ax[x][y].grid('on',color = '#E4E1E1')
    ax[x][y].legend(['Normalized terminus ablation','Normalized discharge'],loc = 'best')
# ax2 = plt.gca()
for row in ax:
    for subplot in row:
        # Set parameters for each subplot
        subplot.tick_params(labelsize=14)

fig.savefig(root_path + 'Results/Images/Regional_Figures/'+'Normalized_discharge_term-ablation_updated'+'.png',dpi=300)
plt.show()

#This plots the actual magnitude of terminus ablation and discharge with the same y-axis
#with same y-axis
root_path = '/Users/amankc/Terminus_Ablation/'
regions_name = ['NW', 'NE', 'CW', 'CE', 'SE', 'SW']
output_path = os.path.join(root_path, 'Results/Output/Term_Ablation_Latest/')
total_output_path = os.path.join(root_path, 'Results/Output/Total_Terminus_Ablation/')
fig, ax = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(14, 12))
# plt.subplots_adjust(wspace=0.1, hspace=0.1)

plt.tight_layout()

color_codes = {
    'CE': ('#B7FCFF', '#00B1D8', 1, 1),
    'SE': ('#E7D6FA', '#7900FF', 2, 1),
    'SW': ('#F9DBBB', '#FF8403', 2, 0),
    'CW': ('#FFB8BF', '#ff41d6', 1, 0),
    'NW': ('#A1FC77', '#014E05', 0, 0),
    'NE': ('#848484', '#848484', 0, 1)
}

for region, (color_code1, color_code2, x, y) in color_codes.items():
    file_pattern = f'Term_Ablation_latest_updated_*{region}*.csv'
    matching_files = glob.glob(os.path.join(output_path, file_pattern))
    
    csv_files = [pd.read_csv(file) for file in matching_files]

    # Ensure all CSV files have the same length by trimming
    min_length = min(len(df) for df in csv_files)
    csv_files = [df.head(min_length) for df in csv_files]

    date_axis = csv_files[0]['dates']
    
    # Compute the total of term_ablation and discharge
    total_term_abl = np.sum([df['term_abl_year'].values for df in csv_files], axis=0)
    total_discharge = np.sum([df['discharge'].values for df in csv_files], axis=0)

    # Create a DataFrame to save the total data
    total_df = pd.DataFrame({
        'date': date_axis,
        'total_term_abl': total_term_abl,
        'total_discharge': total_discharge
    })
    # total_df.to_csv(os.path.join(total_output_path, f'Total_Time_series_for_{region}.csv'), index=False)
    
    # Plotting
    ax[x, y].plot(pd.to_datetime(date_axis), total_term_abl, color=color_code2, linestyle = '--',label='Total terminus ablation')
    ax[x, y].plot(pd.to_datetime(date_axis), total_discharge, color=color_code2, label='Total discharge')
    ax[x, y].set_ylim(0,140)
    ax[x, y].set_xlim(np.datetime64('2013-01-01'), np.datetime64('2023-03-31'))
    ax[x, y].grid(True, color='#E4E1E1')
    # ax[x, y].legend(loc='upper left')

for subplot in ax.flatten():
    subplot.tick_params(labelsize=16)

for row in ax:
    for subplot in row:
        # Set parameters for each subplot
        subplot.tick_params(labelsize=14)
        for tick in subplot.get_xticklabels():
            tick.set_fontsize(14)
            tick.set_fontweight('bold')
        for tick in subplot.get_yticklabels():
            tick.set_fontsize(14)
            tick.set_fontweight('bold')

fig.savefig(root_path + 'Results/Images/Regional_Figures/'+'Total_discharge_term-ablation_updated2'+'.png',dpi=300)
plt.show()

#This section plots the montly averaged terminus ablation and discharge with standard deviation as error bounds.
root_path = '/Users/amankc/Terminus_Ablation/'
regions_name = ['NW', 'NE', 'CW', 'CE', 'SE', 'SW']
output_path = os.path.join(root_path, 'Results/Output/Term_Ablation_Latest/')
total_output_path = os.path.join(root_path, 'Results/Output/Total_Terminus_Ablation/')
plt.figure(figsize=(10, 8))
color_codes = {
    'CE': ('#B7FCFF', '#00B1D8', 1, 1),
    'SE': ('#E7D6FA', '#7900FF', 2, 1),
    'SW': ('#F9DBBB', '#FF8403', 2, 0),
    'CW': ('#FFB8BF', '#ff41d6', 1, 0),
    'NW': ('#A1FC77', '#014E05', 0, 0),
    'NE': ('#848484', '#848484', 0, 1)
}

def get_season(month):
    '''Function to map month to season'''
    if month in [12, 1, 2]:
        return 'Winter'
    elif month in [3, 4, 5]:
        return 'Spring'
    elif month in [6, 7, 8]:
        return 'Summer'
    elif month in [9, 10, 11]:
        return 'Fall'

for region, (color_code1, color_code2, x, y) in color_codes.items():
    file_pattern = f'Term_Ablation_latest_updated_*{region}*.csv'
    matching_files = glob.glob(os.path.join(output_path, file_pattern))
    
    csv_files = [pd.read_csv(file) for file in matching_files]

    # Ensure all CSV files have the same length by trimming
    min_length = min(len(df) for df in csv_files)
    csv_files = [df.head(min_length) for df in csv_files]

    date_axis = csv_files[0]['dates']
    
    # Compute the total of term_ablation and discharge
    total_term_abl = np.sum([df['term_abl_year'].values for df in csv_files], axis=0)
    total_discharge = np.sum([df['discharge'].values for df in csv_files], axis=0)

    # Create a DataFrame to save the total data
    total_df = pd.DataFrame({
        'date': pd.to_datetime(date_axis, format='%d-%b-%Y'),
        'total_term_abl': total_term_abl,
        'total_discharge': total_discharge
    })

    # Add a season column
    total_df['season'] = total_df['date'].dt.month.apply(get_season)
    total_df['year'] = total_df['date'].dt.year
    total_df['month'] = total_df['date'].dt.month
    monthly_avgs = total_df.groupby('month').mean(numeric_only=True)

    seasonal_sums = total_df.groupby('season').mean(numeric_only=True)
    sd = total_df.groupby('month').std(numeric_only=True)
    
    ta = monthly_avgs['total_term_abl'];
    error_bound = sd['total_term_abl'];
    # print(region)
    # print(monthly_avgs)
    # Group by year and season and calculate the sum for each season
    # seasonal_totals = total_df.groupby(['year', 'season']).sum().reset_index()

    # seasonal_totals.to_csv(os.path.join(total_output_path, f'Seasonal_Totals_for_{region}.csv'), index=False)
    
    # Plotting
    a = monthly_avgs['total_term_abl']
    plt.plot(ta, color=color_code2, label='Total terminus ablation',linewidth = 3)
    # ax[x, y].plot(monthly_avgs['total_discharge'], color='#ff0000ff', label='Total discharge')
    # ax[x, y].set_xlim(np.datetime64('2013-01-01'), np.datetime64('2023-03-31'))
    
    # plt.plot(x, y, 'k-')
    
    list_of_months = np.arange(1, 13, 1)
    plt.fill_between(list_of_months,ta-error_bound, ta+error_bound,color=color_code1,alpha= 0.3)
    # ax[x, y].legend(loc='upper left')
    del seasonal_sums, ta, error_bound
plt.grid(True, color='#E4E1E1')
plt.ylabel('average terminus ablation (Gt/year)',fontsize = 14)
months = list(calendar.month_name)[1:]
name = []
for month in months:
    name.append(month[:3])
plt.xlabel('months',fontsize = 14)
plt.xticks(list_of_months,name,fontsize = 14,weight = 'bold')
plt.yticks (fontsize = 14,weight = 'bold')

plt.tight_layout()
# plt.legend()
plt.savefig(root_path + 'Results/Images/Regional_Figures/'+'average_term-ablation_updated'+'.png', dpi=300)
plt.show() 
