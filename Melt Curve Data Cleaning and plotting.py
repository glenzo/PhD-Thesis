## Melt Curve Data Cleaning and plotting.
import pandas as pd

import matplotlib.pyplot as plt
import numpy as np

# Read the Excel file
file_name = "Melt_curve.xlsx"
sheet_name = "Melt Curve Raw Data"
df1 = pd.read_excel(file_name, sheet_name=sheet_name, engine='openpyxl')

# Group the data by 'Target Name' and 'Reading' and calculate the average fluorescence and derivative
grouped_df = df1.groupby(['Target Name', 'Reading'])[['Fluorescence', 'Derivative']].mean().reset_index()

# Save the result to a new Excel file
grouped_df.to_excel("Grouped_data.xlsx", index=False)

print("The data has been grouped by 'Target Name' and 'Reading', and the average fluorescence and derivative values have been calculated.")


# Read the Excel file
file_name = "Grouped_data.xlsx"
df = pd.read_excel(file_name, engine='openpyxl')

# Create the temperature variable
temperature = np.linspace(4.0, 95.0, 911)

# Get unique target names
target_names = df['Target Name'].unique()

# Set the style to use
plt.style.use('default')

# Create a figure
fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))

# Set the font size for the entire plot
font_size = 14
plt.rc('font', size=font_size, family='Helvetica')

# Define the colors to use
colors = plt.cm.viridis(np.linspace(0, 1, len(target_names)))
colors = ['#E24A33', '#348ABD', '#988ED5', '#777777', '#FBC15E', '#8EBA42', '#FFB5B8', '#FF8101', '#16A085', '#0B486B']

colors = plt.cm.nipy_spectral(np.linspace(0, 1, len(target_names)))

# Plot Temperature versus Fluorescence and Derivative for each target name
for target, color in zip(target_names, colors):
    target_df = df[df['Target Name'] == target]
    fluorescence = target_df['Fluorescence'].values
    derivative = target_df['Derivative'].values

    # Make sure there are 911 values corresponding to the 911 temperature values
    if len(fluorescence) == 911 and len(derivative) == 911:
        ax1.plot(temperature, fluorescence, label=target, linewidth=2, color=color)
        ax2.plot(temperature, derivative, label=target, linewidth=2, color=color)
    else:
        print(f"Skipping {target} due to insufficient data points.")

# Customize the plot
ax1.set_xlabel('Temperature (°C)')
ax1.set_ylabel('Fluorescence (a.u.)')
# ax1.set_title('Temperature vs Fluorescence')

ax2.set_xlabel('Temperature (°C)')
ax2.set_ylabel('Derivative (a.u.)')
# ax2.set_title('Temperature vs Derivative')

# Turn off the grid
ax1.grid(False)
ax2.grid(False)

# Set tick marks
ax1.xaxis.set_ticks(np.arange(0, 100, 10))
ax2.xaxis.set_ticks(np.arange(0, 100, 10))

# Set the legend outside the plot
legend = ax1.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=font_size-2, frameon=True, framealpha=1)
legend.get_frame().set_edgecolor('black')

legend2 = ax2.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=font_size-2, frameon=True, framealpha=1)
legend2.get_frame().set_edgecolor('black')

# Set tick labels and sizes
ax1.tick_params(axis='both', direction='in', length=6, width=1, labelsize=font_size-2)
ax2.tick_params(axis='both', direction='in', length=6, width=1, labelsize=font_size-2)


plt.tight_layout()

# Save the plot as a high-resolution image
fig.savefig('plot_no_grid.png', dpi=600, bbox_inches='tight')

# Show the plot
plt.show()

####

## Assuming df1 is already defined

# Get unique target names and readings
target_names = df1['Target Name'].unique()
readings = df1['Reading'].unique()

# Calculate the number of readings per target
num_readings_per_target = len(df1) // len(target_names)

# Create empty dataframes for each group of readings
df1_A = pd.DataFrame(columns=df1.columns)
df1_B = pd.DataFrame(columns=df1.columns)
df1_C = pd.DataFrame(columns=df1.columns)

# Split the dataframe into three dataframes based on the number of readings per target name
for target in target_names:
    target_data = df1[df1['Target Name'] == target]
    for i in range(0, num_readings_per_target, 911):
        if i == 0:
            df1_A = pd.concat([df1_A, target_data[i:i+911]], ignore_index=True)
        elif i == 911:
            df1_B = pd.concat([df1_B, target_data[i:i+911]], ignore_index=True)
        elif i == 1822:
            df1_C = pd.concat([df1_C, target_data[i:i+911]], ignore_index=True)

# Display the three new dataframes
print(df1_A)
print(df1_B)
print(df1_C)


###

def get_min_max_temp(df, target_names, temperature_range, derivative_extrema):
    result_df = pd.DataFrame(columns=['Target Name', 'Temperature', 'Derivative'])

    for target in target_names:
        target_data = df[df['Target Name'] == target]
        range_df = target_data[(target_data['Temperature'] >= temperature_range[0]) & (target_data['Temperature'] <= temperature_range[1])]

        if not range_df.empty:
            if derivative_extrema == 'min':
                extrema_row = range_df.loc[range_df['Derivative'].idxmin()]
            elif derivative_extrema == 'max':
                extrema_row = range_df.loc[range_df['Derivative'].idxmax()]
            result_df = result_df.append(extrema_row, ignore_index=True)

    return result_df

# Assuming df1_A, df1_B, and df1_C are already defined
unique_target_names = df1['Target Name'].unique()

# Get the temperature values for the minimum and maximum derivative values for each target name in each dataframe
min_temp_A = get_min_max_temp(df1_A, unique_target_names, (15, 30), 'min')
max_temp_A = get_min_max_temp(df1_A, unique_target_names, (50, 70), 'max')

min_temp_B = get_min_max_temp(df1_B, unique_target_names, (15, 30), 'min')
max_temp_B = get_min_max_temp(df1_B, unique_target_names, (50, 70), 'max')

min_temp_C = get_min_max_temp(df1_C, unique_target_names, (15, 30), 'min')
max_temp_C = get_min_max_temp(df1_C, unique_target_names, (50, 70), 'max')

# Export to Excel
with pd.ExcelWriter('Min_Max_Derivative_Temperatures.xlsx') as writer:
    min_temp_A.to_excel(writer, sheet_name='Dataset A - Min Temp', index=False)
    max_temp_A.to_excel(writer, sheet_name='Dataset A - Max Temp', index=False)

    min_temp_B.to_excel(writer, sheet_name='Dataset B - Min Temp', index=False)
    max_temp_B.to_excel(writer, sheet_name='Dataset B - Max Temp', index=False)

    min_temp_C.to_excel(writer, sheet_name='Dataset C - Min Temp', index=False)
    max_temp_C.to_excel(writer, sheet_name='Dataset C - Max Temp', index=False)



# Function to calculate standard deviation across datasets for min/max values
def calculate_std_dev(min_temp_A, min_temp_B, min_temp_C):
    target_names = min_temp_A['Target Name'].unique()
    std_dev_list = []

    for target_name in target_names:
        temp_A = min_temp_A.loc[min_temp_A['Target Name'] == target_name]['Temperature'].values[0]
        temp_B = min_temp_B.loc[min_temp_B['Target Name'] == target_name]['Temperature'].values[0]
        temp_C = min_temp_C.loc[min_temp_C['Target Name'] == target_name]['Temperature'].values[0]

        std_dev = np.std([temp_A, temp_B, temp_C])

        std_dev_list.append({'Target Name': target_name, 'Std Dev': std_dev})

    return pd.DataFrame(std_dev_list)

# Calculate standard deviation for min and max temperature values
min_temp_std_dev = calculate_std_dev(min_temp_A, min_temp_B, min_temp_C)
max_temp_std_dev = calculate_std_dev(max_temp_A, max_temp_B, max_temp_C)

# Rename columns
min_temp_std_dev.rename(columns={'Std Dev': 'Tm* Std Dev'}, inplace=True)
max_temp_std_dev.rename(columns={'Std Dev': 'Tm Std Dev'}, inplace=True)

# Merge the dataframes
temp_std_dev = min_temp_std_dev.merge(max_temp_std_dev, on='Target Name')

# Save the standard deviation values in a new sheet in the same Excel file
with pd.ExcelWriter('Min_Max_Derivative_Temperatures.xlsx', mode='a') as writer:
    temp_std_dev.to_excel(writer, sheet_name='Temperature Std Dev', index=False)

###