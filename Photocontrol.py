# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns

# Replace 'your_excel_file.xlsx' with the actual name of your Excel file
# excel_file = 'MB-RC 1 min exposures.xlsx'
# excel_file = 'MB-RC 1 min exposures.xlsx'
# Sheet2="Python"
# Sheet="Control"

# # Read the Excel file into a pandas DataFrame
# df = pd.read_excel(excel_file, sheet_name=Sheet, engine='openpyxl')
# df2 = pd.read_excel(excel_file, sheet_name=Sheet2, engine='openpyxl')

# # Initialize the cumulative time variable and gap between experiments
# cumulative_time = 0

# gap_between_experiments = 1

# # Assuming there are 5 experiments and their data are stored in pairs of columns (Time [min], Int)
# for i in range(5):
#     time_col_idx = 2 * i
#     int_col_idx = 2 * i + 1

#     # Check if the required columns exist in the DataFrame
#     if time_col_idx < len(df.columns) and int_col_idx < len(df.columns):
#         time_col = df.columns[time_col_idx]
#         int_col = df.columns[int_col_idx]

#         # Calculate the cumulative time for the current experiment
#         cumulative_time_experiment = df[time_col] + cumulative_time

#         # Update the cumulative time for the next experiment
#         previous_cumulative_time = cumulative_time
#         cumulative_time += df[time_col].max() + gap_between_experiments

#         plt.plot(cumulative_time_experiment, df[int_col], label=f'Experiment {i+1}', linewidth=2)

#         # Add blue lines for every second of the 1-minute gap at the end of the previous experiment
#         if i > 0:
#             for t in range(int(previous_cumulative_time), int(previous_cumulative_time + gap_between_experiments)):
#                 plt.axvline(t, color='blue', linestyle='--', linewidth=0.7, alpha=0.3)
#     else:
#         print(f'Could not find required columns for Experiment {i+1}')

# # Customize plot
# plt.xlabel('Cumulative Time [min]', fontsize=14)
# plt.ylabel('Int', fontsize=14)
# plt.title('Cumulative Time vs Int for 5 Experiments with 1-min gaps', fontsize=16)
# plt.legend(fontsize=12)

# # Remove the top and right spines for a cleaner look
# sns.despine()

# # Display plot
# plt.show()

# # Save the plot as a high-resolution image
# plt.savefig('MB-RC_New-Control_Photocontrol_concatenated.svg', dpi=600, bbox_inches='tight')

# # Apply a professional style to the plot
# sns.set(style="ticks", font_scale=1.2)
# plt.rcParams['figure.figsize'] = [20, 4]

# # Create the subplot with 1 row and 5 columns
# fig, axes = plt.subplots(1, 5, sharey=True)

# # Initialize the cumulative time variable and gap between experiments
# cumulative_time = 0
# gap_between_experiments = 1

# # Assuming there are 5 experiments and their data are stored in pairs of columns (Time [min], Int)
# for i in range(5):
#     time_col_idx = 2 * i
#     int_col_idx = 2 * i + 1

#     # Check if the required columns exist in the DataFrame
#     if time_col_idx < len(df.columns) and int_col_idx < len(df.columns):
#         time_col = df.columns[time_col_idx]
#         int_col = df.columns[int_col_idx]

#         # Calculate the cumulative time for the current experiment
#         cumulative_time_experiment = df[time_col] + cumulative_time

#         # Update the cumulative time for the next experiment
#         cumulative_time += df[time_col].max() + gap_between_experiments

#         axes[i].plot(cumulative_time_experiment, df[int_col], linewidth=2)
#         axes[i].set_xlabel('Cumulative Time [min]', fontsize=12)
#         axes[i].set_title(f'Experiment {i+1}', fontsize=14)

#         if i == 0:
#             axes[i].set_ylabel('Int', fontsize=12)

#     else:
#         print(f'Could not find required columns for Experiment {i+1}')

# # Adjust the space between subplots
# plt.subplots_adjust(wspace=0.3)

# # Overall title
# fig.suptitle('Cumulative Time vs Int for 5 Experiments', fontsize=16, y=1.05)

# # Remove the top and right spines for a cleaner look
# sns.despine()

# #Display plot
# plt.show()

# # Save the plot as a high-resolution image
# fig.savefig('MB-RC_Control_Photocontrol_Subplots.svg', dpi=600, bbox_inches='tight')
# ##############################################################################################################

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Replace 'your_excel_file.xlsx' with the actual name of your Excel file
excel_file = 'MB-RC 1 min exposures.xlsx'
sheet1 = "Python"
sheet2 = "Control"

# Read the Excel file into pandas DataFrames
df1 = pd.read_excel(excel_file, sheet_name=sheet1, engine='openpyxl')
df2 = pd.read_excel(excel_file, sheet_name=sheet2, engine='openpyxl')

def plot_experiments(df, title):
    # Initialize the cumulative time variable and gap between experiments
    cumulative_time = 0
    gap_between_experiments = 1

    # Assuming there are 5 experiments and their data are stored in pairs of columns (Time [min], Int)
    for i in range(5):
        time_col_idx = 2 * i
        int_col_idx = 2 * i + 1

        # Check if the required columns exist in the DataFrame
        if time_col_idx < len(df.columns) and int_col_idx < len(df.columns):
            time_col = df.columns[time_col_idx]
            int_col = df.columns[int_col_idx]

            # Calculate the cumulative time for the current experiment
            cumulative_time_experiment = df[time_col] + cumulative_time

            # Update the cumulative time for the next experiment
            previous_cumulative_time = cumulative_time
            cumulative_time += df[time_col].max() + gap_between_experiments

            plt.plot(cumulative_time_experiment, df[int_col], label=f'Experiment {i+1}', linewidth=2)

            # Add blue lines for every second of the 1-minute gap at the end of the previous experiment
            if i > 0:
                for t in range(int(previous_cumulative_time), int(previous_cumulative_time + gap_between_experiments)):
                    plt.axvline(t, color='blue', linestyle='--', linewidth=0.7, alpha=0.3)
        else:
            print(f'Could not find required columns for Experiment {i+1}')



def plot_subplots(df, title):
    # Apply a professional style to the plot
    sns.set(style="ticks", font_scale=1.2)
    plt.rcParams['figure.figsize'] = [20, 4]

    # Create the subplot with 1 row and 5 columns
    fig, axes = plt.subplots(1, 5, sharey=True)

    # Initialize the cumulative time variable and gap between experiments
    cumulative_time = 0
    gap_between_experiments = 1

    # Assuming there are 5 experiments and their data are stored in pairs of columns (Time [min], Int)
    for i in range(5):
        time_col_idx = 2 * i
        int_col_idx = 2 * i + 1

        # Check if the required columns exist in the DataFrame
        if time_col_idx < len(df.columns) and int_col_idx < len(df.columns):
            time_col = df.columns[time_col_idx]
            int_col = df.columns[int_col_idx]

            # Calculate the cumulative time for the current experiment
            cumulative_time_experiment = df[time_col] + cumulative_time

            # Update the cumulative time for the next experiment
            cumulative_time += df[time_col].max() + gap_between_experiments

            axes[i].plot(cumulative_time_experiment, df[int_col], linewidth=2)
            axes[i].set_xlabel('Cumulative Time [min]', fontsize=12)
            axes[i].set_title(f'Experiment {i+1}', fontsize=14)

            if i == 0:
                axes[i].set_ylabel('Int', fontsize=12)

        else:
            print(f'Could not find required columns for Experiment {i+1}')

    # Adjust the space between subplots
    plt.subplots_adjust(wspace=0.3)

    # Overall title
    fig.suptitle(title, fontsize=16, y=1.05)

    # Remove the top and right spines for a cleaner look
    sns.despine()

    # Display plot
    plt.show()

# Plot the first set of experiments
plot_experiments(df1, f'Cumulative Time vs Int for 5 Experiments with 1-min gaps - {sheet1}')
plot_subplots(df1, f'Cumulative Time vs Int for 5 Experiments - {sheet1}')

# # Plot the second set of experiments
# plot_experiments(df2, f'Cumulative Time vs Int for 5 Experiments with 1-min gaps - {sheet2}')
# plot_subplots(df2, f'Cumulative Time vs Int for 5 Experiments - {sheet2}')

# #################


# import pandas as pd
# import matplotlib.pyplot as plt
# import seaborn as sns

# # Replace 'your_excel_file.xlsx' with the actual name of your Excel file
# excel_file = 'MB-RC 1 min exposures.xlsx'
# Sheet = "Python"
# Sheet2 = "Control"

# # Read the Excel file into a pandas DataFrame
# df = pd.read_excel(excel_file, sheet_name=Sheet, engine='openpyxl')
# df2 = pd.read_excel(excel_file, sheet_name=Sheet2, engine='openpyxl')

# # Initialize the cumulative time variable and gap between experiments
# cumulative_time = 0
# gap_between_experiments = 1

# # Apply a professional style to the plot
# sns.set(style="ticks", font_scale=1.2)
# plt.rcParams['figure.figsize'] = [20, 4]

# # Create the subplot with 1 row and 5 columns
# fig, axes = plt.subplots(1, 5, sharey=True)

# for i in range(5):
#     time_col_idx = 2 * i
#     int_col_idx = 2 * i + 1

#     # Check if the required columns exist in the DataFrame
#     if time_col_idx < len(df.columns) and int_col_idx < len(df.columns):
#         time_col1 = df.columns[time_col_idx]
#         int_col1 = df.columns[int_col_idx]
#         time_col2 = df2.columns[time_col_idx]
#         int_col2 = df2.columns[int_col_idx]

#         # Calculate the cumulative time for the current experiment
#         cumulative_time_experiment1 = df[time_col1] + cumulative_time
#         cumulative_time_experiment2 = df2[time_col2] + cumulative_time

#         # Update the cumulative time for the next experiment
#         cumulative_time += max(df[time_col1].max(), df2[time_col2].max()) + gap_between_experiments

#         axes[i].plot(cumulative_time_experiment1, df[int_col1], linewidth=2, label=f'Sheet1 - Experiment {i+1}')
#         axes[i].plot(cumulative_time_experiment2, df2[int_col2], linewidth=2, linestyle='--', label=f'Sheet2 - Experiment {i+1}')
#         axes[i].set_xlabel('Cumulative Time [min]', fontsize=12)
#         axes[i].set_title(f'Experiment {i+1}', fontsize=14)

#         if i == 0:
#             axes[i].set_ylabel('Int', fontsize=12)

#     else:
#         print(f'Could not find required columns for Experiment {i+1}')

# # Adjust the space between subplots
# plt.subplots_adjust(wspace=0.3)

# # Overall title
# fig.suptitle('Cumulative Time vs Int for 5 Experiments from Both Sheets', fontsize=16, y=1.05)

# # Add the legend
# axes[-1].legend(fontsize=12, loc='upper right')

# # Remove the top and right spines for a cleaner look
# sns.despine()

# # Display plot
# plt.show()

# # Save the plot as a high-resolution image
# fig.savefig('MB-RC_Both_Photocontrol_Subplots.svg', dpi=600, bbox_inches='tight')
