import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# Input data as a dictionary
data = {
    "Sample": [
        "Thermal Equilibrium",
        "UV (First Irradiation)",
        "Blue (First Irradiation)",
        "UV (Second Irradiation)",
        "Blue (Second Irradiation)",
    ],
    "RC": [545, 1068, 899, 1068, 887],
    "M8": [546, 568, 664, 1051, 863],
    "M9": [431, 839, 763, 941, 756],
    "M11": [399, 834, 763, 898, 754],
}

# Create a DataFrame from the data
df = pd.DataFrame(data)

# Calculate the relative change for each column
df_relative_change = df.iloc[:, 1:].div(df.iloc[0, 1:], axis=1)

# Set the Sample column as the index for the DataFrame
df_relative_change["Sample"] = df["Sample"]
df_relative_change.set_index("Sample", inplace=True)

# Exclude the Thermal Equilibrium data (remove first row)
df_relative_change = df_relative_change.iloc[1:]

# Customize the appearance of the graph
plt.style.use("default")  # Reset the style to the default one

# Set the font to Helvetica Neue
plt.rcParams["font.family"] = "Helvetica Neue"

# Define the colors for the gradient (bright fushia and royal navy blue)
color_gradient = LinearSegmentedColormap.from_list("custom_gradient", ["#FF00FF", "#002147"], N=256)

# Plot the data
ax = df_relative_change.T.plot(kind="bar", figsize=(12, 6), colormap=color_gradient, edgecolor="black", linewidth=1)

plt.title("Relative Fluorescence Intensities")
plt.ylabel("Relative Change")
plt.xticks(rotation=0)
plt.legend(title="Sample", bbox_to_anchor=(1.05, 1), loc="upper left")

# Make the y-axis start at 0
plt.ylim(0, max(df_relative_change.max())+0.2)

# Remove background color
ax.set_facecolor("none")

# Remove the right and top lines of the graph
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)

# Add gridlines
ax.yaxis.grid(True, linestyle="--", linewidth=1)

plt.tight_layout()
plt.show()
