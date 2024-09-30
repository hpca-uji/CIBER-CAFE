import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Read the CSV data from a file
csv_file_path = 'time-256_noPlain_newLabel.csv'  # Replace with the path to your CSV file
data = pd.read_csv(csv_file_path, sep=';')

# Convert relevant columns to numeric, coercing errors to NaN, then fill NaN with 0
for column in ["Generate HE context", "Generate HE keys", "HE Encrypt", "HE Multiplication", "HE Decrypt"]:
    data[column] = pd.to_numeric(data[column].str.replace(',', '.'), errors='coerce').fillna(0)

# Define colors for the bars
colors = {
    "Generate HE context": "#1f77b4",
    "Generate HE keys": "#ff7f0e",
    "HE Encrypt": "#2ca02c",
    "HE Multiplication": "#d62728",
    "HE Decrypt": "#9467bd"
}

# Prepare the plot
plt.figure(figsize=(16, 8))
bar_width = 0.35

# Determine the positions of the bars
unique_dimensions = sorted(data['Matrix Dimension'].unique())
num_dimensions = len(unique_dimensions)
positions = []
current_position = 0

# Generate the positions for each group of bars
for dimension in unique_dimensions:
    sub_data = data[data['Matrix Dimension'] == dimension]
    num_subgroups = len(sub_data)
    pos = np.arange(num_subgroups) + current_position
    positions.extend(pos)
    current_position += num_subgroups + 1  # Add space between groups

# Plot each component in the stacked bar
bottom = np.zeros(len(data))
for column in ["Generate HE context", "Generate HE keys", "HE Encrypt", "HE Multiplication", "HE Decrypt"]:
    plt.bar(positions, data[column], bar_width, label=column, bottom=bottom, color=colors[column])
    bottom += data[column]
#    print(data[column], " ")

# Log scale for y-axis
plt.yscale('log')

# Set labels and title
plt.xlabel('Library and matrix dimension', fontsize=18)
plt.ylabel('Execution time (ms)', fontsize=18)
plt.title('Execution time break down (key=256)', pad=26, fontsize=24)

# Create custom x-tick labels without the "Dim X" part
xtick_labels = []
for dimension in unique_dimensions:
    sub_data = data[data['Matrix Dimension'] == dimension]
    labels = sub_data['Library'].tolist()
    xtick_labels.extend(labels)

# Set x-ticks and labels
plt.xticks(positions, xtick_labels, rotation=90, fontsize=17)
plt.yticks(fontsize=17)
#plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
plt.legend(loc='upper center', fontsize=14, ncol=5, columnspacing=2)

# Add vertical lines to separate groups visually
offset = 0.5  # Adjust this value as needed to move the line one position to the right
for i, dimension in enumerate(unique_dimensions[:-1]):
    num_subgroups = len(data[data['Matrix Dimension'] == dimension])
    line_position = current_position - (num_subgroups + 1) * (num_dimensions - i - 1) - 0.5
    plt.axvline(line_position + offset, color='black', linewidth=0.8)


# Add group labels at the top of the plot
group_labels = ["1x1", "2x2", "4x4", "8x8", "16x16", "32x32", "64x64", "128x128"]
group_positions = []
for i, dimension in enumerate(unique_dimensions):
    num_subgroups = len(data[data['Matrix Dimension'] == dimension])
    line_position = current_position - (num_subgroups + 1) * (num_dimensions - i - 1) - 0.5
    mean_pos = line_position - 2
    group_positions.append(mean_pos)

plt.ylim(bottom=0.01, top=1000)
for pos, label in zip(group_positions, group_labels):
    plt.text(pos, plt.ylim()[1] * 1.05, label, ha='center', va='bottom', fontsize=17, color='black')

# Show plot
plt.tight_layout()
plt.savefig('Time256.png')
plt.show()
