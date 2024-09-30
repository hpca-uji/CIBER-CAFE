import os
import pandas as pd

# Define the base directory where the files are stored
base_dir = "../Results_plain"  
outdir = "../Parsed_Results/"  
output_file = outdir+"plain_results.csv"  # Name of the output file

# Define the range of X and rep values
X_values = [1, 2, 4, 8, 16, 32, 64, 128]
rep_values = [1, 2, 3, 4, 5]

# List to store all dataframes
df_list = []

# Loop through each combination of X and rep
for X in X_values:
    for rep in rep_values:
        # Construct the filename
        file_name = f"res_P{X}_{rep}.csv"
        file_path = os.path.join(base_dir, file_name)

        # Check if the file exists before reading
        if os.path.exists(file_path):
            # Read the CSV into a DataFrame, selecting only the required column
            df = pd.read_csv(file_path, sep=';', usecols=['t_matr1_matr2_matrix_multiply_total'])

            # Add metadata columns for X and rep
            df['P_value'] = X
            df['rep'] = rep

            # Append to the list of dataframes
            df_list.append(df)
        else:
            print(f"File {file_path} does not exist!")

# Combine all DataFrames into one
combined_df = pd.concat(df_list, ignore_index=True)

# Compute the mean of the selected column for each X (P_value), ignoring the rep column
mean_df = combined_df.groupby('P_value')['t_matr1_matr2_matrix_multiply_total'].mean().reset_index()

# Rename the columns to 'dimension' and 'Time'
mean_df.columns = ['Dimension', 'Time']

# Save the resulting DataFrame with the mean values to a new CSV file
mean_df.to_csv(output_file, sep=';', index=False)

print(f"[PLAIN] Mean combined CSV saved as {output_file}")
