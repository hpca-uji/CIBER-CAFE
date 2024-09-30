import pandas as pd

# Load both CSV files
mean_seal_path = 'Parsed_Results/sorted_mean_seal.csv'
mean_eva_path = 'Parsed_Results/mean_eva.csv'

df_seal = pd.read_csv(mean_seal_path, sep=';')
df_eva = pd.read_csv(mean_eva_path, sep=';')

# Remove the 'REP' column from the EVA data as it's not needed for merging
df_eva = df_eva.drop(columns=['REP'])

# Ensure the EVA dataframe has the same column names as SEAL for consistent merging
df_eva.columns = df_seal.columns

# Define the mapping for the DIM values
dim_mapping = {
    1: 1,
    2: 2,
    4: 3,
    8: 4,
    16: 5,
    32: 6,
    64: 7,
    128: 8
}

# Apply the DIM mapping to both dataframes
df_seal['DIM'] = df_seal['DIM'].map(dim_mapping)
df_eva['DIM'] = df_eva['DIM'].map(dim_mapping)

# Concatenate the corresponding EVA rows before each SEAL (EVA) row
combined_rows = []
for index, row in df_seal.iterrows():
    # If the row corresponds to SEAL (EVA), insert the matching EVA row first
    if row['lib'] == 'SEAL (EVA)':
        matching_eva_row = df_eva[(df_eva['key'] == row['key']) & (df_eva['DIM'] == row['DIM'])]
        if not matching_eva_row.empty:
            combined_rows.append(matching_eva_row.iloc[0])  # Add the corresponding EVA row
    combined_rows.append(row)  # Then add the current SEAL row

# Create a new DataFrame from the combined rows
df_combined = pd.DataFrame(combined_rows)

# Rename the columns as requested
df_combined = df_combined.rename(columns={
    'lib': 'Library',
    'DIM': 'Matrix Dimension'
})

# Save the reordered DataFrame to a new CSV file
output_file_path = 'Parsed_Results/combined_mean_seal_eva.csv'
df_combined.to_csv(output_file_path, sep=';', index=False)

print("Rows combined, DIM values modified, headers updated, and saved to", output_file_path)

