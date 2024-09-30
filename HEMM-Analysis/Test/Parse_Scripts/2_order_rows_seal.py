import pandas as pd

# Load the CSV file
file_path = 'Parsed_Results/mean_seal.csv'
df = pd.read_csv(file_path, sep=';')

# Define the specific order for 'lib'
lib_order = ['SEAL (EVA)', 'SEAL (8K)', 'SEAL (16K)', 'SEAL (32K)']

# Convert 'lib' to a categorical type with the defined order
df['lib'] = pd.Categorical(df['lib'], categories=lib_order, ordered=True)

# Sort the DataFrame first by 'key', then by 'DIM', and finally by 'lib' using the specified order
df_sorted = df.sort_values(by=['key', 'DIM', 'lib'])

# Save the sorted DataFrame to a new CSV file
output_file_path = 'Parsed_Results/sorted_mean_seal.csv'
df_sorted.to_csv(output_file_path, sep=';', index=False)

print("Rows reordered and saved to", output_file_path)

