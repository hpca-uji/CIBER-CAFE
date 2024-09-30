import pandas as pd

# Load the CSV file
# Adjust the path to where your CSV file is located
file_path = '0_PARSED-Results/seal_results.csv'
df = pd.read_csv(file_path, sep=';', header=0)

# Define the columns for computation
value_columns = ['Generate HE context', 'Generate HE keys', 'Encrypt', 'HE Multiplication', 'HE Decrypt']

# Group by 'lib', 'key', and 'DIM', and compute the mean for the specified columns
mean_values = df.groupby(['lib', 'key', 'DIM'])[value_columns].mean().reset_index()

# Optional: If you want to save the results to a new CSV file
output_file_path = '1_MEAN-Results/mean_seal.csv'
mean_values.to_csv(output_file_path, sep=';', index=False)

print("Mean values computed and saved to", output_file_path)

