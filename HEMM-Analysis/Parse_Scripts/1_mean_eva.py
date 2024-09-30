import pandas as pd

# Load the CSV file into a pandas DataFrame
file_path = "../Parsed_Results/eva_results.csv"  # Replace with the correct path to your CSV file
df = pd.read_csv(file_path, delimiter=';')

# Group by 'lib', 'key', 'DIM' and calculate the mean of each numeric column
mean_df = df.groupby(['lib', 'key', 'DIM']).mean().reset_index()

# Round the numeric columns to a reasonable number of decimal places (optional)
mean_df = mean_df.round(10)

# Write the result to a new CSV file
outdir = "../Parsed_Results/"
output_file_path = outdir+'mean_eva.csv'
mean_df.to_csv(output_file_path, sep=';', index=False)

print(f"[EVA] Mean values CSV file saved to: {output_file_path}")

