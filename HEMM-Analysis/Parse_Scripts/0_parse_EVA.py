import os
import csv
import pandas as pd

# Directory where your CSV files are located
directory = "../Results/"  
outdir = "../Parsed_Results/"
output_file = outdir+"eva_results.csv"

# Columns for output CSV
output_columns = [
    "lib", "key", "DIM", "REP", 
    "Generate HE context", "Generate HE keys", 
    "HE Encrypt", "HE Multiplication", "HE Decrypt"
]

# Create an empty dataframe to store results
output_data = []

# Process each CSV file in the directory
for filename in os.listdir(directory):
    if filename.endswith(".csv"):
        # Extract KEY, DIM, REP from filename
        parts = filename.split('_')
        key = int(parts[1])
        dim = int(parts[3])
        rep = int(parts[4].split('.')[0])

        # Read the CSV file
        file_path = os.path.join(directory, filename)
        df = pd.read_csv(file_path, sep=';')

        # Filter for rows where "lib" column is "EVA"
        eva_rows = df[df['lib'] == 'EVA']

        for _, row in eva_rows.iterrows():
            # Extract required fields
            generate_he_context = row["t_program"] + row["t_compiler"] + row["t_compile"]
            generate_he_keys = row["t_generate_keys"]
            he_encrypt = row["t_encrypt"]
            he_multiplication = row["t_execute"]
            he_decrypt = row["t_decrypt"]

            # Append the result to output_data
            output_data.append([
                row["lib"], key, dim, rep, 
                generate_he_context, generate_he_keys, 
                he_encrypt, he_multiplication, he_decrypt
            ])

# Create a DataFrame from the collected data
output_df = pd.DataFrame(output_data, columns=output_columns)

# Sort by key (ascending), then by DIM (ascending), and then by REP (ascending)
output_df.sort_values(by=["key", "DIM", "REP"], ascending=[True, True, True], inplace=True)

# Write the output to a CSV file using semicolon as the separator
output_df.to_csv(output_file, index=False, sep=';')

print(f"[EVA] Results saved to {output_file}")

