import os
import csv
import pandas as pd

# Directory where your CSV files are located
directory = "./../Results/"  
outdir = "./../Parsed_Results/"
output_file = outdir+"seal_results.csv"

# Columns for output CSV
output_columns = [
    "lib", "key", "DIM", "REP", 
    "Generate HE context", "Generate HE keys", 
    "Encrypt", "HE Multiplication", "HE Decrypt"
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

        # Filter for rows where "lib" column is "SEAL"
        seal_rows = df[df['lib'] == 'SEAL']

        if len(seal_rows) > 0:
            # Process the first row as "SEAL (EVA)"
            first_seal_row = seal_rows.iloc[0]

            generate_he_context = (
                first_seal_row["t_prep_coeff_modulus"] +
                first_seal_row["t_prep_params"] +
                first_seal_row["t_prep_context"] +
                first_seal_row["t_prep_evaluator"] +
                first_seal_row["t_prep_encoder"]
            )
            generate_he_keys = (
                first_seal_row["t_prep_keygen"] +
                first_seal_row["t_prep_secret_key"] +
                first_seal_row["t_prep_public_key"] +
                first_seal_row["t_prep_relin_keys"] +
                first_seal_row["t_prep_gal_keys"] +
                first_seal_row["t_prep_encryptor"] +
                first_seal_row["t_prep_decryptor"]
            )
            he_encrypt = first_seal_row["t_matr1_matr2_matrix_multiply_encrypt"]
            he_multiplication = (
                first_seal_row["t_matr1_matr2_matrix_multiply_prep"] +
                first_seal_row["t_matr1_matr2_matrix_multiply"]
            )
            he_decrypt = first_seal_row["t_matr1_matr2_matrix_multiply_decrypt"]

            # Append "SEAL (EVA)" row
            output_data.append([
                "SEAL (EVA)", key, dim, rep, 
                generate_he_context, generate_he_keys, 
                he_encrypt, he_multiplication, he_decrypt
            ])

            # Process the rest of the SEAL rows
            for i, (_, row) in enumerate(seal_rows.iloc[1:].iterrows(), start=1):
                # Determine the lib type based on poly_modulus_degree
                poly_modulus_degree = row["poly_modulus_degree"]
                if poly_modulus_degree == 8192:
                    lib_type = "SEAL (8K)"
                elif poly_modulus_degree == 16384:
                    lib_type = "SEAL (16K)"
                elif poly_modulus_degree == 32768:
                    lib_type = "SEAL (32K)"
                else:
                    continue  # Skip if poly_modulus_degree is not one of the specified ones

                # Calculate values based on the formula provided
                generate_he_context = (
                    row["t_prep_coeff_modulus"] +
                    row["t_prep_params"] +
                    row["t_prep_context"] +
                    row["t_prep_evaluator"] +
                    row["t_prep_encoder"]
                )
                generate_he_keys = (
                    row["t_prep_keygen"] +
                    row["t_prep_secret_key"] +
                    row["t_prep_public_key"] +
                    row["t_prep_relin_keys"] +
                    row["t_prep_gal_keys"] +
                    row["t_prep_encryptor"] +
                    row["t_prep_decryptor"]
                )
                he_encrypt = row["t_matr1_matr2_matrix_multiply_encrypt"]
                he_multiplication = (
                    row["t_matr1_matr2_matrix_multiply_prep"] +
                    row["t_matr1_matr2_matrix_multiply"]
                )
                he_decrypt = row["t_matr1_matr2_matrix_multiply_decrypt"]

                # Append the result to output_data
                output_data.append([
                    lib_type, key, dim, rep, 
                    generate_he_context, generate_he_keys, 
                    he_encrypt, he_multiplication, he_decrypt
                ])

# Create a DataFrame from the collected data
output_df = pd.DataFrame(output_data, columns=output_columns)

# Sort by key (ascending), then by DIM (ascending), and then by REP (ascending)
output_df.sort_values(by=["key", "DIM", "REP"], ascending=[True, True, True], inplace=True)

# Write the output to a CSV file using semicolon as the separator
output_df.to_csv(output_file, index=False, sep=';')

print(f"[SEAL] Results saved to {output_file}")

