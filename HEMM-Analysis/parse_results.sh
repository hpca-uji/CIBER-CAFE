echo "-----------------------------------------"
echo "[1] EXECUTING PARSE SCRIPTS"
echo "-----------------------------------------"

# Create directory to store CSV files of parsed results
mkdir Parsed_Results

# Execute parse scripts
cd Parse_Scripts
python3 0_parse_plain.py
python3 0_parse_EVA.py
python3 1_mean_eva.py
python3 0_parse_SEAL.py
python3 1_mean_seal.py
cd ..

# Generate the final CSV files
echo "-----------------------------------------"
echo "[2] GENERATING THE FINAL CSV FILES"
echo "-----------------------------------------"
python3 Parse_Scripts/2_order_rows_seal.py
python3 Parse_Scripts/2_order_rows_eva.py
echo "Library;key;Matrix Dimension;Generate HE context;Generate HE keys;HE Encrypt;HE Multiplication;HE Decrypt" > times_128.csv
cat Parsed_Results/combined_mean_seal_eva.csv | grep ";128;" >> times_128.csv
echo "Library;key;Matrix Dimension;Generate HE context;Generate HE keys;HE Encrypt;HE Multiplication;HE Decrypt" > times_256.csv
cat Parsed_Results/combined_mean_seal_eva.csv | grep ";256;" >> times_256.csv
echo "File times_128.csv and times_256.csv generated"

# Generate plots
echo "-----------------------------------------"
echo "[3] GENERATING PLOTS"
echo "-----------------------------------------"
python3 plot_128_bars.py
echo "Plot for key=128 generated"
python3 plot_256_bars.py
echo "Plot for key=256 generated"

# Remove unnecessary files
echo "-----------------------------------------"
echo "[4] REMOVING UNNECESSARY FILES"
echo "-----------------------------------------"
echo "rm -rf Parsed_Results"
rm -rf Parsed_Results
echo "-----------------------------------------"
echo "[5] DONE"
echo "-----------------------------------------"
