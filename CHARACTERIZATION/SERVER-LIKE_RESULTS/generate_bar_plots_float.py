import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Color mapping for libraries
color_map = {
    'HElib': 'orange',
    'TFHE': 'blue',
    'palisade': 'red',
    'SEAL': 'green'
}

dims = ['4', '8', '16', '32']
for dimension in dims:

    # Read the CSV file into a Pandas DataFrame
    df = pd.read_csv('dimension_' + str(dimension) + '_times.csv', delimiter=';')
    
    # Extract the 'Library' column for labels
    libraries = df['Library']

    # Extract Orin and Xavier columns for int and float
    float_columns = [col for col in df.columns if col.endswith('-float')]

    float_data = df[float_columns]

    # Remove TFHE from float data
    if 'TFHE' in libraries.values:
        float_data = float_data[libraries != 'TFHE']
        libraries_filtered_float = libraries[libraries != 'TFHE']
    else:
        libraries_filtered_float = libraries

    # Function to create bar plots
    def create_bar_plot(data, title, filename, filtered_libraries):
        width = 0.2  # the width of the bars

        fig, ax = plt.subplots(figsize=(10, 6))

        # Ensure that the data rows correspond to the filtered libraries
        filtered_data = data.loc[libraries.isin(filtered_libraries)].reset_index(drop=True)

        # Check if filtered_data has the right number of rows
        num_bars = len(filtered_libraries)
        if filtered_data.shape[0] != num_bars:
            raise ValueError("Filtered data does not match the number of libraries.")

        x = np.arange(len(x_ticks))  # the label locations

        # Plot each library
        for i, library in enumerate(filtered_libraries):
            data_for_library = filtered_data.iloc[i].values
            color = color_map.get(library, 'gray')  # Default to gray if library is not in color_map
            ax.bar(x + i * width, data_for_library, width, label=library, color=color)

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_xlabel('System')
        ax.set_ylabel('Time (ms)')
        ax.set_title(title)
        ax.set_xticks(x + width * (num_bars - 1) / 2)
        ax.set_xticklabels(x_ticks)

        # Adjusting legend position
        ax.legend(loc='upper center', ncol=4)

        # Save the plot to a PDF file
        plt.savefig(filename, format='png')
        plt.close()

    # Common X tick labels for both plots
    x_ticks = ['Orin', 'Xavier', 'Nano']
    
    # Plotting the Orin-* and Xavier-* float values and saving to PDF
    create_bar_plot(float_data, 'Matrix multiplication, dim ' + str(dimension) + 'x' + str(dimension) + ' (CKKS, float)', str(dimension) + 'x' + str(dimension) + '_float_times_barplot.png', libraries_filtered_float.tolist())

    print(f"Bar plots for dimension {dimension}x{dimension} have been saved to 2 PDF files.")

