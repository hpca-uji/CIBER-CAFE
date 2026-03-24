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

    int_columns = [col for col in df.columns if col.endswith('-int')]

    int_data = df[int_columns]

    # Common X tick labels for both plots
    x_ticks = ['Orin', 'Xavier', 'Nano']
    x = np.arange(len(x_ticks))  # the label locations

    # Function to create bar plots
    def create_bar_plot(data, title, filename):
        width = 0.2  # the width of the bars

        fig, ax = plt.subplots(figsize=(10, 6))

        # Loop through each library to plot its data
        for i, library in enumerate(libraries):
            color = color_map.get(library, 'gray')  # Default to gray if library is not in color_map
            ax.bar(x + i * width, data.iloc[i], width, label=library, color=color)

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_xlabel('System')
        ax.set_ylabel('Time (ms)')
        ax.set_title(title)
        ax.set_xticks(x + width * (len(libraries) - 1) / 2)
        ax.set_xticklabels(x_ticks)

        # Adjusting legend position
        ax.legend(loc='upper center', ncol=4)

        # Save the plot to a PDF file
        plt.savefig(filename, format='png')
        plt.close()

    # Plotting the Orin-* and Xavier-* int values and saving to PDF
    create_bar_plot(int_data, 'Matrix multiplication, dim ' + str(dimension) + 'x' + str(dimension) + ' (BFV, int)', str(dimension) + 'x' + str(dimension) + '_int_times_barplot.png')


    print(f"Bar plots for dimension {dimension}x{dimension} have been saved to 2 PDF files.")

