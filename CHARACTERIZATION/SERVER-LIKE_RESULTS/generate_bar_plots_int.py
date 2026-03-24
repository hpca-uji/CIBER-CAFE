import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Dimensions and configurations
dims = ['4', '8', '16', '32']
configs = ['c00', 'c01', 'c02', 'c03']

# Color mapping for libraries
color_map = {
    'HElib': 'orange',
    'TFHE': 'blue',
    'palisade': 'red',
    'SEAL': 'green'
}

for dimension in dims:
    # Read the CSV file into a Pandas DataFrame
    df = pd.read_csv(f'dimension_{dimension}_times.csv', delimiter=';')

    # Extract the 'Library' column for labels
    libraries = df['Library']
    
    # Extract int columns
    int_columns = [col for col in df.columns if col.startswith('int-')]
    int_data = df[int_columns]

    # Common X tick labels for plots
    x_ticks = configs
    x = np.arange(len(x_ticks))  # the label locations

    # Function to create bar plots
    def create_bar_plot(data, title, filename):
        width = 0.2  # the width of the bars

        fig, ax = plt.subplots(figsize=(12, 8))

        # Initialize a DataFrame for plotting to handle missing values
        plot_data = pd.DataFrame(index=libraries, columns=x_ticks)

        # Fill in the DataFrame with the data
        for library in libraries:
            library_data = data.loc[df['Library'] == library].values.flatten()
            if len(library_data) == len(x_ticks):
                plot_data.loc[library] = library_data
            else:
                # Debug: Print data length for each library
                print(f"Warning: Data length for library {library} does not match x_ticks length. Filling missing values with NaN.")
                # Fill missing values with NaN
                if len(library_data) < len(x_ticks):
                    library_data = np.concatenate([library_data, [np.nan] * (len(x_ticks) - len(library_data))])
                plot_data.loc[library] = library_data

        # Debug: Print out the plot_data DataFrame
        print(f"Plot Data for dimension {dimension}x{dimension}:\n{plot_data}")

        # Plot each library
        for i, library in enumerate(plot_data.index):
            row_data = plot_data.loc[library].values
            # Remove NaNs for plotting
            row_data = np.nan_to_num(row_data, nan=np.nan)
            color = color_map.get(library, 'gray')  # Default to gray if library is not in color_map
            ax.bar(x + i * width, row_data, width, label=library, color=color)

        # Add some text for labels, title and custom x-axis tick labels, etc.
        ax.set_xlabel('Server')
        ax.set_ylabel('Time (ms)')
        ax.set_title(title)
        ax.set_xticks(x + width * (len(plot_data.index) - 1) / 2)
        ax.set_xticklabels(x_ticks)

        # Adjusting legend position
        ax.legend(loc='upper center', ncol=4)

        # Save the plot to a PNG file
        plt.savefig(filename, format='png')
        plt.close()

    # Plotting the c00-* to c03-* int values and saving to PNG
    create_bar_plot(int_data, f'Matrix multiplication, dim {dimension}x{dimension} (BFV, int)', f'{dimension}x{dimension}_int_times_barplot.png')

    print(f"Bar plot for dimension {dimension}x{dimension} has been saved to PNG file.")
