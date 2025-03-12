# codigo que analiza la agrupación óptima respecto a los macro-grupos resultantes
# muestra un heatmap con la frecuencia que aparece cada combinación de grupo en la solución final

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import json
import numpy.ma as ma

# Print the loaded JSON data
verbose = True # print outputs

election = "2021_11_Presidencial"

# función que pasa dos indices a la notación de macro-grupos
def macro_group_notation(index_from, index_to):
    # index_from : index with values from 0 to 7
    # index_to   : index with values from 0 to 7
    #age_ranges = ["18-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"]
    age_ranges_from = ["18", "20", "30", "40", "50", "60", "70", "80"]
    age_ranges_to = ["-19", "-29", "-39", "-49", "-59", "-69", "-79", "+"]
    return f"{age_ranges_from[index_from]}{age_ranges_to[index_to]}"


def heatmap_agg_groups_test():
    
    # dict_freq = {
    #     "18-19_HM": 0.1,
    #     "20-29_HM": 0.2,
    #     "30-39_HM": 0.3,
    #     "40-49_HM": 0.4,
    #     "50-59_HM": 0.5,
    #     "60-69_HM": 0.6,
    #     "70-79_HM": 0.7,
    #     "80+_HM": 0.8,
    #     "18-29_HM": 0.1,
    #     "20-39_HM": 0.2,
    #     "30-49_HM": 0.3,
    #     "40-59_HM": 0.4,
    #     "50-69_HM": 0.5,
    #     "60-79_HM": 0.6,
    #     "70+_HM": 0.7,
    #     "18-39_HM": 0.1,
    #     "20-49_HM": 0.2,
    #     "30-59_HM": 0.3,
    #     "40-69_HM": 0.4,
    #     "50-79_HM": 0.5,
    #     "60+_HM": 0.6,
    #     "18-49_HM": 0.1,
    #     "20-59_HM": 0.2,
    #     "30-69_HM": 0.3,
    #     "40-79_HM": 0.4,
    #     "50+_HM": 0.5,
    #     "18-59_HM": 0.1,
    #     "20-69_HM": 0.2,
    #     "30-79_HM": 0.3,
    #     "40+_HM": 0.4,
    #     "18-69_HM": 0.1,
    #     "20-79_HM": 0.2,
    #     "30+_HM": 0.3,
    #     "18-79_HM": 0.1,
    #     "20+_HM": 0.2,
    #     "18+_HM": 0.1,
    # }

    age_ranges = ["18-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+"]

    # read JSON file group_count_dict
    with open(f'{election}/output_L100/group_count_dict.json', 'r') as file:
        # Load the JSON data
        dict_freq = json.load(file)
    file.close()
    # complete with zeros
    for i in range(len(age_ranges)):
        for j in range(i, len(age_ranges)):
            if macro_group_notation(i, j) not in dict_freq:
                dict_freq[macro_group_notation(i, j)] = 0
    # normalize
    for key in dict_freq:
        dict_freq[key] = 100 * dict_freq[key] / 679
        
    data = np.nan * np.zeros((8, 20))
    
    data_column_start = 0
    rectangles = []
    texts = []
    columns_groups_sizes = [] # tamaño de los macro-grupos
    for i in range(8): # recorre los tamaños de los macro-grupos
        n_columns = i + 1 if i < 4 else 8 - i
        for _ in range(n_columns):
            columns_groups_sizes.append(i + 1)
        from_to_tuples = [(j, j + i) for j in range(8) if j + i < 8]
        #print(f"i = {i}; n_columns = {n_columns}")
        #print(f"from_to_tuples = {from_to_tuples}")
        #print(f"data_column_start = {data_column_start}")
        for from_to in from_to_tuples: # recorre las combinaciones de macro-grupos de tamaño i + 1
            #print(f"\tfrom_to = {from_to}")
            index_from, index_to = from_to
            for j in range(index_from, index_to + 1): # recorre las celdas de un rectángulo asociado a un macro-grupo
                #print(f"\t\tj = {j}; data_column_start + index_from % n_columns = {data_column_start + index_from % n_columns}")
                heatmap_column = data_column_start + index_from % n_columns
                data[j, heatmap_column] = dict_freq[macro_group_notation(index_from, index_to)]
            rectangles.append((index_from, heatmap_column, 1, i + 1)) # (row, column, width, height)
            texts.append((heatmap_column, index_from + (i + 0) / 2, f'{dict_freq[macro_group_notation(index_from, index_to)]:.1f}', 'black' if dict_freq[macro_group_notation(index_from, index_to)] < 12.0 else 'white'))
        data_column_start += n_columns        

    #masked_data = ma.masked_invalid(data)
    #print(f"masked_data = {masked_data}")
    # Create a custom colormap from white to blue
    cmap_colors = [(1, 1, 1), (0, 0, 1)]  # White to blue
    cmap = LinearSegmentedColormap.from_list('custom_cmap', cmap_colors)

    # Create the heatmap
    #print(f"data.nanmax = {np.nanmax(data)}")
    # plt.imshow(data, cmap=cmap, interpolation='nearest', vmin=0, vmax=np.nanmax(data))
    plt.figure(figsize=(10.5, 4.5))  # Change width and height as needed
    #nan_color = 'gray'
    #data_with_nan_color = np.where(np.isnan(data), nan_color, data)
    #plt.imshow(data, cmap=cmap, interpolation='nearest')
    reversed_map = plt.cm.get_cmap('viridis').reversed() # set_bad(color='gray')
    extent = [-0.5, data.shape[1]-0.5, data.shape[0]-0.5, -0.5]
    plt.imshow(data, cmap=reversed_map, interpolation='nearest', extent=extent)
    for text in texts:
        plt.text(text[0], text[1], text[2], ha='center', va='center', color=text[3])
        #plt.text(j, i, f'{data[i, j]:.2f}', ha='center', va='center', color='white')
    #plt.imshow(data_with_nan_color, cmap='viridis', interpolation='nearest')
    #plt.imshow(masked_data, cmap='gray', interpolation='nearest', alpha=0.5)  # Set alpha for transparency
    for cell in rectangles:
        plt.gca().add_patch(plt.Rectangle((cell[1]-0.5, cell[0]-0.5), cell[2], cell[3], fill=False, edgecolor='black', linewidth=1.5))
    plt.colorbar(shrink=0.7, format='%.0f%%', aspect=10)
    plt.yticks(np.arange(len(age_ranges)), age_ranges)
    plt.xticks(np.arange(len(columns_groups_sizes)), columns_groups_sizes)
    # x-axis label
    plt.xlabel('Number of groups in each macro-group')
    plt.ylabel('Age ranges')
    plt.title('Percentage of districts with each macro-group\ncombination in the optimal solution')
    plt.xlim(extent[:2])  # Limit x-axis to extent
    plt.ylim(extent[2:])  # Limit y-axis to extent
    plt.tight_layout()  # Adjust layout to reduce empty space
    
    plt.savefig("images/2021_11_Presidencial_MACRO_GROUPS_HEAT.pdf")
    plt.show()

    
if __name__ == '__main__':
    heatmap_agg_groups_test()
    

# /opt/homebrew/bin/python3.9 heatmap_agg_groups.py
