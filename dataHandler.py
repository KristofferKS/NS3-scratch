import pandas as pd
import pickle
from collections import Counter  # Add this import at the top if not already present
import matplotlib.pyplot as plt

nNodes = 50
connections = 258
foldername = "50-simidk"


nodeData = []

for i in range(nNodes):
    df = pd.read_csv(f"simulations/{foldername}/gossip_log_node{i}.csv")
    nodeData.append(df)

totalData = pd.read_csv(f"simulations/{foldername}/final_knowledge.csv")


combined_df = pd.concat(nodeData, ignore_index=True)

# Precompute gen info (timestamp_ms and message_id) and group updates by message_id
gen_df = combined_df[combined_df['event'] == 'gen']
gen_list = gen_df[['timestamp_ms', 'node_from']].to_dict('records')  # List of dicts: [{'timestamp_ms': val, 'message_id': val}, ...]
updates_df = combined_df[combined_df['event'] == 'update']
updates_grouped = updates_df.groupby('node_from')

def get_path(gen_number):
    if gen_number % 100 == 0:
        percent = (gen_number + 1) / genEvents * 100
        print(f"Processing gen event {gen_number + 1} / {genEvents} ({percent:.1f}%)")
    gen_info = gen_list[gen_number]
    message_id = gen_info['node_from']
    if message_id in updates_grouped.groups:
        stepByStep_df = updates_grouped.get_group(message_id)
        # Filter to only updates with the same timestamp_ms as the gen event
        stepByStep_df = stepByStep_df[stepByStep_df['timestamp_ms'] == gen_info['timestamp_ms']]
        # Sort by 'time' (arrival time) to get the propagation order
        stepByStep_df = stepByStep_df.sort_values(by='time')
    else:
        stepByStep_df = pd.DataFrame()  # No updates for this message
    return stepByStep_df

def print_all_paths(paths):
    for i in range(len(paths)):
        print(f"Path for gen event {i}:")
        print(paths[i])

def get_path_lengths(paths):
    lengths = {}
    for path in paths:
        if not path.empty:
            send_to_value = int(path['node_to'].iloc[0])
            if send_to_value not in lengths:
                lengths[send_to_value] = []
            lengths[send_to_value].append(len(path)+1)
    return lengths

genEvents = len(gen_list)  # Update to use gen_list length
def get_paths_lengths():
    paths = [get_path(i) for i in range(genEvents)]
    return paths

# only if paths.pkl does not exist
import os
updates_grouped = updates_df.groupby(['node_from', 'timestamp_ms'])

from multiprocessing import Pool

def get_path_fast(gen_number):
    gen_info = gen_list[gen_number]
    message_id = gen_info['node_from']
    timestamp_ms = gen_info['timestamp_ms']
    
    try:
        stepByStep_df = updates_grouped.get_group((message_id, timestamp_ms))
        stepByStep_df = stepByStep_df.sort_values(by='time', ignore_index=True)
        return stepByStep_df
    except KeyError:
        return pd.DataFrame()

def get_paths_lengths_parallel():
    num_cores = os.cpu_count()
    print(f"Using {num_cores} cores for parallel processing")
    total = len(gen_list)
    with Pool(num_cores) as pool:
        paths = []
        for i, path in enumerate(pool.imap(get_path_fast, range(total)), 1):
            if i % 1000 == 0:
                print(f"Processed {i}/{total} ({i/total*100:.1f}%)")
            paths.append(path)
    return paths

genEvents = len(gen_list)

if not os.path.exists(f"simulations/{foldername}/paths.pkl"):
    paths = get_paths_lengths_parallel()
    with open(f"simulations/{foldername}/paths.pkl", "wb") as f:
        pickle.dump(paths, f)
else:
    with open(f"simulations/{foldername}/paths.pkl", "rb") as f:
        paths = pickle.load(f)
def get_average_path_length(lengths_dict):
    averages = {}
    for key, lengths_list in lengths_dict.items():
        if lengths_list:
            averages[key] = sum(lengths_list) / len(lengths_list)
        else:
            averages[key] = 0  # Or handle empty lists as needed
    return averages

def get_procentage_covered(average_lengths, total_nodes=nNodes):
    procentages = {}
    for key, avg_length in average_lengths.items():
        procentages[key] = (avg_length / total_nodes) * 100
    return procentages

# Example usage: Compute averages for the lengths of the first 10 paths
# lengths_for_first_10 = get_path_lengths(paths[0:10])
# average_lengths = get_average_path_length(lengths_for_first_10)
# print(f"Average path lengths per node_to for first 10: {average_lengths}")

# Or for all paths
lengths_all = get_path_lengths(paths)
average_lengths_all = get_average_path_length(lengths_all)
sorted_averages = dict(sorted(average_lengths_all.items()))
#print(f"Average path lengths per node_to overall: {sorted_averages}")

procentage_covered_all = get_procentage_covered(average_lengths_all)
sorted_procentages = dict(sorted(procentage_covered_all.items()))
# print(f"Procentage covered per node_to overall: {sorted_procentages}")

# Combine all lengths and count frequencies
def get_length_frequencies(lengths_dict):
    all_lengths = []
    for lengths_list in lengths_dict.values():
        all_lengths.extend(lengths_list)
    return dict(Counter(all_lengths))

def get_procentage_for_succesfull_jump(length_frequencies):
    total_paths = sum(length_frequencies.values())
    cumulative = 0
    succesfull_jumps = {'1': 100, '2': 100}  # 100% for length 1 and 2
    for length, freq in length_frequencies.items():
        succesfull_jumps[str(length - 1)] = ((total_paths - cumulative) / total_paths) * 100
        cumulative += freq
    return succesfull_jumps

def plot_frequency_distribution(length_frequencies):

    lengths = list(length_frequencies.keys())
    frequencies = list(length_frequencies.values())

    plt.bar(lengths, frequencies)
    plt.xlabel('Path Length')
    plt.ylabel('Frequency')
    plt.title(f'Path Length Frequency Distribution for {nNodes} Nodes')
    plt.savefig(f'simulations/{foldername}/pathLengthFreqDist{nNodes}{connections}.png')
    #plt.show()



lengths_all = get_path_lengths(paths)
frequencies = get_length_frequencies(lengths_all)
sorted_frequencies = dict(sorted(frequencies.items()))
procentage_jumps = get_procentage_for_succesfull_jump(sorted_frequencies)

def save_csv_files():
    # CSV 1: Average Path Lengths
    avg_df = pd.DataFrame(
        [(int(k), v) for k, v in sorted_averages.items()],
        columns=['Node', 'Average Path Length (Jumps)']
    )
    avg_df.to_csv(f'simulations/{foldername}/avg_path_lengths.csv', index=False)
    
    # CSV 2: Procentage Covered
    perc_cov_df = pd.DataFrame(
        [(int(k), v) for k, v in sorted_procentages.items()],
        columns=['Node', 'Average transmission Coverage (%)']
    )
    perc_cov_df.to_csv(f'simulations/{foldername}/avg_transmission_coverage.csv', index=False)
    
    # CSV 3: Length Frequencies
    freq_df = pd.DataFrame(
        [[1, 0], [2, 0]] + [(int(k), v) for k, v in sorted_frequencies.items()],
        columns=['Jumps before transmission Death', 'quantity']
    )
    freq_df.to_csv(f'simulations/{foldername}/point_of_transmission_death.csv', index=False)
    
    # CSV 4: Procentage for Successful Jump
    jump_df = pd.DataFrame(
        [(int(k), v) for k, v in procentage_jumps.items()],
        columns=['Path Length', 'Probability for Successful Jump (%)']
    )
    jump_df.to_csv(f'simulations/{foldername}/probability_successful_jump.csv', index=False)
    
    print(f"CSV files saved to simulations/{foldername}/")

def save_xlsx_file():
    with pd.ExcelWriter(f'simulations/{foldername}/path_analysis_{nNodes}_{connections}.xlsx') as writer:
        # Sheet 1: Average Path Lengths (keys as integers)
        avg_df = pd.DataFrame(
            [(int(k), v) for k, v in sorted_averages.items()],
            columns=['Node', 'Average Path Length (Jumps)']
        )
        avg_df.to_excel(writer, sheet_name='Average Path Lengths', index=False)
        
        # Sheet 2: Procentage Covered (keys as integers)
        perc_cov_df = pd.DataFrame(
            [(int(k), v) for k, v in sorted_procentages.items()],
            columns=['Node', 'Average transmission Coverage (%)']
        )
        perc_cov_df.to_excel(writer, sheet_name='Average Transmission Coverage', index=False)
        
        # Sheet 3: Length Frequencies (keys as integers)
        freq_df = pd.DataFrame(
            [[1, 0], [2, 0]] + [(int(k), v) for k, v in sorted_frequencies.items()],
            columns=['Jumps before transmission Death', 'quantity']
        )
        freq_df.to_excel(writer, sheet_name='Point of Transmission Death', index=False)
        
        # Sheet 4: Procentage for Successful Jump (keys as integers)
        jump_df = pd.DataFrame(
            [(int(k), v) for k, v in procentage_jumps.items()],
            columns=['Path Length', 'Density ']
        )
        jump_df.to_excel(writer, sheet_name='Probability for Successful Jump', index=False)

def save_xlsx_file_with_graphs():
    with pd.ExcelWriter(f'simulations/{foldername}/path_analysis_{nNodes}_{connections}_withGraphs.xlsx', engine='xlsxwriter') as writer:
        # Sheet 1: Average Path Lengths
        avg_df = pd.DataFrame(
            [(int(k), v) for k, v in sorted_averages.items()],
            columns=['Node', 'Average Path Length (Jumps)']
        )
        avg_df.to_excel(writer, sheet_name='Average Path Lengths', index=False)
        
        # Add chart for Sheet 1
        workbook = writer.book
        worksheet = writer.sheets['Average Path Lengths']
        chart1 = workbook.add_chart({'type': 'column'})
        chart1.add_series({
            'name': 'Average Path Length',
            'categories': ['Average Path Lengths', 1, 0, len(avg_df), 0],
            'values': ['Average Path Lengths', 1, 1, len(avg_df), 1],
        })
        chart1.set_title({'name': 'Average Path Length per Node'})
        chart1.set_x_axis({'name': 'Node ID'})
        chart1.set_y_axis({'name': 'Average Path Length (Jumps)'})
        worksheet.insert_chart('D2', chart1)
        
        # Sheet 2: Percentage Covered
        perc_cov_df = pd.DataFrame(
            [(int(k), v) for k, v in sorted_procentages.items()],
            columns=['Node', 'Average transmission Coverage (%)']
        )
        perc_cov_df.to_excel(writer, sheet_name='Average Transmission Coverage', index=False)
        
        # Add chart for Sheet 2
        worksheet2 = writer.sheets['Average Transmission Coverage']
        chart2 = workbook.add_chart({'type': 'column'})
        chart2.add_series({
            'name': 'Coverage %',
            'categories': ['Average Transmission Coverage', 1, 0, len(perc_cov_df), 0],
            'values': ['Average Transmission Coverage', 1, 1, len(perc_cov_df), 1],
        })
        chart2.set_title({'name': 'Transmission Coverage per Node'})
        chart2.set_x_axis({'name': 'Node ID'})
        chart2.set_y_axis({'name': 'Coverage (%)'})
        worksheet2.insert_chart('D2', chart2)
        
        # Sheet 3: Length Frequencies
        freq_df = pd.DataFrame(
            [[1, 0], [2, 0]] + [(int(k), v) for k, v in sorted_frequencies.items()],
            columns=['Jumps before transmission Death', 'quantity']
        )
        freq_df.to_excel(writer, sheet_name='Point of Transmission Death', index=False)
        
        # Add chart for Sheet 3
        worksheet3 = writer.sheets['Point of Transmission Death']
        chart3 = workbook.add_chart({'type': 'column'})
        chart3.add_series({
            'name': 'Frequency',
            'categories': ['Point of Transmission Death', 1, 0, len(freq_df), 0],
            'values': ['Point of Transmission Death', 1, 1, len(freq_df), 1],
        })
        chart3.set_title({'name': 'Path Length Frequency Distribution'})
        chart3.set_x_axis({'name': 'Path Length (Jumps)'})
        chart3.set_y_axis({'name': 'Frequency'})
        worksheet3.insert_chart('D2', chart3)
        
        # Sheet 4: Probability for Successful Jump
        jump_df = pd.DataFrame(
            [(int(k), v) for k, v in procentage_jumps.items()],
            columns=['Path Length', 'Probability for Successful Jump (%)']
        )
        jump_df.to_excel(writer, sheet_name='Probability for Successful Jump', index=False)
        
        # Add chart for Sheet 4 (line chart)
        worksheet4 = writer.sheets['Probability for Successful Jump']
        chart4 = workbook.add_chart({'type': 'line'})
        chart4.add_series({
            'name': 'Success Probability',
            'categories': ['Probability for Successful Jump', 1, 0, len(jump_df), 0],
            'values': ['Probability for Successful Jump', 1, 1, len(jump_df), 1],
            'marker': {'type': 'circle'},
        })
        chart4.set_title({'name': 'Probability for Successful Jump vs Path Length'})
        chart4.set_x_axis({'name': 'Path Length (Jumps)'})
        chart4.set_y_axis({'name': 'Probability (%)'})
        worksheet4.insert_chart('D2', chart4)
    
    print(f"Excel file with graphs saved!")

#print(sorted_frequencies)
print(f"Procentage for successful jump: {procentage_jumps}")
#print(f"Frequency of each path length: {sorted_frequencies}")
plot_frequency_distribution(procentage_jumps)
save_xlsx_file()
#save_xlsx_file_with_graphs()
#save_csv_files()
