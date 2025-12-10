import pandas as pd
import pickle
from collections import Counter  # Add this import at the top if not already present
import matplotlib.pyplot as plt

nNodes = 50
foldername = "2025-12-10_13-19-38"

nodeData = []

for i in range(nNodes):
    df = pd.read_csv(f"test/{foldername}/gossip_log_node{i}.csv")
    nodeData.append(df)

totalData = pd.read_csv(f"test/{foldername}/final_knowledge.csv")


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
if not os.path.exists(f"test/{foldername}/paths.pkl"):
    paths = get_paths_lengths()
    with open(f"test/{foldername}/paths.pkl", "wb") as f:
        pickle.dump(paths, f)
else:
    with open(f"test/{foldername}/paths.pkl", "rb") as f:
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
    plt.title('Path Length Frequency Distribution')
    plt.show()

lengths_all = get_path_lengths(paths)
frequencies = get_length_frequencies(lengths_all)
sorted_frequencies = dict(sorted(frequencies.items()))
procentage_jumps = get_procentage_for_succesfull_jump(sorted_frequencies)
#print(sorted_frequencies)
print(f"Procentage for successful jump: {procentage_jumps}")
#print(f"Frequency of each path length: {sorted_frequencies}")
plot_frequency_distribution(procentage_jumps)

