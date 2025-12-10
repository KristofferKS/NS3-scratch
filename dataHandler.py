import pandas as pd
import pickle

nNodes = 50
foldername = "2025-12-09_15-15-12"

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
        send_to_value = int(path['node_to'].iloc[0])
        if send_to_value not in lengths:
            lengths[send_to_value] = []
        lengths[send_to_value].append(len(path))
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

print(paths[0:10])