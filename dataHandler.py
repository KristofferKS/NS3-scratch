import pandas as pd
import numpy as np

nNodes = 50
foldername = "2025-12-09_11-50-40"

nodeData = []

for i in range(nNodes):
    df = pd.read_csv(f"test/{foldername}/gossip_log_node{i}.csv")
    nodeData.append(df)

totalData = pd.read_csv(f"test/{foldername}/final_knowledge.csv")

total_gen_count = 0
for df in nodeData:
    filtered_df = df[df['event'] == 'gen']
    total_gen_count += len(filtered_df)
#print(f"Total number of 'gen' events across all nodes: {total_gen_count}")

missedEvents = 0

for i, df in enumerate(nodeData):
    gen_events = df[df['event'] == 'gen']
    for index, row in gen_events.iterrows():
        next_index = index + 1
        if next_index < len(df):
            next_event = df.loc[next_index, 'event']
            if next_event != 'send':
                missedEvents += 1

#print(f"Total missed events (gen not followed by send): {missedEvents}")

missedEvents = 0

for i, df in enumerate(nodeData):
    gen_events = df[df['event'] == 'update']
    for index, row in gen_events.iterrows():
        next_index = index + 1
        if next_index < len(df):
            next_event = df.loc[next_index, 'event']
            if next_event != 'send':
                missedEvents += 1

#print(f"Total missed events (update not followed by send): {missedEvents}")

#print(totalData.iloc[:, :11].to_string(na_rep=''))

matrix = pd.DataFrame(index=range(nNodes), columns=range(nNodes), dtype=float)
for i in range(nNodes):
    if i in totalData.index:
        row = totalData.loc[i]
        for j in range(nNodes):
            times = []
            for k in range(5):  # Assuming max 5 entries per node
                col_name = f'node{j}time{k}'
                if col_name in row.index and not pd.isna(row[col_name]):
                    times.append(row[col_name])
            if times:
                matrix.loc[i, j] = max(times)
            else:
                matrix.loc[i, j] = np.nan

#print("Knowledge matrix (newest timestamp each node has on others):")
#print(matrix.to_string(na_rep=''))
matrix.to_csv(f"test/{foldername}/knowledge_matrix.csv")
#print(f"Matrix saved to test/{foldername}/knowledge_matrix.csv")

def get_n_event(event_name: str, node_table: pd.DataFrame, event_number:int=0, match_value=None):
    events = node_table[node_table['event'] == event_name]
    if not events.empty and event_number < len(events):
        event = events.iloc[event_number]
    else:
        return None
    if match_value is None:
        return event
    else:
        if event['timestamp_ms'] == match_value:
            return event
        else:
            return None

def get_entries_by_timestamp(node_table: pd.DataFrame, timestamp_value):
    return node_table[node_table['timestamp_ms'] == timestamp_value]

def convert_node_to(node_to_str):
    if isinstance(node_to_str, str) and ':' in node_to_str:
        parts = node_to_str.split(':')
        last = parts[-1]
        return int(last, 16) - 1
    else:
        return int(node_to_str) - 1
class path:
    def __init__(self, genNumber, nodeIndex):
        self.genNumber = genNumber
        self.nodeIndex = nodeIndex
        self.entries = []  # List to collect entries
        self.get_nth_gen_from_node(genNumber, nodeIndex)
    
    def print_entries(self):
        if self.entries:
            df = pd.DataFrame(self.entries)
            print(f"Pretty-printed entries for gen {self.genNumber}, node {self.nodeIndex}:")
            print(df.to_string(index=False))
        else:
            print("No entries to print")

    def get_nth_gen_from_node(self, n, node_index):
        if node_index < len(nodeData):
            event = get_n_event('gen', nodeData[node_index], event_number=n)
            if event is not None:
                self.entries.append(event.to_dict())  # Collect the gen event
                # Get all entries with the same timestamp_ms from ALL nodeData
                timestamp = event['timestamp_ms']
                for i in range(nNodes):
                    all_entries_with_timestamp = get_entries_by_timestamp(nodeData[i], timestamp)
                    if not all_entries_with_timestamp.empty:
                        # Collect the entries
                        self.entries.extend(all_entries_with_timestamp.to_dict('records'))

test = path(0, 0)
test.print_entries()  # Pretty-print the entries