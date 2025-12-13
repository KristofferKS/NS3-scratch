import pandas as pd
import pickle
from collections import Counter
import matplotlib.pyplot as plt
import os
import re
from multiprocessing import Pool

nNodes = 50
folder_name = 'test/Temp_test'

# Global variables for multiprocessing
_gen_list = None
_updates_grouped = None

def init_worker(gen_list, updates_grouped):
    """Initialize worker process with shared data"""
    global _gen_list, _updates_grouped
    _gen_list = gen_list
    _updates_grouped = updates_grouped

def get_path_fast(gen_number):
    """Worker function for multiprocessing"""
    gen_info = _gen_list[gen_number]
    message_id = gen_info['node_from']
    timestamp_ms = gen_info['timestamp_ms']
    
    try:
        stepByStep_df = _updates_grouped.get_group((message_id, timestamp_ms))
        stepByStep_df = stepByStep_df.sort_values(by='time', ignore_index=True)
        return stepByStep_df
    except KeyError:
        return pd.DataFrame()

def get_connections_from_folder(folder_name):
    """Extract connection number from folder name like '50-sim3' or '258.5'"""
    match = re.search(r'(\d+(?:\.\d+)?)', folder_name)
    if match:
        return float(match.group(1))
    return None

def process_folder(foldername, total, current):
    """Process a single simulation folder and return the probability jumps data, coverage, and drop rate"""
    print(f"\nProcessing folder: {foldername} | {current}/{total}")
    
    # Extract connection number from folder name
    connections = get_connections_from_folder(foldername)
    if connections is None:
        print(f"Warning: Could not extract connection number from {foldername}")
        connections = 0
    
    # Check if cached metrics exist
    metrics_cache_path = f"{folder_name}/{foldername}/metrics_cache.pkl"
    if os.path.exists(metrics_cache_path):
        print(f"Loading cached metrics from {metrics_cache_path}")
        with open(metrics_cache_path, "rb") as f:
            cached_metrics = pickle.load(f)
        return (cached_metrics['procentage_jumps'], 
                cached_metrics['percentage_covered'], 
                cached_metrics['drop_rate'], 
                cached_metrics['avg_coverage'], 
                connections)
    
    nodeData = []
    
    # Load all node data
    for i in range(nNodes):
        csv_path = f"{folder_name}/{foldername}/gossip_log_node{i}.csv"
        if not os.path.exists(csv_path):
            print(f"Warning: {csv_path} not found, skipping folder")
            return None, None, None, None, connections
        df = pd.read_csv(csv_path)
        nodeData.append(df)
    
    combined_df = pd.concat(nodeData, ignore_index=True)
    
    # Calculate drop rate from combined data
    total_send_events = len(combined_df[combined_df['event'] == 'send'])
    total_drop_events = len(combined_df[combined_df['event'] == 'drop'])
    
    if total_send_events > 0:
        drop_rate = (total_drop_events / total_send_events) * 100
    else:
        drop_rate = 0
    
    print(f"  Total send events: {total_send_events}")
    print(f"  Total drop events: {total_drop_events}")
    print(f"  Drop rate: {drop_rate:.2f}%")
    
    # Precompute gen info and group updates
    gen_df = combined_df[combined_df['event'] == 'gen']
    gen_list = gen_df[['timestamp_ms', 'node_from']].to_dict('records')
    updates_df = combined_df[combined_df['event'] == 'update']
    updates_grouped = updates_df.groupby(['node_from', 'timestamp_ms'])
    
    # Check if paths.pkl exists
    pkl_path = f"{folder_name}/{foldername}/paths.pkl"
    if os.path.exists(pkl_path):
        print(f"Loading cached paths from {pkl_path}")
        with open(pkl_path, "rb") as f:
            paths = pickle.load(f)
    else:
        print(f"Computing paths for {len(gen_list)} gen events...")
        
        # Use multiprocessing
        num_cores = os.cpu_count()
        print(f"Using {num_cores} cores for parallel processing")
        
        with Pool(num_cores, initializer=init_worker, initargs=(gen_list, updates_grouped)) as pool:
            paths = []
            total = len(gen_list)
            for i, path in enumerate(pool.imap(get_path_fast, range(total)), 1):
                if i % 1000 == 0:
                    print(f"Processed {i}/{total} ({i/total*100:.1f}%)")
                paths.append(path)
        
        print(f"Finished processing {total} gen events")
        
        # Save for future use
        with open(pkl_path, "wb") as f:
            pickle.dump(paths, f)
    
    # Calculate path lengths - OPTIMIZED VERSION
    def get_path_lengths(paths):
        """Calculate how many nodes each message reached - OPTIMIZED"""
        lengths = {}
        for path in paths:
            if not path.empty:
                # Get the originating node from the first row
                originating_node = int(path['node_from'].iloc[0])
                
                # Count all unique nodes that received the message (FAST - vectorized)
                unique_nodes_reached = set(path['node_to'].unique())
                unique_nodes_reached.add(originating_node)  # Include sender
                
                if originating_node not in lengths:
                    lengths[originating_node] = []
                
                # Store how many nodes this message reached
                lengths[originating_node].append(len(unique_nodes_reached))
        
        return lengths
    
    def get_average_path_length(lengths_dict):
        """Average number of nodes reached per originating node"""
        averages = {}
        # Ensure all nodes are represented
        for node_id in range(nNodes):
            if node_id in lengths_dict and lengths_dict[node_id]:
                averages[node_id] = sum(lengths_dict[node_id]) / len(lengths_dict[node_id])
            else:
                averages[node_id] = 0
        return averages
    
    def get_percentage_covered(average_lengths, total_nodes=nNodes):
        """Percentage of network reached on average"""
        percentages = {}
        for node_id in range(total_nodes):
            avg = average_lengths.get(node_id, 0)
            if avg > 0:
                percentages[node_id] = (avg / total_nodes) * 100
            else:
                percentages[node_id] = 0
        return percentages
    
    def get_length_frequencies(lengths_dict):
        all_lengths = []
        for lengths_list in lengths_dict.values():
            all_lengths.extend(lengths_list)
        return dict(Counter(all_lengths))
    
    def get_procentage_for_succesfull_jump(length_frequencies):
        total_paths = sum(length_frequencies.values())
        cumulative = 0
        succesfull_jumps = {'1': 100, '2': 100}
        for length, freq in length_frequencies.items():
            succesfull_jumps[str(length - 1)] = ((total_paths - cumulative) / total_paths) * 100
            cumulative += freq
        return succesfull_jumps
    
    print("Calculating path lengths...")
    lengths_all = get_path_lengths(paths)
    print("Calculating averages...")
    average_lengths = get_average_path_length(lengths_all)
    print("Calculating coverage...")
    percentage_covered = get_percentage_covered(average_lengths)
    
    # Calculate average coverage across all nodes
    avg_coverage = sum(percentage_covered.values()) / len(percentage_covered)
    print(f"  Average coverage: {avg_coverage:.2f}%")
    
    print("Calculating frequencies...")
    frequencies = get_length_frequencies(lengths_all)
    sorted_frequencies = dict(sorted(frequencies.items()))
    print("Calculating jump probabilities...")
    procentage_jumps = get_procentage_for_succesfull_jump(sorted_frequencies)
    
    # Cache all computed metrics
    metrics = {
        'procentage_jumps': procentage_jumps,
        'percentage_covered': percentage_covered,
        'drop_rate': drop_rate,
        'avg_coverage': avg_coverage
    }
    
    print(f"Saving metrics cache to {metrics_cache_path}")
    with open(metrics_cache_path, "wb") as f:
        pickle.dump(metrics, f)
    
    return procentage_jumps, percentage_covered, drop_rate, avg_coverage, connections

def create_combined_excel():
    """Create Excel file with combined data from all simulation folders"""
    
    # Find all simulation folders
    sim_folders = []
    base_path = folder_name
    
    if not os.path.exists(base_path):
        print(f"Error: {base_path} directory not found")
        return
    
    for folder in os.listdir(base_path):
        folder_path = os.path.join(base_path, folder)
        if os.path.isdir(folder_path):
            sim_folders.append(folder)
    
    if not sim_folders:
        print("No simulation folders found")
        return
    
    print(f"Found {len(sim_folders)} simulation folders: {sim_folders}")
    
    # Process all folders
    all_jump_data = {}
    all_coverage_data = {}
    all_drop_rates = {}
    all_avg_coverage = {}
    folder_connections = {}
    
    folder_amount = len(sim_folders)
    folders_started = 1
    for folder in sim_folders:
        procentage_jumps, percentage_covered, drop_rate, avg_coverage, connections = process_folder(folder,folder_amount,folders_started)
        folders_started += 1
        if procentage_jumps is not None:
            all_jump_data[folder] = procentage_jumps
            all_coverage_data[folder] = percentage_covered
            all_drop_rates[folder] = drop_rate
            all_avg_coverage[folder] = avg_coverage
            folder_connections[folder] = connections
    
    if not all_jump_data:
        print("No data processed successfully")
        return
    
    # Sort folders by connection number
    sorted_folders = sorted(all_jump_data.keys(), key=lambda x: folder_connections[x])
    
    # Create combined DataFrame for Probability jumps (Sheet 1)
    path_lengths = list(range(1, nNodes))  # 1 to 49
    jump_data = {'Path Length': path_lengths}
    
    for folder in sorted_folders:
        procentage_jumps = all_jump_data[folder]
        connections = folder_connections[folder]
        col_name = f'prob {connections}'
        
        probabilities = []
        for length in path_lengths:
            prob = procentage_jumps.get(str(length), 0)
            probabilities.append(prob)
        
        jump_data[col_name] = probabilities
    
    jump_df = pd.DataFrame(jump_data)
    
    # Create combined DataFrame for Coverage (Sheet 2)
    nodes = list(range(nNodes))  # 0 to 49
    coverage_data = {'Node': nodes}
    
    for folder in sorted_folders:
        percentage_covered = all_coverage_data[folder]
        connections = folder_connections[folder]
        col_name = f'prob {connections}'
        
        coverages = []
        for node in nodes:
            cov = percentage_covered.get(node, 0)
            coverages.append(cov)
        
        coverage_data[col_name] = coverages
    
    coverage_df = pd.DataFrame(coverage_data)
    
    # Create DataFrame for Coverage/Drop Rate Ratio (Sheet 3)
    ratio_data = {
        'Gossip Probability (%)': [folder_connections[folder] * 100 for folder in sorted_folders],
        'Average Coverage (%)': [all_avg_coverage[folder] for folder in sorted_folders],
        'Drop Rate (%)': [all_drop_rates[folder] for folder in sorted_folders],
    }
    
    # Calculate ratio (handle division by zero)
    ratios = []
    for folder in sorted_folders:
        coverage = all_avg_coverage[folder]
        drop_rate = all_drop_rates[folder]
        if drop_rate > 0:
            ratio = coverage / drop_rate
        else:
            ratio = coverage  # If no drops, ratio is just coverage
        ratios.append(ratio)
    
    ratio_data['Coverage/Drop Rate Ratio'] = ratios
    ratio_df = pd.DataFrame(ratio_data)
    
    # Save to Excel with graphs
    output_file = f'{folder_name}/combined_analysis4.xlsx'
    
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        # Sheet 1: Probability Comparison
        jump_df.to_excel(writer, sheet_name='Probability Comparison', index=False)
        
        workbook = writer.book
        worksheet1 = writer.sheets['Probability Comparison']
        
        # Create line chart for probabilities
        chart1 = workbook.add_chart({'type': 'line'})
        
        for col_idx, folder in enumerate(sorted_folders, start=1):
            connections = folder_connections[folder]
            chart1.add_series({
                'name': f'prob {connections}',
                'categories': ['Probability Comparison', 1, 0, len(path_lengths), 0],
                'values': ['Probability Comparison', 1, col_idx, len(path_lengths), col_idx],
                'marker': {'type': 'circle', 'size': 5},
            })
        
        chart1.set_title({'name': 'Probability for Successful Jump Comparison',
                         'name_font': {'size': 20, 'bold': False}})
        chart1.set_x_axis({'name': 'Path Length (Jumps)',
                          'name_font': {'size': 14, 'bold': False}})
        chart1.set_y_axis({'name': 'Probability (%)',
                          'name_font': {'size': 14, 'bold': False},
                          'major_unit': 10,
                          'min': 0,
                          'max': 100})
        chart1.set_size({'width': 720, 'height': 480})
        chart1.set_legend({'position': 'right'})
        
        worksheet1.insert_chart('E2', chart1)
        
        # Sheet 2: Coverage Comparison
        coverage_df.to_excel(writer, sheet_name='Coverage Comparison', index=False)
        
        worksheet2 = writer.sheets['Coverage Comparison']
        
        # Create line chart for coverage
        chart2 = workbook.add_chart({'type': 'line'})
        
        for col_idx, folder in enumerate(sorted_folders, start=1):
            connections = folder_connections[folder]
            chart2.add_series({
                'name': f'prob {connections}',
                'categories': ['Coverage Comparison', 1, 0, len(nodes), 0],
                'values': ['Coverage Comparison', 1, col_idx, len(nodes), col_idx],
                'marker': {'type': 'circle', 'size': 5},
            })
        
        chart2.set_title({'name': 'Average Transmission Coverage Comparison',
                         'name_font': {'size': 20, 'bold': False}})
        chart2.set_x_axis({'name': 'Node ID',
                          'name_font': {'size': 14, 'bold': False}})
        chart2.set_y_axis({'name': 'Coverage (%)',
                          'name_font': {'size': 14, 'bold': False},
                          'major_unit': 10,
                          'min': 0,
                          'max': 100})
        chart2.set_size({'width': 720, 'height': 480})
        chart2.set_legend({'position': 'right'})
        
        worksheet2.insert_chart('E2', chart2)
        
        # Sheet 3: Coverage/Drop Rate Ratio
        ratio_df.to_excel(writer, sheet_name='Efficiency Analysis', index=False)
        
        worksheet3 = writer.sheets['Efficiency Analysis']
        
        # Create line chart for efficiency ratio
        chart3 = workbook.add_chart({'type': 'line'})
        
        chart3.add_series({
            'name': 'Coverage/Drop Rate Ratio',
            'categories': ['Efficiency Analysis', 1, 0, len(sorted_folders), 0],
            'values': ['Efficiency Analysis', 1, 3, len(sorted_folders), 3],
            'marker': {'type': 'circle', 'size': 7},
            'line': {'width': 2.5},
        })
        
        chart3.set_title({'name': 'Gossip Efficiency: Coverage/Drop Rate Ratio',
                         'name_font': {'size': 20, 'bold': False}})
        chart3.set_x_axis({'name': 'Gossip Probability (%)',
                          'name_font': {'size': 14, 'bold': False}})
        chart3.set_y_axis({'name': 'Coverage / Drop Rate Ratio',
                          'name_font': {'size': 14, 'bold': False}})
        chart3.set_size({'width': 720, 'height': 480})
        chart3.set_legend({'position': 'none'})
        
        worksheet3.insert_chart('F2', chart3)
    
    print(f"\nCombined analysis saved to {output_file}")
    print(f"Processed {len(all_jump_data)} folders: {sorted_folders}")
    
    # Print summary
    print("\n=== Efficiency Summary ===")
    for folder in sorted_folders:
        prob = folder_connections[folder] * 100
        coverage = all_avg_coverage[folder]
        drop_rate = all_drop_rates[folder]
        if drop_rate > 0:
            ratio = coverage / drop_rate
        else:
            ratio = coverage
        print(f"Gossip {prob}%: Coverage={coverage:.2f}%, Drop Rate={drop_rate:.2f}%, Ratio={ratio:.2f}")

# Run the combined analysis
if __name__ == "__main__":
    create_combined_excel()