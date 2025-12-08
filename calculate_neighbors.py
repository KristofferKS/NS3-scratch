#!/usr/bin/env python3
import json
import math
import sys
from pathlib import Path


def load_placement(placement_file):
    """Load node positions from JSON file."""
    if not Path(placement_file).exists():
        print(f"Error: Placement file '{placement_file}' not found")
        return None
    
    with open(placement_file, 'r') as f:
        placement = json.load(f)
    
    return {int(node_id): pos for node_id, pos in placement.items()}


def calculate_distance(pos1, pos2):
    """Calculate Euclidean distance between two positions."""
    return math.sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)


def calculate_neighbors(placement, comm_range):
    """Calculate neighbor relationships for all nodes."""
    neighbors = {}
    node_ids = sorted(placement.keys())
    
    for node_id in node_ids:
        neighbors[node_id] = []
        node_pos = placement[node_id]
        
        for other_id in node_ids:
            if node_id == other_id:
                continue
            
            other_pos = placement[other_id]
            distance = calculate_distance(node_pos, other_pos)
            
            if distance <= comm_range:
                neighbors[node_id].append(other_id)
    
    return neighbors


def save_neighbors(neighbors, output_file):
    """Save neighbor relationships to JSON file."""
    with open(output_file, 'w') as f:
        json.dump(neighbors, f, indent=2)
    return True


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 calculate_neighbors.py <comm_range> [placement_file] [output_file]")
        return 1
    
    try:
        comm_range = float(sys.argv[1])
    except ValueError:
        print(f"Error: comm_range must be a number")
        return 1
    
    placement_file = sys.argv[2] if len(sys.argv) > 2 else "placement.json"
    output_file = sys.argv[3] if len(sys.argv) > 3 else "neighbors.json"
    
    placement = load_placement(placement_file)
    if placement is None:
        return 1
    
    print(f"Loaded {len(placement)} nodes from {placement_file}")
    
    neighbors = calculate_neighbors(placement, comm_range)
    save_neighbors(neighbors, output_file)
    print(f"Saved neighbor relationships to {output_file}")
    
    total_connections = sum(len(n) for n in neighbors.values())
    avg_neighbors = total_connections / len(neighbors) if neighbors else 0
    
    print(f"Communication range: {comm_range}")
    print(f"Total nodes: {len(neighbors)}")
    print(f"Total connections: {total_connections}")
    print(f"Average neighbors per node: {avg_neighbors:.2f}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
