""" Create a gui where nodes are able to be placed, 
the bounds of the space should be written from the c++ 
file that calls this python program as well as the number 
of nodes that should be place, and the range of the nodes. 
The nodes should be diplayed as red dots, and there should 
be a circle illustrating range. The program should be able 
to detect if all nodes have atleast two neighbour nodes in 
range, and if there is a redundant node. Then it need to output 
the coordinates, so the c++ program can use them """

import sys
import tkinter as tk
from tkinter import messagebox
import json
import math

num_nodes = int(sys.argv[1])
wifi_range = int(float(sys.argv[2]))
centerX = float(sys.argv[3])
centerY = float(sys.argv[4])
grid_range = int(float(sys.argv[5]))
nodes = {}
range_circles = {}  # Store range circle IDs for each node

# Calculate scale factor: meters to pixels
# If grid_range is 200m and we want it to fit in a reasonable window
SCALE = 2  # 2 pixels per meter (adjust as needed for your screen)
CANVAS_SIZE = grid_range * 2 * SCALE
range_pixels = wifi_range * SCALE  # Convert range to pixels

def add_node(event):
    x, y = event.x, event.y
    radius = 5
    
    # Check if node is being placed on top of another node
    overlapping = m.find_overlapping(x-radius, y-radius, x+radius, y+radius)
    # Filter to only check ovals (nodes)
    node_overlap = [item for item in overlapping if "oval" in m.gettags(item)]
    
    if node_overlap:
        remove_node(event, node_overlap)
    else:
        if len(nodes) < num_nodes:
            # Draw the range circle first (so it's behind the node)
            range_circle = m.create_oval(
                x - range_pixels, y - range_pixels,
                x + range_pixels, y + range_pixels,
                outline="lightblue", width=2, dash=(4, 4)
            )
            m.addtag_withtag("range_circle", range_circle)
            
            # Draw a small circle to visualize the node
            node_ID = m.create_oval(x-radius, y-radius, x+radius, y+radius, fill="red")
            
            # Save the coordinates and tags
            node_tag = f"{len(nodes)}"  # Start from 0
            m.addtag_withtag(node_tag, node_ID)
            m.addtag_withtag("oval", node_ID)
            
            # Add label
            oval_label = m.create_text((x, y+12), text=node_tag, font=("Arial", 10, "bold"))
            m.addtag_withtag("label", oval_label)
            m.addtag_withtag(f"label_{node_tag}", oval_label)
            
            # Store node position and range circle ID
            nodes[node_tag] = (x, y)
            range_circles[node_tag] = range_circle
            
            # Update neighbor count display
            update_neighbor_info()
        else:
            tk.messagebox.showinfo("Pop-Up", "No more nodes to place")

def remove_node(event, node_IDs=None):
    if node_IDs is None:
        return
    
    node_ID = node_IDs[0]
    tags = m.gettags(node_ID)
    
    # Find the node tag (it's a number)
    node_tag = None
    for tag in tags:
        if tag.isdigit():
            node_tag = tag
            break
    
    if node_tag:
        # Remove from nodes dictionary
        nodes.pop(node_tag, None)
        
        # Delete range circle
        if node_tag in range_circles:
            m.delete(range_circles[node_tag])
            del range_circles[node_tag]
        
        # Delete label
        label_items = m.find_withtag(f"label_{node_tag}")
        for item in label_items:
            m.delete(item)
        
        # Delete the node itself
        m.delete(node_ID)
        
        # Update neighbor count display
        update_neighbor_info()

def calculate_distance(pos1, pos2):
    """Calculate Euclidean distance between two positions"""
    return math.sqrt((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)

def update_neighbor_info():
    """Update the info display showing neighbor counts"""
    # Count neighbors for each node
    neighbor_counts = {}
    for node_tag, pos in nodes.items():
        count = 0
        for other_tag, other_pos in nodes.items():
            if node_tag != other_tag:
                dist = calculate_distance(pos, other_pos)
                if dist <= range_pixels:
                    count += 1
        neighbor_counts[node_tag] = count
    
    # Update info text
    info_text = f"Nodes placed: {len(nodes)}/{num_nodes}\n"
    info_text += f"WiFi Range: {wifi_range}m\n\n"
    
    if nodes:
        info_text += "Neighbor counts:\n"
        for node_tag in sorted(nodes.keys(), key=int):
            count = neighbor_counts[node_tag]
            status = "✓" if count >= 2 else "✗"
            info_text += f"Node {node_tag}: {count} neighbors {status}\n"
    
    info_label.config(text=info_text)
    
    # Update the scroll region
    info_canvas.configure(scrollregion=info_canvas.bbox("all"))

def finish_placement():
    """Finish placement and output coordinates"""
    if len(nodes) < num_nodes:
        result = tk.messagebox.askyesno(
            "Incomplete Placement",
            f"Only {len(nodes)}/{num_nodes} nodes placed. Continue anyway?"
        )
        if not result:
            return
    
    # Convert pixel coordinates back to meters (relative to center)
    output_nodes = {}
    for node_tag, (px, py) in nodes.items():
        # Convert from canvas coordinates to meters
        x_meters = (px / SCALE) - grid_range
        y_meters = (py / SCALE) - grid_range
        output_nodes[node_tag] = [x_meters, y_meters]
    
    # Sort by node ID (numerically) before saving
    sorted_output = {k: output_nodes[k] for k in sorted(output_nodes.keys(), key=int)}

    
    # Save to file
    with open("placement.json", "w") as f:
        f.write(json.dumps(sorted_output))
    
    # Print to stdout for C++ to read
    print(json.dumps(sorted_output))
    
    m.quit()

def _on_mousewheel(event):
    """Handle mousewheel scrolling"""
    info_canvas.yview_scroll(int(-1*(event.delta/120)), "units")

# Create main window
root = tk.Tk()
root.title("Node Placement - Click to place nodes")

# Create frame for canvas and info
main_frame = tk.Frame(root)
main_frame.pack()

# Create canvas
m = tk.Canvas(main_frame, width=CANVAS_SIZE, height=CANVAS_SIZE, bg="white")
m.grid(row=0, column=0)

# Draw grid border
m.create_rectangle(2, 2, CANVAS_SIZE-2, CANVAS_SIZE-2, outline="gray", width=2)

# Add center point reference
center_px = CANVAS_SIZE / 2
m.create_line(center_px - 10, center_px, center_px + 10, center_px, fill="gray", width=2)
m.create_line(center_px, center_px - 10, center_px, center_px + 10, fill="gray", width=2)
m.create_text(center_px + 20, center_px - 20, text="Center (0,0)", fill="gray")

# Create scrollable info panel
info_outer_frame = tk.Frame(main_frame, bg="lightgray")
info_outer_frame.grid(row=0, column=1, sticky="ns")

# Add canvas for scrolling
info_canvas = tk.Canvas(info_outer_frame, bg="lightgray", width=250, height=CANVAS_SIZE, highlightthickness=0)
info_canvas.pack(side="left", fill="both", expand=True)

# Add scrollbar
scrollbar = tk.Scrollbar(info_outer_frame, orient="vertical", command=info_canvas.yview)
scrollbar.pack(side="right", fill="y")
info_canvas.configure(yscrollcommand=scrollbar.set)

# Create frame inside canvas
info_frame = tk.Frame(info_canvas, bg="lightgray", padx=10, pady=10)
info_canvas_frame = info_canvas.create_window((0, 0), window=info_frame, anchor="nw")

info_label = tk.Label(
    info_frame,
    text=f"Nodes placed: 0/{num_nodes}\nWiFi Range: {wifi_range}m",
    font=("Arial", 10),
    bg="lightgray",
    justify="left",
    anchor="w"
)
info_label.pack()

# Add instructions
instructions = tk.Label(
    info_frame,
    text="\nInstructions:\n"
         "• Left-click to place node\n"
         "• Right-click to remove node\n"
         "• Blue circles show range\n"
         "• Need ≥2 neighbors per node\n",
    font=("Arial", 9),
    bg="lightgray",
    justify="left",
    anchor="w"
)
instructions.pack(pady=10)

# Add finish button
finish_btn = tk.Button(
    info_frame,
    text="Finish Placement",
    command=finish_placement,
    bg="green",
    fg="white",
    font=("Arial", 10, "bold"),
    padx=20,
    pady=10
)
finish_btn.pack(pady=10)

# Bind mousewheel for scrolling
info_canvas.bind_all("<MouseWheel>", _on_mousewheel)

# Update scroll region after frame is created
info_frame.update_idletasks()
info_canvas.configure(scrollregion=info_canvas.bbox("all"))

# Bind mouse events
m.bind("<Button-1>", add_node)
m.bind("<Button-3>", lambda event: remove_node(event, m.find_overlapping(event.x, event.y, event.x, event.y)))

root.mainloop()