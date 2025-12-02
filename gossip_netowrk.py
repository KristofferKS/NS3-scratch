#!/usr/bin/env python3
"""
NS-3 3.46 Peer-to-Peer Gossip Protocol Network Simulation
Place in: ns-3.46/scratch/gossip_network.py
Run with: ./ns3 run scratch/gossip_network.py
"""

import sys
import json
import random

# NS-3 imports
try:
    from ns import ns
except ImportError:
    print("Error: NS-3 Python bindings not found!")
    print("Make sure you've built NS-3 with Python bindings enabled.")
    sys.exit(1)


class DensityData:
    """Represents a density data vector"""
    def __init__(self, timestamp, node_id, density):
        self.timestamp = timestamp
        self.node_id = node_id
        self.density = density
    
    def to_dict(self):
        return {
            'timestamp': self.timestamp,
            'node_id': self.node_id,
            'density': self.density
        }
    
    @staticmethod
    def from_dict(data):
        return DensityData(data['timestamp'], data['node_id'], data['density'])
    
    def __repr__(self):
        return f"DensityData(t={self.timestamp:.2f}, node={self.node_id}, d={self.density})"


class GossipApplication(ns.applications.Application):
    """NS-3 Application that implements gossip protocol for density data sharing"""
    
    # Class variables to track all applications
    all_apps = []
    node_data_stores = {}
    
    def __init__(self):
        super().__init__()
        self.node_id = -1
        self.socket = None
        self.neighbors = []  # List of IPv4 addresses
        self.data_store = {}  # Maps node_id -> DensityData
        self.gossip_fanout = 3
        self.gossip_interval = 1.0
        self.density_update_interval = 3.0
        self.port = 9999
        
        self.gossip_event = None
        self.density_event = None
        
    def Setup(self, node_id, neighbors):
        """Configure the application"""
        self.node_id = node_id
        self.neighbors = neighbors
        GossipApplication.all_apps.append(self)
        GossipApplication.node_data_stores[node_id] = self.data_store
        
    def StartApplication(self):
        """Called when application starts"""
        # Create UDP socket
        tid = ns.core.TypeId.LookupByName("ns3::UdpSocketFactory")
        self.socket = ns.network.Socket.CreateSocket(self.GetNode(), tid)
        
        # Bind to any address on the specified port
        local_addr = ns.network.InetSocketAddress(ns.network.Ipv4Address.GetAny(), self.port)
        self.socket.Bind(local_addr)
        self.socket.SetRecvCallback(self.HandleRead)
        
        # Enable broadcast
        self.socket.SetAllowBroadcast(True)
        
        print(f"[{ns.core.Simulator.Now().GetSeconds():.2f}s] Node {self.node_id} started with {len(self.neighbors)} neighbors")
        
        # Schedule first operations
        self.ScheduleDensityUpdate()
        self.ScheduleGossip()
        
    def StopApplication(self):
        """Called when application stops"""
        if self.socket:
            self.socket.Close()
            self.socket = None
        
        # Cancel pending events
        if self.gossip_event and not self.gossip_event.IsExpired():
            ns.core.Simulator.Cancel(self.gossip_event)
        if self.density_event and not self.density_event.IsExpired():
            ns.core.Simulator.Cancel(self.density_event)
    
    def ScheduleDensityUpdate(self):
        """Schedule periodic density measurement updates"""
        self.UpdateOwnDensity()
        
        # Schedule next update
        next_time = self.density_update_interval + random.uniform(-0.5, 0.5)
        self.density_event = ns.core.Simulator.Schedule(
            ns.core.Seconds(next_time),
            self.ScheduleDensityUpdate
        )
    
    def UpdateOwnDensity(self):
        """Generate new density measurement"""
        density = random.randint(10, 100)
        timestamp = ns.core.Simulator.Now().GetSeconds()
        
        new_data = DensityData(timestamp, self.node_id, density)
        self.UpdateData(new_data)
        
    def ScheduleGossip(self):
        """Schedule periodic gossip to neighbors"""
        self.GossipToNeighbors()
        
        # Schedule next gossip
        self.gossip_event = ns.core.Simulator.Schedule(
            ns.core.Seconds(self.gossip_interval),
            self.ScheduleGossip
        )
    
    def UpdateData(self, new_data):
        """Update data store with new data if it's newer"""
        node_id = new_data.node_id
        
        if node_id not in self.data_store or \
           self.data_store[node_id].timestamp < new_data.timestamp:
            self.data_store[node_id] = new_data
            print(f"[{ns.core.Simulator.Now().GetSeconds():.2f}s] Node {self.node_id} updated: {new_data}")
            return True
        return False
    
    def GossipToNeighbors(self):
        """Send data to randomly selected neighbors"""
        if not self.neighbors or not self.socket:
            return
        
        # Select random subset of neighbors
        num_to_select = min(self.gossip_fanout, len(self.neighbors))
        selected_neighbors = random.sample(self.neighbors, num_to_select)
        
        # Prepare message
        message = {
            'type': 'gossip',
            'sender': self.node_id,
            'data': [data.to_dict() for data in self.data_store.values()]
        }
        
        message_str = json.dumps(message)
        message_bytes = message_str.encode('utf-8')
        
        # Send to selected neighbors
        for neighbor_addr in selected_neighbors:
            packet = ns.network.Packet(len(message_bytes))
            packet.AddHeader(ns.network.SeqTsHeader())
            
            # Copy data into packet
            buffer_data = ns.network.Buffer(len(message_bytes))
            buffer_data.AddAtStart(len(message_bytes))
            iterator = buffer_data.Begin()
            for byte in message_bytes:
                iterator.WriteU8(byte)
                iterator.Next()
            
            packet = ns.network.Packet(buffer_data)
            
            # Send packet
            remote_addr = ns.network.InetSocketAddress(neighbor_addr, self.port)
            self.socket.SendTo(packet, 0, remote_addr)
        
        if selected_neighbors:
            neighbor_ids = [self.GetNodeIdFromAddress(addr) for addr in selected_neighbors]
            print(f"[{ns.core.Simulator.Now().GetSeconds():.2f}s] Node {self.node_id} gossiping to nodes: {neighbor_ids}")
    
    def HandleRead(self, socket):
        """Handle received packets"""
        address = ns.network.Address()
        packet = socket.RecvFrom(address)
        
        if not packet:
            return
        
        # Extract data from packet
        size = packet.GetSize()
        buffer_data = ns.network.Buffer(size)
        packet.CopyData(buffer_data, size)
        
        # Read bytes from buffer
        message_bytes = bytearray()
        iterator = buffer_data.Begin()
        for i in range(size):
            message_bytes.append(iterator.ReadU8())
            if not iterator.IsEnd():
                iterator.Next()
        
        try:
            message_str = message_bytes.decode('utf-8')
            message = json.loads(message_str)
            
            if message['type'] == 'gossip':
                sender_id = message['sender']
                print(f"[{ns.core.Simulator.Now().GetSeconds():.2f}s] Node {self.node_id} received gossip from Node {sender_id}")
                
                # Process received data
                updated = False
                for data_dict in message['data']:
                    new_data = DensityData.from_dict(data_dict)
                    if self.UpdateData(new_data):
                        updated = True
                
                # Propagate new data immediately (epidemic spreading)
                if updated:
                    delay = random.uniform(0.1, 0.3)
                    ns.core.Simulator.Schedule(
                        ns.core.Seconds(delay),
                        self.GossipToNeighbors
                    )
        
        except Exception as e:
            print(f"Error processing packet: {e}")
    
    def GetNodeIdFromAddress(self, ipv4_addr):
        """Helper to get node ID from IP address (assumes 10.1.1.x mapping)"""
        # This is a simplified mapping - adjust based on your addressing scheme
        addr_str = str(ipv4_addr)
        parts = addr_str.split('.')
        if len(parts) == 4:
            return int(parts[3]) - 1  # Assuming 10.1.1.1 is node 0
        return -1


def CalculateNeighbors(nodes, interfaces, transmission_range):
    """Calculate which nodes are neighbors based on distance"""
    num_nodes = nodes.GetN()
    neighbor_lists = [[] for _ in range(num_nodes)]
    
    for i in range(num_nodes):
        node_i = nodes.Get(i)
        mob_i = node_i.GetObject(ns.mobility.MobilityModel.GetTypeId())
        pos_i = mob_i.GetPosition()
        
        for j in range(num_nodes):
            if i == j:
                continue
            
            node_j = nodes.Get(j)
            mob_j = node_j.GetObject(ns.mobility.MobilityModel.GetTypeId())
            pos_j = mob_j.GetPosition()
            
            # Calculate distance
            dx = pos_j.x - pos_i.x
            dy = pos_j.y - pos_i.y
            distance = (dx*dx + dy*dy) ** 0.5
            
            if distance <= transmission_range:
                # Get IP address of neighbor
                ip_addr = interfaces.GetAddress(j)
                neighbor_lists[i].append(ip_addr)
    
    return neighbor_lists


def PrintFinalStatistics():
    """Print final convergence statistics"""
    print("\n" + "="*80)
    print("FINAL STATISTICS - Data Convergence Analysis")
    print("="*80)
    
    num_nodes = len(GossipApplication.all_apps)
    
    # Collect all unique data entries
    all_timestamps = set()
    for node_id, data_store in GossipApplication.node_data_stores.items():
        for data in data_store.values():
            all_timestamps.add((data.node_id, data.timestamp))
    
    total_unique_entries = len(all_timestamps)
    
    print(f"\nTotal nodes: {num_nodes}")
    print(f"Total unique data entries generated: {total_unique_entries}")
    
    # Check convergence for each node
    print(f"\n{'Node ID':<10} {'Entries':<10} {'Completeness':<15} {'Status'}")
    print("-" * 60)
    
    all_converged = True
    for node_id in sorted(GossipApplication.node_data_stores.keys()):
        data_store = GossipApplication.node_data_stores[node_id]
        num_entries = len(data_store)
        completeness = (num_entries / total_unique_entries * 100) if total_unique_entries > 0 else 0
        status = "✓ CONVERGED" if num_entries == total_unique_entries else "✗ INCOMPLETE"
        
        if num_entries != total_unique_entries:
            all_converged = False
        
        print(f"{node_id:<10} {num_entries:<10} {completeness:>6.1f}%{'':<7} {status}")
    
    print("\n" + "="*80)
    if all_converged:
        print("SUCCESS: All nodes have converged to identical global state!")
    else:
        print("PARTIAL: Some nodes have not fully converged yet.")
    print("="*80 + "\n")


def main():
    """Main simulation function"""
    
    # Simulation parameters
    num_nodes = 10
    simulation_time = 30.0
    transmission_range = 150.0
    
    print("="*80)
    print("NS-3 Gossip Protocol P2P Network Simulation")
    print("="*80)
    print(f"Number of nodes: {num_nodes}")
    print(f"Simulation time: {simulation_time}s")
    print(f"Transmission range: {transmission_range}m")
    print("="*80 + "\n")
    
    # Enable logging (optional)
    # ns.core.LogComponentEnable("UdpSocketImpl", ns.core.LOG_LEVEL_INFO)
    
    # Create nodes
    nodes = ns.network.NodeContainer()
    nodes.Create(num_nodes)
    
    # Set up WiFi Ad-hoc network
    wifi = ns.wifi.WifiHelper()
    wifi.SetStandard(ns.wifi.WIFI_STANDARD_80211b)
    
    # PHY layer
    phy = ns.wifi.YansWifiPhyHelper()
    channel = ns.wifi.YansWifiChannelHelper.Default()
    phy.SetChannel(channel.Create())
    
    # MAC layer - Ad-hoc mode for P2P
    mac = ns.wifi.WifiMacHelper()
    mac.SetType("ns3::AdhocWifiMac")
    
    # Install WiFi devices
    devices = wifi.Install(phy, mac, nodes)
    
    # Install Internet stack
    stack = ns.internet.InternetStackHelper()
    stack.Install(nodes)
    
    # Assign IP addresses
    address = ns.internet.Ipv4AddressHelper()
    address.SetBase(ns.network.Ipv4Address("10.1.1.0"), 
                    ns.network.Ipv4Mask("255.255.255.0"))
    interfaces = address.Assign(devices)
    
    # Set up mobility (static positions in a grid)
    mobility = ns.mobility.MobilityHelper()
    
    grid_size = int(num_nodes ** 0.5) + 1
    spacing = 100.0
    
    mobility.SetPositionAllocator(
        "ns3::GridPositionAllocator",
        "MinX", ns.core.DoubleValue(0.0),
        "MinY", ns.core.DoubleValue(0.0),
        "DeltaX", ns.core.DoubleValue(spacing),
        "DeltaY", ns.core.DoubleValue(spacing),
        "GridWidth", ns.core.UintegerValue(grid_size),
        "LayoutType", ns.core.StringValue("RowFirst")
    )
    
    mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel")
    mobility.Install(nodes)
    
    # Calculate neighbors based on transmission range
    neighbor_lists = CalculateNeighbors(nodes, interfaces, transmission_range)
    
    # Install gossip applications on all nodes
    for i in range(num_nodes):
        app = GossipApplication()
        app.Setup(i, neighbor_lists[i])
        nodes.Get(i).AddApplication(app)
        app.SetStartTime(ns.core.Seconds(1.0))
        app.SetStopTime(ns.core.Seconds(simulation_time))
        
        neighbor_count = len(neighbor_lists[i])
        print(f"Node {i} configured with {neighbor_count} neighbors")
    
    print("\nStarting simulation...\n")
    
    # Schedule final statistics
    ns.core.Simulator.Schedule(
        ns.core.Seconds(simulation_time),
        PrintFinalStatistics
    )
    
    # Run simulation
    ns.core.Simulator.Stop(ns.core.Seconds(simulation_time + 1.0))
    ns.core.Simulator.Run()
    ns.core.Simulator.Destroy()
    
    print("\nSimulation completed successfully!")


if __name__ == '__main__':
    main()