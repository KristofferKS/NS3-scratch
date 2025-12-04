# Distributed Network System Simulation Description

## 1. Simulation Overview

This report presents a comprehensive ns-3 based simulation framework designed to evaluate various network technologies and protocols in distributed network systems. The simulation environment enables systematic comparison of different communication paradigms, physical layer technologies, transport protocols, and dissemination strategies suitable for distributed sensing and monitoring applications.

## 2. Simulation Architecture

### 2.1 Network Topology

The simulation implements a hybrid network topology consisting of:

- **Static Base Stations**: A set of stationary nodes (typically 5-10) positioned in strategic formations (e.g., pentagon or grid patterns) serving as infrastructure nodes
- **Mobile Nodes**: Multiple mobile entities (3-5 nodes) following predetermined waypoint paths through the network area
- **Network Coverage Area**: Simulation area of 100m × 100m with configurable dimensions
\section    
The base stations are positioned using geometric patterns (e.g., pentagon formation with 50m radius) to ensure overlapping coverage zones while mobile nodes traverse the network following custom waypoint paths at configurable speeds (typically 2-10 m/s).

### 2.2 Node Capabilities

Each node in the simulation is equipped with:
- Network interface card with configurable physical layer
- IPv4 network stack with UDP/TCP support
- Application layer for protocol implementation
- Mobility model (static for base stations, waypoint-based for mobile nodes)
- Position tracking and range detection capabilities

## 3. Technology Evaluation Parameters

### 3.1 Physical Layer Technologies

The simulation framework supports evaluation of multiple physical layer options:

#### WiFi (IEEE 802.11)
- **Standard**: 802.11a/b/g/n/ac configurable
- **Mode**: Ad-hoc (IBSS) or Infrastructure (AP-STA)
- **MAC Layer**: CSMA/CA with configurable parameters
- **Data Rates**: 1-600 Mbps depending on standard
- **Remote Station Manager**: AARF, Minstrel, or Constant Rate
- **Typical Range**: 50-100m depending on power and propagation model

#### LoRa Physical Layer (Planned Integration)
The framework is designed to accommodate LoRa physical layer with parameters:
- **Frequency Bands**: 433 MHz, 868 MHz (EU), 915 MHz (US)
- **Spreading Factor (SF)**: SF7-SF12 (configurable)
- **Bandwidth**: 125 kHz, 250 kHz, 500 kHz
- **Coding Rate**: 4/5, 4/6, 4/7, 4/8
- **Transmission Power**: 2-20 dBm
- **Expected Range**: 2-15 km depending on SF and environment

### 3.2 Transport Layer Protocols

#### UDP (User Datagram Protocol)
- **Characteristics**: Connectionless, low overhead, no reliability guarantees
- **Use Cases**: Real-time beaconing, periodic sensor data, gossip dissemination
- **Packet Size**: 64-1500 bytes
- **Advantages**: Low latency, broadcast support
- **Configuration**: Socket-based with broadcast capability enabled

#### TCP (Transmission Control Protocol)
- **Characteristics**: Connection-oriented, reliable, ordered delivery
- **Use Cases**: Critical data transfer, file transmission, control messages
- **Congestion Control**: NewReno, Cubic, or custom algorithms
- **Advantages**: Reliability, flow control, error recovery

### 3.3 Propagation and Channel Models

The simulation incorporates realistic wireless channel characteristics:

#### Propagation Loss Models
- **Log-Distance Path Loss**: Basic distance-dependent attenuation
  - Path loss exponent: 2-4 (configurable for environment)
  - Reference distance: 1m
  - Reference loss: Calculated from Friis equation

- **Interference Zones**: Custom areas with additional signal degradation
  - Zone definition: Center position, radius, attenuation (dB)
  - Multiple zones: Support for heterogeneous interference patterns
  - Typical attenuation: 10-25 dB in interference zones

#### Fading and Delay Models
- **Propagation Delay**: Constant speed model (speed of light)
- **Optional Fading**: Nakagami, Rayleigh for realistic channel effects

## 4. Communication Protocols

### 4.1 Beacon-Response Protocol

Mobile nodes periodically broadcast beacon messages containing:
- **Node Identifier**: Unique node ID
- **Position Information**: Current (x, y) coordinates
- **Timestamp**: Message generation time
- **Beacon Interval**: 1-5 seconds (configurable)

Static base stations listen for beacons and respond with acknowledgment messages, enabling:
- Mobile node detection and tracking
- Coverage area mapping
- Range-based connectivity analysis

### 4.2 Gossip/Epidemic Protocol

The simulation implements a probabilistic gossip protocol for message dissemination:

#### Protocol Parameters
- **Gossip Probability (p)**: 0.6-0.9 (percentage chance of forwarding)
- **Maximum Hop Count**: 5-10 hops (prevents infinite propagation)
- **Message Deduplication**: Track seen message IDs
- **Random Forwarding Delay**: 10-100ms (prevents collision storms)

#### Message Structure
```
- Message Type: GOSSIP
- Origin Node ID: Source of the message
- Message ID: Unique identifier (deduplication)
- Hop Count: Current hop number
- Timestamp: Creation time
- Position: Origin node location
- Payload: 64-128 bytes application data
```

#### Gossip Metrics
- **Coverage**: Percentage of nodes receiving message
- **Latency**: Time from origin to reception
- **Redundancy**: Number of duplicate receptions
- **Network Overhead**: Total packets vs. minimum required

## 5. Evaluation Plan and Parameters

### 5.1 Simulation Scenarios

The evaluation plan includes the following scenario matrix:

| Scenario | Physical Layer | Transport | Protocol | Mobile Nodes | Static Nodes | Duration |
|----------|---------------|-----------|----------|--------------|--------------|----------|
| S1       | WiFi          | UDP       | Beacon   | 3            | 5            | 60s      |
| S2       | WiFi          | UDP       | Gossip   | 3            | 5            | 60s      |
| S3       | WiFi          | UDP       | Both     | 3            | 5            | 60s      |
| S4       | LoRa          | UDP       | Beacon   | 3            | 5            | 120s     |
| S5       | LoRa          | UDP       | Gossip   | 3            | 5            | 120s     |
| S6       | WiFi+TCP      | TCP       | Gossip   | 3            | 5            | 60s      |

### 5.2 Variable Parameters for Sensitivity Analysis

#### Network Density
- Number of static nodes: 5, 10, 15
- Number of mobile nodes: 1, 3, 5, 10
- Network area: 100×100m, 200×200m, 500×500m

#### Mobility Parameters
- Mobile node speed: 2, 5, 10, 20 m/s
- Movement pattern: Random waypoint, predetermined paths, random walk
- Pause time at waypoints: 0, 5, 10 seconds

#### Physical Layer Parameters (WiFi)
- Transmission power: 10, 15, 20 dBm
- Data rate: 1, 6, 11, 54 Mbps
- MAC retry limit: 3, 7, 15 attempts

#### Physical Layer Parameters (LoRa)
- Spreading factor: SF7, SF9, SF12
- Bandwidth: 125 kHz, 250 kHz
- Coding rate: 4/5, 4/8
- Transmission power: 2, 10, 20 dBm

#### Gossip Protocol Parameters
- Gossip probability: 0.5, 0.7, 0.8, 0.9, 1.0
- Maximum hops: 3, 5, 7, 10, unlimited
- Forwarding delay: 10ms, 50ms, 100ms, 500ms

#### Beacon Protocol Parameters
- Beacon interval: 0.5s, 1s, 2s, 5s
- Beacon packet size: 64, 128, 256, 512 bytes
- Response timeout: 100ms, 500ms, 1s

#### Environmental Parameters
- Path loss exponent: 2.0, 2.5, 3.0, 4.0
- Interference zones: 0, 2, 5 zones
- Interference attenuation: 10dB, 20dB, 30dB

### 5.3 Performance Metrics

#### Connectivity Metrics
- **Detection Rate**: Percentage of time mobile nodes are detected by base stations
- **Detection Latency**: Time from beacon transmission to response reception
- **Coverage Area**: Geographic area with successful communication
- **Handover Events**: Number of transitions between base stations

#### Gossip Protocol Metrics
- **Message Delivery Ratio**: (Nodes receiving message) / (Total nodes)
- **Gossip Latency**: Time for message to reach all nodes
- **Redundancy Factor**: (Total transmissions) / (Nodes - 1)
- **Network Load**: Total bytes transmitted / simulation time

#### Network Performance Metrics
- **Packet Delivery Ratio (PDR)**: Successfully delivered / Total sent
- **End-to-End Delay**: Average packet transmission time
- **Throughput**: Successful data rate (Mbps or kbps)
- **Packet Loss Rate**: 1 - PDR

#### Energy Efficiency Metrics (Future Work)
- **Energy per Bit**: Joules per successfully delivered bit
- **Network Lifetime**: Time until first node depletes energy
- **Energy Balance**: Standard deviation of remaining energy

### 5.4 Data Collection and Analysis

#### Trace Files
- **PCAP Files**: Packet-level capture for detailed analysis
- **ASCII Traces**: Human-readable event logs
- **Custom Logs**: Application-level events and statistics

#### Visualization
- **NetAnim**: Real-time network animation showing node movement and packet flow
- **Position Traces**: CSV files with node trajectories
- **Connectivity Graphs**: Time-series plots of network topology changes

#### Statistical Analysis
- Mean, median, standard deviation for all metrics
- Confidence intervals (95%) with multiple simulation runs
- Comparative analysis across scenarios
- Correlation analysis between parameters and performance

## 6. Simulation Workflow

### 6.1 Initialization Phase
1. Create node containers (static and mobile)
2. Configure physical layer (WiFi or LoRa)
3. Install network stack (IPv4, UDP/TCP)
4. Position static nodes in formation
5. Configure mobile node waypoint paths
6. Assign IP addresses

### 6.2 Execution Phase
1. Start applications on all nodes
2. Mobile nodes begin waypoint traversal
3. Beacon transmission (if enabled)
4. Gossip message injection at specified times
5. Continuous packet trace collection
6. Event logging throughout simulation

### 6.3 Analysis Phase
1. Parse trace files and logs
2. Calculate performance metrics
3. Generate visualizations and plots
4. Compare scenarios and parameter variations
5. Statistical significance testing

## 7. Expected Outcomes

The simulation framework will provide insights into:

1. **Technology Comparison**: Quantitative comparison of WiFi vs. LoRa for distributed sensing
2. **Protocol Suitability**: Effectiveness of beacon-response vs. gossip for different scenarios
3. **Parameter Optimization**: Optimal values for gossip probability, beacon interval, etc.
4. **Scalability Analysis**: Performance degradation with increasing network size
5. **Mobility Impact**: Effect of node speed and movement patterns on connectivity
6. **Environmental Factors**: Impact of interference and path loss on communication
7. **Trade-offs**: Reliability vs. latency vs. energy efficiency relationships

## 8. Simulation Extensions

The modular framework allows for future extensions:

- **Advanced Routing**: AODV, OLSR, DSR routing protocols
- **Energy Models**: Battery simulation and energy-aware protocols
- **Security**: Authentication and encryption overhead analysis
- **Real Traces**: Integration of real-world mobility patterns
- **Machine Learning**: Adaptive parameter tuning based on network conditions
- **Hybrid Networks**: Combination of multiple physical layer technologies

## 9. Conclusion

This simulation framework provides a comprehensive testbed for evaluating distributed network systems under various configurations and environmental conditions. The systematic variation of physical layer technologies, transport protocols, and application-layer strategies enables data-driven decision-making for real-world deployments. The modular architecture ensures extensibility for future research directions while maintaining reproducibility of results.
