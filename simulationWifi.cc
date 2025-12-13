#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/applications-module.h"
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <deque>
#include <random>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include <nlohmann/json.hpp>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("GossipSimulation");

// DATA STRUCTURES
struct DensityReading {
    float value;
    uint32_t timestamp;  // milliseconds since simulation start
    
    DensityReading() : value(0.0f), timestamp(0) {}
    DensityReading(float v, uint32_t t) : value(v), timestamp(t) {}
};

struct GossipMessage {
    uint32_t originNodeId;   // Node that generated this data
    uint32_t senderNodeId;   // Node that sent this packet (for hop tracking)
    DensityReading data;
    
    GossipMessage() : originNodeId(0), senderNodeId(0) {}
    GossipMessage(uint32_t origin, uint32_t sender, const DensityReading& reading)
        : originNodeId(origin), senderNodeId(sender), data(reading) {}
};

// GOSSIP PACKET HEADER
class GossipHeader : public Header {
public:
    GossipHeader() : m_originId(0), m_senderId(0), m_density(0.0f), m_timestamp(0) {}
    
    GossipHeader(uint32_t origin, uint32_t sender, float density, uint32_t timestamp)
        : m_originId(origin), m_senderId(sender), m_density(density), m_timestamp(timestamp) {}
    
    static TypeId GetTypeId() {
        static TypeId tid = TypeId("GossipHeader")
            .SetParent<Header>()
            .AddConstructor<GossipHeader>();
        return tid;
    }
    
    virtual TypeId GetInstanceTypeId() const override { return GetTypeId(); }
    
    virtual uint32_t GetSerializedSize() const override {
        return 16;  // 4 uint32_t fields
    }
    
    virtual void Serialize(Buffer::Iterator start) const override {
        start.WriteHtonU32(m_originId);
        start.WriteHtonU32(m_senderId);
        start.WriteHtonU32(m_timestamp);
        start.WriteHtonU32(*reinterpret_cast<const uint32_t*>(&m_density));
    }
    
    virtual uint32_t Deserialize(Buffer::Iterator start) override {
        m_originId = start.ReadNtohU32();
        m_senderId = start.ReadNtohU32();
        m_timestamp = start.ReadNtohU32();
        uint32_t densityBits = start.ReadNtohU32();
        m_density = *reinterpret_cast<float*>(&densityBits);
        return GetSerializedSize();
    }
    
    virtual void Print(std::ostream& os) const override {
        os << "Origin=" << m_originId << " Sender=" << m_senderId 
           << " Density=" << m_density << " Time=" << m_timestamp;
    }
    
    // Getters
    uint32_t GetOriginId() const { return m_originId; }
    uint32_t GetSenderId() const { return m_senderId; }
    float GetDensity() const { return m_density; }
    uint32_t GetTimestamp() const { return m_timestamp; }

private:
    uint32_t m_originId;
    uint32_t m_senderId;
    float m_density;
    uint32_t m_timestamp;
};

NS_OBJECT_ENSURE_REGISTERED(GossipHeader);

// MESSAGE GENERATOR
class MessageGenerator : public Object {
public:
    static TypeId GetTypeId() {
        static TypeId tid = TypeId("MessageGenerator")
            .SetParent<Object>()
            .SetGroupName("Applications")
            .AddConstructor<MessageGenerator>();
        return tid;
    }
    
    MessageGenerator()
        : m_nodeId(0)
        , m_running(false)
        , m_minDensity(0.0f)
        , m_maxDensity(100.0f)
        , m_deltaDensity(10.0f)
        , m_currentDensity(50.0f)
    {
        m_rng = CreateObject<UniformRandomVariable>();
    }
    
    virtual ~MessageGenerator() { Stop(); }
    
    void Configure(uint32_t nodeId, float minDensity, float maxDensity, 
                   float deltaDensity, Time minInterval, Time maxInterval) {
        m_nodeId = nodeId;
        m_minDensity = minDensity;
        m_maxDensity = maxDensity;
        m_deltaDensity = deltaDensity;
        m_minInterval = minInterval;
        m_maxInterval = maxInterval;
        
        // Initialize to random value in range
        m_currentDensity = m_rng->GetValue(minDensity, maxDensity);
    }
    
    void SetCallback(Callback<void, const DensityReading&> callback) {
        m_callback = callback;
    }
    
    void Start() {
        if (m_running) return;
        m_running = true;
        
        // Schedule first generation with small random delay
        Time delay = MilliSeconds(m_rng->GetInteger(0, 500));
        m_event = Simulator::Schedule(delay, &MessageGenerator::Generate, this);
    }
    
    void Stop() {
        m_running = false;
        if (m_event.IsPending()) {
            Simulator::Cancel(m_event);
        }
    }

private:
    void Generate() {
        if (!m_running) return;
        
        // Generate new density value (random walk within delta)
        float minAllowed = std::max(m_minDensity, m_currentDensity - m_deltaDensity);
        float maxAllowed = std::min(m_maxDensity, m_currentDensity + m_deltaDensity);
        m_currentDensity = m_rng->GetValue(minAllowed, maxAllowed);
        
        // Create reading
        DensityReading reading(m_currentDensity, Simulator::Now().GetMilliSeconds());
        
        NS_LOG_INFO("Node " << m_nodeId << " generated density=" << m_currentDensity 
                    << " at t=" << reading.timestamp << "ms");
        
        // Invoke callback
        if (!m_callback.IsNull()) {
            m_callback(reading);
        }
        
        // Schedule next generation
        Time nextInterval = Seconds(m_rng->GetValue(m_minInterval.GetSeconds(), 
                                                     m_maxInterval.GetSeconds()));
        m_event = Simulator::Schedule(nextInterval, &MessageGenerator::Generate, this);
    }
    
    uint32_t m_nodeId;
    bool m_running;
    float m_minDensity;
    float m_maxDensity;
    float m_deltaDensity;
    float m_currentDensity;
    Time m_minInterval;
    Time m_maxInterval;
    EventId m_event;
    Callback<void, const DensityReading&> m_callback;
    Ptr<UniformRandomVariable> m_rng;
};

NS_OBJECT_ENSURE_REGISTERED(MessageGenerator);

// GOSSIP APPLICATION
class GossipApp : public Application {
public:
    static TypeId GetTypeId() {
        static TypeId tid = TypeId("GossipApp")
            .SetParent<Application>()
            .AddConstructor<GossipApp>();
        return tid;
    }
    
    GossipApp()
        : m_socket(nullptr)
        , m_port(9999)
        , m_commRange(50.0)
        , m_gossipRatio(0.5)
        , m_nNodes(0)
    {
        m_generator = CreateObject<MessageGenerator>();
    }
    
    virtual ~GossipApp() {}
    
    void Configure(uint16_t port, double commRange, double gossipRatio, uint32_t nNodes) {
        m_port = port;
        m_commRange = commRange;
        m_gossipRatio = gossipRatio;
        m_nNodes = nNodes;
    }
    
    void SetNodeIpMap(const std::map<uint32_t, Ipv4Address>& nodeToIp) {
        m_nodeToIp = nodeToIp;
    }
    
    void SetLogDirectory(const std::string& logDir) {
        m_logDir = logDir;
    }
    
    void ConfigureGenerator(float minDensity, float maxDensity, float deltaDensity,
                           Time minInterval, Time maxInterval) {
        m_genConfigured = true;
        m_genMinDensity = minDensity;
        m_genMaxDensity = maxDensity;
        m_genDeltaDensity = deltaDensity;
        m_genMinInterval = minInterval;
        m_genMaxInterval = maxInterval;
    }

protected:
    virtual void StartApplication() override {
        uint32_t nodeId = GetNode()->GetId();
        
        // Open log file
        std::string logFile = m_logDir + "node_" + std::to_string(nodeId) + ".csv";
        m_logStream.open(logFile);
        m_logStream << "event,time_s,origin_node,sender_node,density,timestamp_ms\n";
        
        // Initialize knowledge base
        for (uint32_t i = 0; i < m_nNodes; ++i) {
            m_knowledge[i] = std::deque<DensityReading>();
        }
        
        // Discover neighbors
        DiscoverNeighbors();
        
        // Setup socket
        m_socket = Socket::CreateSocket(GetNode(), UdpSocketFactory::GetTypeId());
        m_socket->Bind(InetSocketAddress(Ipv4Address::GetAny(), m_port));
        m_socket->SetRecvCallback(MakeCallback(&GossipApp::HandleReceive, this));
        
        // Configure and start generator
        if (m_genConfigured) {
            m_generator->Configure(nodeId, m_genMinDensity, m_genMaxDensity, 
                                  m_genDeltaDensity, m_genMinInterval, m_genMaxInterval);
            m_generator->SetCallback(MakeCallback(&GossipApp::OnNewReading, this));
            m_generator->Start();
        }
    }
    
    virtual void StopApplication() override {
        m_generator->Stop();
        
        if (m_socket) {
            m_socket->Close();
        }
        
        // Write final knowledge to file
        WriteFinalKnowledge();
        
        if (m_logStream.is_open()) {
            m_logStream.close();
        }
    }

private:
    void DiscoverNeighbors() {
        uint32_t myId = GetNode()->GetId();
        Vector myPos = GetNode()->GetObject<MobilityModel>()->GetPosition();
        
        m_neighbors.clear();
        
        for (uint32_t i = 0; i < m_nNodes; ++i) {
            if (i == myId) continue;
            
            Ptr<Node> otherNode = NodeList::GetNode(i);
            Vector otherPos = otherNode->GetObject<MobilityModel>()->GetPosition();
            
            double distance = CalculateDistance(myPos, otherPos);
            if (distance <= m_commRange) {
                m_neighbors.push_back(i);
            }
        }
        
        NS_LOG_INFO("Node " << myId << " discovered " << m_neighbors.size() << " neighbors");
    }
    
    void OnNewReading(const DensityReading& reading) {
        uint32_t myId = GetNode()->GetId();
        
        // Log generation
        m_logStream << "generate," << Simulator::Now().GetSeconds() << "," 
                    << myId << "," << myId << "," 
                    << reading.value << "," << reading.timestamp << "\n";
        
        // Update own knowledge
        UpdateKnowledge(myId, reading);
        
        // Gossip to neighbors
        GossipMessage msg(myId, myId, reading);
        GossipToNeighbors(msg, myId);  // Exclude self
    }
    
    void HandleReceive(Ptr<Socket> socket) {
        Ptr<Packet> packet;
        Address from;
        
        while ((packet = socket->RecvFrom(from))) {
            GossipHeader header;
            packet->RemoveHeader(header);
            
            uint32_t myId = GetNode()->GetId();
            uint32_t originId = header.GetOriginId();
            uint32_t senderId = header.GetSenderId();
            
            // Ignore self-messages
            if (originId == myId) continue;
            
            DensityReading reading(header.GetDensity(), header.GetTimestamp());
            
            // Log reception
            m_logStream << "receive," << Simulator::Now().GetSeconds() << ","
                       << originId << "," << senderId << ","
                       << reading.value << "," << reading.timestamp << "\n";
            
            // Check if this is new information
            if (IsNewerThan(originId, reading)) {
                // Update knowledge
                UpdateKnowledge(originId, reading);
                
                // Log update
                m_logStream << "update," << Simulator::Now().GetSeconds() << ","
                           << originId << "," << senderId << ","
                           << reading.value << "," << reading.timestamp << "\n";
                
                // Propagate further
                GossipMessage msg(originId, myId, reading);
                GossipToNeighbors(msg, senderId);  // Exclude sender
            } else {
                // Log drop
                m_logStream << "drop," << Simulator::Now().GetSeconds() << ","
                           << originId << "," << senderId << ","
                           << reading.value << "," << reading.timestamp << "\n";
            }
        }
    }
    
    bool IsNewerThan(uint32_t nodeId, const DensityReading& reading) {
        auto& knowledge = m_knowledge[nodeId];
        
        if (knowledge.empty()) return true;
        
        // Check if timestamp is newer than all existing entries
        for (const auto& entry : knowledge) {
            if (reading.timestamp <= entry.timestamp) {
                return false;
            }
        }
        return true;
    }
    
    void UpdateKnowledge(uint32_t nodeId, const DensityReading& reading) {
        auto& knowledge = m_knowledge[nodeId];
        knowledge.push_back(reading);
        
        // Keep only last 5 entries
        while (knowledge.size() > 5) {
            knowledge.pop_front();
        }
    }
    
    void GossipToNeighbors(const GossipMessage& msg, uint32_t excludeNode) {
        std::vector<uint32_t> candidates;
        
        for (uint32_t neighbor : m_neighbors) {
            if (neighbor != excludeNode && neighbor != msg.originNodeId) {
                candidates.push_back(neighbor);
            }
        }
        
        if (candidates.empty()) return;
        
        // Select subset based on gossip ratio
        size_t numToSend = std::max(1, (int)std::ceil(m_gossipRatio * candidates.size()));
        
        // Shuffle and select
        std::random_device rd;
        std::mt19937 gen(rd());
        std::shuffle(candidates.begin(), candidates.end(), gen);
        
        for (size_t i = 0; i < std::min(numToSend, candidates.size()); ++i) {
            SendTo(candidates[i], msg);
        }
    }
    
    void SendTo(uint32_t targetNodeId, const GossipMessage& msg) {
        Ipv4Address targetIp = m_nodeToIp[targetNodeId];
        
        GossipHeader header(msg.originNodeId, GetNode()->GetId(), 
                           msg.data.value, msg.data.timestamp);
        
        Ptr<Packet> packet = Create<Packet>();
        packet->AddHeader(header);
        
        m_socket->SendTo(packet, 0, InetSocketAddress(targetIp, m_port));
        
        // Log send
        m_logStream << "send," << Simulator::Now().GetSeconds() << ","
                   << msg.originNodeId << "," << targetNodeId << ","
                   << msg.data.value << "," << msg.data.timestamp << "\n";
    }
    
    void WriteFinalKnowledge() {
        std::ofstream knowledgeFile(m_logDir + "final_knowledge.csv", std::ios::app);
        
        // Write header if file is new
        knowledgeFile.seekp(0, std::ios::end);
        if (knowledgeFile.tellp() == 0) {
            knowledgeFile << "node_id";
            for (uint32_t i = 0; i < m_nNodes; ++i) {
                for (int j = 0; j < 5; ++j) {
                    knowledgeFile << ",n" << i << "_d" << j << ",n" << i << "_t" << j;
                }
            }
            knowledgeFile << "\n";
        }
        
        // Write this node's knowledge
        knowledgeFile << GetNode()->GetId();
        
        for (uint32_t i = 0; i < m_nNodes; ++i) {
            const auto& knowledge = m_knowledge[i];
            size_t count = 0;
            
            for (const auto& reading : knowledge) {
                knowledgeFile << "," << reading.value << "," << reading.timestamp;
                count++;
            }
            
            // Pad with empty values
            for (size_t j = count; j < 5; ++j) {
                knowledgeFile << ",,";
            }
        }
        
        knowledgeFile << "\n";
        knowledgeFile.close();
    }
    
    // Member variables
    Ptr<Socket> m_socket;
    uint16_t m_port;
    double m_commRange;
    double m_gossipRatio;
    uint32_t m_nNodes;
    std::string m_logDir;
    std::ofstream m_logStream;
    
    std::vector<uint32_t> m_neighbors;
    std::map<uint32_t, Ipv4Address> m_nodeToIp;
    std::map<uint32_t, std::deque<DensityReading>> m_knowledge;
    
    Ptr<MessageGenerator> m_generator;
    bool m_genConfigured = false;
    float m_genMinDensity;
    float m_genMaxDensity;
    float m_genDeltaDensity;
    Time m_genMinInterval;
    Time m_genMaxInterval;
};

NS_OBJECT_ENSURE_REGISTERED(GossipApp);

// GOSSIP APP HELPER
class GossipAppHelper {
public:
    GossipAppHelper(uint16_t port = 9999, double commRange = 50.0, 
                   double gossipRatio = 0.5, uint32_t nNodes = 50)
        : m_port(port)
        , m_commRange(commRange)
        , m_gossipRatio(gossipRatio)
        , m_nNodes(nNodes)
    {}
    
    void SetLogDirectory(const std::string& logDir) {
        m_logDir = logDir;
    }
    
    void SetNodeIpMap(const std::map<uint32_t, Ipv4Address>& nodeToIp) {
        m_nodeToIp = nodeToIp;
    }
    
    void SetMessageGenerator(float minDensity, float maxDensity, float deltaDensity,
                            Time minInterval, Time maxInterval) {
        m_useGenerator = true;
        m_minDensity = minDensity;
        m_maxDensity = maxDensity;
        m_deltaDensity = deltaDensity;
        m_minInterval = minInterval;
        m_maxInterval = maxInterval;
    }
    
    ApplicationContainer Install(NodeContainer nodes) {
        ApplicationContainer apps;
        for (auto it = nodes.Begin(); it != nodes.End(); ++it) {
            apps.Add(InstallPriv(*it));
        }
        return apps;
    }

private:
    Ptr<Application> InstallPriv(Ptr<Node> node) {
        Ptr<GossipApp> app = CreateObject<GossipApp>();
        app->Configure(m_port, m_commRange, m_gossipRatio, m_nNodes);
        app->SetNodeIpMap(m_nodeToIp);
        app->SetLogDirectory(m_logDir);
        
        if (m_useGenerator) {
            app->ConfigureGenerator(m_minDensity, m_maxDensity, m_deltaDensity,
                                   m_minInterval, m_maxInterval);
        }
        
        node->AddApplication(app);
        return app;
    }
    
    uint16_t m_port;
    double m_commRange;
    double m_gossipRatio;
    uint32_t m_nNodes;
    std::string m_logDir;
    std::map<uint32_t, Ipv4Address> m_nodeToIp;
    
    bool m_useGenerator = false;
    float m_minDensity;
    float m_maxDensity;
    float m_deltaDensity;
    Time m_minInterval;
    Time m_maxInterval;
};

// WIFI STATISTICS TRACKING - APPLICATION LEVEL
struct WifiStatistics {
    // Application-level metrics (calculated from logs at end)
    uint32_t appSends = 0;
    uint32_t appReceives = 0;
    uint32_t appUpdates = 0;
    uint32_t appDropsDuplicate = 0;
    
    // PHY-level metrics (for debugging only)
    uint32_t phyTxCount = 0;
    uint32_t phyRxCount = 0;
    uint32_t phyDropCount = 0;
    
    std::ofstream dropLog;
    
    void OpenLog(const std::string& filename) {
        dropLog.open(filename);
        dropLog << "time_s,packet_uid,drop_reason,packet_size\n";
    }
    
    void CloseLog() {
        if (dropLog.is_open()) {
            dropLog.close();
        }
    }
    
    void PrintSummary(const std::string& filename, const std::string& logDir) const {
        // Calculate from application logs
        std::cout << "\n=== Application-Level Statistics ===\n";
        std::cout << "App Sends: " << appSends << "\n";
        std::cout << "App Receives: " << appReceives << "\n";
        std::cout << "App Updates (new info): " << appUpdates << "\n";
        std::cout << "App Drops (duplicate): " << appDropsDuplicate << "\n";
        
        double deliveryRate = appSends > 0 ? 100.0 * appReceives / appSends : 0.0;
        double updateRate = appReceives > 0 ? 100.0 * appUpdates / appReceives : 0.0;
        
        std::cout << "Delivery Rate: " << std::fixed << std::setprecision(2) << deliveryRate << "%\n";
        std::cout << "Update Rate (of received): " << std::fixed << std::setprecision(2) << updateRate << "%\n";
        
        std::cout << "\n=== PHY-Level Statistics (for debugging) ===\n";
        std::cout << "PHY TX: " << phyTxCount << "\n";
        std::cout << "PHY RX: " << phyRxCount << "\n";
        std::cout << "PHY Drops (all nodes): " << phyDropCount << "\n";
        
        std::ofstream summary(filename);
        summary << "Application-Level Statistics\n";
        summary << "============================\n";
        summary << "App Sends: " << appSends << "\n";
        summary << "App Receives: " << appReceives << "\n";
        summary << "App Updates (new info): " << appUpdates << "\n";
        summary << "App Drops (duplicate): " << appDropsDuplicate << "\n";
        summary << "Delivery Rate: " << std::fixed << std::setprecision(2) << deliveryRate << "%\n";
        summary << "Update Rate (of received): " << std::fixed << std::setprecision(2) << updateRate << "%\n";
        summary << "\nPHY-Level Statistics (for debugging)\n";
        summary << "=====================================\n";
        summary << "PHY TX: " << phyTxCount << "\n";
        summary << "PHY RX: " << phyRxCount << "\n";
        summary << "PHY Drops (all nodes): " << phyDropCount << "\n";
        summary << "\nNote: PHY drops include all nodes that couldn't decode the signal.\n";
        summary << "This is expected behavior - only neighbors should successfully receive.\n";
        summary.close();
    }
};

WifiStatistics g_wifiStats;

// Simple PHY-level tracking (just for debugging, not meaningful metrics)
void PhyTxCallback(std::string context, Ptr<const Packet> packet, double txPowerW) {
    g_wifiStats.phyTxCount++;
}

void PhyRxOkCallback(std::string context, Ptr<const Packet> packet, uint16_t channelFreqMhz, 
                    WifiTxVector txVector, MpduInfo aMpdu, SignalNoiseDbm signalNoise, uint16_t staId) {
    g_wifiStats.phyRxCount++;
}

void PhyRxDropCallback(std::string context, Ptr<const Packet> packet, WifiPhyRxfailureReason reason) {
    g_wifiStats.phyDropCount++;
    // Don't log every drop - there are too many and they're expected
}

void MacTxDropCallback(std::string context, Ptr<const Packet> packet) {
    // MAC drops are more interesting - log these
    g_wifiStats.dropLog << Simulator::Now().GetSeconds() << ","
                       << packet->GetUid() << ",MAC_TX_DROP,"
                       << packet->GetSize() << "\n";
}

void ConnectWifiCallbacks(NodeContainer& nodes) {
    for (uint32_t i = 0; i < nodes.GetN(); ++i) {
        Ptr<Node> node = nodes.Get(i);
        
        for (uint32_t j = 0; j < node->GetNDevices(); ++j) {
            Ptr<WifiNetDevice> wifiDev = DynamicCast<WifiNetDevice>(node->GetDevice(j));
            if (wifiDev) {
                Ptr<WifiPhy> phy = wifiDev->GetPhy();
                
                // Use correct trace sources for NS-3.46
                phy->TraceConnect("PhyTxBegin", "Node" + std::to_string(i), 
                                 MakeCallback(&PhyTxCallback));
                // Use MonitorSnifferRx which has the right signature
                phy->TraceConnect("MonitorSnifferRx", "Node" + std::to_string(i), 
                                 MakeCallback(&PhyRxOkCallback));
                phy->TraceConnect("PhyRxDrop", "Node" + std::to_string(i), 
                                 MakeCallback(&PhyRxDropCallback));
                
                Ptr<WifiMac> mac = wifiDev->GetMac();
                mac->TraceConnect("MacTxDrop", "Node" + std::to_string(i), 
                                 MakeCallback(&MacTxDropCallback));
            }
        }
    }
}

// NODE PLACEMENT
void PlaceNodes(NodeContainer& nodes, uint32_t nNodes, double radius) {
    // Try to load from file first
    std::string filename = "scratch/placement_" + std::to_string(nNodes) + "_32.json";
    
    std::ifstream file(filename);
    std::string jsonStr;
    
    if (file.is_open()) {
        std::stringstream buffer;
        buffer << file.rdbuf();
        jsonStr = buffer.str();
        file.close();
    } else {
        // Generate using Python script
        std::string cmd = "python3 scratch/node_placement.py " + 
                         std::to_string(nNodes) + " 50 0 0 " + std::to_string(radius);
        
        FILE* pipe = popen(cmd.c_str(), "r");
        if (pipe) {
            char buffer[256];
            while (fgets(buffer, sizeof(buffer), pipe)) {
                jsonStr += buffer;
            }
            pclose(pipe);
        }
    }
    
    if (jsonStr.empty()) {
        NS_FATAL_ERROR("Failed to get node placement");
    }
    
    // Parse JSON
    auto placements = nlohmann::json::parse(jsonStr);
    
    // Sort by node ID
    std::vector<std::pair<uint32_t, Vector>> positions;
    for (auto& [key, value] : placements.items()) {
        uint32_t nodeId = std::stoi(key);
        double x = value[0];
        double y = value[1];
        positions.push_back({nodeId, Vector(x, y, 0.0)});
    }
    
    std::sort(positions.begin(), positions.end());
    
    // Apply positions (no scaling - use original coordinates)
    Ptr<ListPositionAllocator> allocator = CreateObject<ListPositionAllocator>();
    
    for (const auto& [id, pos] : positions) {
        allocator->Add(Vector(pos.x, pos.y, 0.0));
    }
    
    MobilityHelper mobility;
    mobility.SetPositionAllocator(allocator);
    mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    mobility.Install(nodes);
}

// NETWORK SETUP
NetDeviceContainer SetupWifi(NodeContainer& nodes, double pathLossExponent) {
    YansWifiChannelHelper channel;
    
    if (pathLossExponent > 0) {
        // Forest/festival propagation model
        // Log-distance with configurable exponent for outdoor with foliage/crowd
        // Reference loss ~46 dB at 1m (standard 2.4 GHz outdoor)
        channel.SetPropagationDelay("ns3::ConstantSpeedPropagationDelayModel");
        channel.AddPropagationLoss("ns3::LogDistancePropagationLossModel",
                                   "Exponent", DoubleValue(pathLossExponent),
                                   "ReferenceLoss", DoubleValue(46.0));
        channel.AddPropagationLoss("ns3::NakagamiPropagationLossModel",
                                   "m0", DoubleValue(1.0),   // Rayleigh fading (NLOS)
                                   "m1", DoubleValue(1.0),
                                   "m2", DoubleValue(1.0));
        NS_LOG_INFO("Using LogDistance + Nakagami propagation, exponent=" << pathLossExponent);
    } else {
        // Default YansWifiChannelHelper (LogDistance exponent=3, no fading)
        channel = YansWifiChannelHelper::Default();
        NS_LOG_INFO("Using default YansWifiChannelHelper propagation");
    }
    
    YansWifiPhyHelper phy;
    phy.SetChannel(channel.Create());
    
    // TX power: 20 dBm (~100 mW) - typical for outdoor ad-hoc scenario
    phy.Set("TxPowerStart", DoubleValue(20.0));
    phy.Set("TxPowerEnd", DoubleValue(20.0));
    
    // WiFi MAC
    WifiHelper wifi;
    wifi.SetStandard(WIFI_STANDARD_80211ac);
    
    WifiMacHelper mac;
    mac.SetType("ns3::AdhocWifiMac", "QosSupported", BooleanValue(true));
    
    // Enable RTS/CTS
    Config::SetDefault("ns3::WifiRemoteStationManager::RtsCtsThreshold", UintegerValue(100));
    
    return wifi.Install(phy, mac, nodes);
}
void PrintProgress(double totalSimTime)
{
    double now = Simulator::Now().GetSeconds();
    double percent = (now / totalSimTime) * 100.0;

    std::cout << "\rProgress: "
              << std::fixed << std::setprecision(1)
              << percent << "% (" << now << " / " << totalSimTime << "s)"
              << std::flush;

    if (percent < 100.0)
    {
        // Schedule next update
        Simulator::Schedule(Seconds(36.0), &PrintProgress, totalSimTime);
    }
}

// MAIN SIMULATION
int main(int argc, char* argv[]) {
    // Default parameters
    uint32_t nNodes = 50;
    double commRange = 70.0;
    double simTime = 3600.0;
    double gossipRatio = 0.1;
    float minDensity = 0.0f;
    float maxDensity = 100.0f;
    float deltaDensity = 15.0f;
    double minInterval = 1.0;
    double maxInterval = 3.0;
    double pathLossExponent = 3.3;  // 0 or negative = use default Friis
    std::string runName = "";      // Optional run name for log directory
    
    // Command line
    CommandLine cmd;
    cmd.AddValue("nNodes", "Number of nodes", nNodes);
    cmd.AddValue("commRange", "Communication range (m)", commRange);
    cmd.AddValue("simTime", "Simulation time (s)", simTime);
    cmd.AddValue("gossipRatio", "Fraction of neighbors to gossip to", gossipRatio);
    cmd.AddValue("minDensity", "Minimum density value", minDensity);
    cmd.AddValue("maxDensity", "Maximum density value", maxDensity);
    cmd.AddValue("deltaDensity", "Maximum density change per step", deltaDensity);
    cmd.AddValue("minInterval", "Minimum message interval (s)", minInterval);
    cmd.AddValue("maxInterval", "Maximum message interval (s)", maxInterval);
    cmd.AddValue("pathLossExp", "Path loss exponent (0 = default Friis)", pathLossExponent);
    cmd.AddValue("runName", "Custom name for this run (optional)", runName);
    cmd.Parse(argc, argv);
    
    // Enable logging
    LogComponentEnable("GossipSimulation", LOG_LEVEL_INFO);
    
    // Create log directory with meaningful name
    std::string logDir;
    if (!runName.empty()) {
        logDir = "simWifi/" + runName + "/";
    } else {
        // Auto-generate name from parameters
        std::stringstream ss;
        ss << "simWifi/";
        if (pathLossExponent > 0) {
            ss << "exp" << std::fixed << std::setprecision(1) << pathLossExponent;
        } else {
            ss << "default";
        }
        ss << "_gr" << std::fixed << std::setprecision(1) << gossipRatio;
        ss << "_t" << (int)simTime << "s/";
        logDir = ss.str();
    }
    std::filesystem::create_directories(logDir);
    
    std::cout << "Log directory: " << logDir << "\n";
    
    // Create nodes
    NodeContainer nodes;
    nodes.Create(nNodes);
    
    // Place nodes
    PlaceNodes(nodes, nNodes, 200.0);
    
    // Setup WiFi
    NetDeviceContainer devices = SetupWifi(nodes, pathLossExponent);
    
    // Install Internet stack
    InternetStackHelper internet;
    internet.Install(nodes);
    
    // Assign IP addresses
    Ipv4AddressHelper ipv4;
    ipv4.SetBase("10.1.1.0", "255.255.255.0");
    Ipv4InterfaceContainer interfaces = ipv4.Assign(devices);
    
    // Create node-to-IP map
    std::map<uint32_t, Ipv4Address> nodeToIp;
    for (uint32_t i = 0; i < nNodes; ++i) {
        nodeToIp[i] = interfaces.GetAddress(i);
    }
    
    // Setup WiFi statistics
    g_wifiStats.OpenLog(logDir + "wifi_drops.csv");
    ConnectWifiCallbacks(nodes);
    
    // Install gossip application
    GossipAppHelper gossipHelper(9999, commRange, gossipRatio, nNodes);
    gossipHelper.SetLogDirectory(logDir);
    gossipHelper.SetNodeIpMap(nodeToIp);
    gossipHelper.SetMessageGenerator(minDensity, maxDensity, deltaDensity,
                                     Seconds(minInterval), Seconds(maxInterval));
    
    ApplicationContainer apps = gossipHelper.Install(nodes);
    apps.Start(Seconds(1.0));
    apps.Stop(Seconds(simTime - 1.0));
    
    // Write configuration
    std::ofstream configFile(logDir + "config.txt");
    configFile << "Configuration\n=============\n";
    configFile << "Nodes: " << nNodes << "\n";
    configFile << "Comm Range: " << commRange << " m\n";
    configFile << "Simulation Time: " << simTime << " s\n";
    configFile << "Gossip Ratio: " << gossipRatio << "\n";
    configFile << "Density Range: [" << minDensity << ", " << maxDensity << "]\n";
    configFile << "Delta Density: " << deltaDensity << "\n";
    configFile << "Message Interval: [" << minInterval << ", " << maxInterval << "] s\n";
    configFile << "Path Loss Exponent: " << (pathLossExponent > 0 ? std::to_string(pathLossExponent) : "default (LogDistance exp=3)") << "\n";
    configFile << "TX Power: 20 dBm\n";
    configFile.close();
    
    // Run simulation
    Simulator::Stop(Seconds(simTime));
    std::cout << "Starting simulation...\n";
    Simulator::Schedule(Seconds(36.0), &PrintProgress, simTime);
    Simulator::Run();
    
    // Calculate application-level statistics from log files
    for (uint32_t i = 0; i < nNodes; ++i) {
        std::string logFile = logDir + "node_" + std::to_string(i) + ".csv";
        std::ifstream file(logFile);
        std::string line;
        std::getline(file, line);  // Skip header
        
        while (std::getline(file, line)) {
            if (line.substr(0, 5) == "send,") g_wifiStats.appSends++;
            else if (line.substr(0, 8) == "receive,") g_wifiStats.appReceives++;
            else if (line.substr(0, 7) == "update,") g_wifiStats.appUpdates++;
            else if (line.substr(0, 5) == "drop,") g_wifiStats.appDropsDuplicate++;
        }
    }
    
    // Print statistics
    g_wifiStats.CloseLog();
    g_wifiStats.PrintSummary(logDir + "wifi_summary.txt", logDir);
    
    Simulator::Destroy();
    
    std::cout << "Simulation complete. Results in: " << logDir << "\n";
    
    return 0;
}