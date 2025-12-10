#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/applications-module.h"
#include "ns3/netanim-module.h"
#include "ns3/buildings-module.h"
#include "ns3/header.h"
#include <cstdint>  // for uint32_t
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <nlohmann/json.hpp>
#include <algorithm>
#include <random>
#include <cmath>
#include <filesystem>
#include <chrono>
#include <sstream>
#include <iomanip> 
#include "ns3/point-to-point-module.h"

using namespace std;

using namespace ns3;
NS_LOG_COMPONENT_DEFINE("First");

map<string, NetDeviceContainer> p2pDevices;
map<pair<uint32_t, uint32_t>, pair<Ipv4Address, Ipv4Address>> linkIPs;  // Add this


/**
 * Structure representing a single density reading
 */
struct DensityReading
{
    float density;
    uint32_t timestamp;  // in milliseconds

    DensityReading() : density(0.0f), timestamp(0) {}
    DensityReading(float d, uint32_t t) : density(d), timestamp(t) {}
};

/**
 * Structure representing a complete message
 */
struct DensityMessage
{
    uint32_t nodeId;
    uint32_t lastNode;
    DensityReading reading;

    DensityMessage() : nodeId(0) {}
    DensityMessage(uint32_t id, float density, uint32_t timestamp)
        : nodeId(id), reading(density, timestamp) {}
    
    // Get serialized size
    uint32_t GetSerializedSize() const
    {
        return sizeof(nodeId) + sizeof(reading.density) + sizeof(reading.timestamp) + sizeof(lastNode);
    }
    
    // Serialize to buffer
    void Serialize(uint8_t* buffer) const
    {
        uint32_t offset = 0;
        
        memcpy(buffer + offset, &nodeId, sizeof(nodeId));
        offset += sizeof(nodeId);

        memcpy(buffer + offset, &lastNode, sizeof(lastNode));
        offset += sizeof(lastNode);
        
        memcpy(buffer + offset, &reading.density, sizeof(reading.density));
        offset += sizeof(reading.density);
        
        memcpy(buffer + offset, &reading.timestamp, sizeof(reading.timestamp));
    }
    
    // Deserialize from buffer
    static DensityMessage Deserialize(const uint8_t* buffer)
    {
        DensityMessage msg;
        uint32_t offset = 0;
        
        memcpy(&msg.nodeId, buffer + offset, sizeof(msg.nodeId));
        offset += sizeof(msg.nodeId);

        memcpy(&msg.lastNode, buffer + offset,sizeof(msg.lastNode));
        offset += sizeof(msg.lastNode);

        memcpy(&msg.reading.density, buffer + offset, sizeof(msg.reading.density));
        offset += sizeof(msg.reading.density);
        
        memcpy(&msg.reading.timestamp, buffer + offset, sizeof(msg.reading.timestamp));
        
        return msg;
    }
};

/**
 * Callback type for when a new message is generated
 */
typedef Callback<void, const DensityMessage&> MessageGeneratedCallback;

/**
 * MessageGenerator class
 */
class MessageGenerator : public Object
{
public:
    static TypeId GetTypeId(void)
    {
        static TypeId tid = TypeId("MessageGenerator")
            .SetParent<Object>()
            .SetGroupName("Applications")
            .AddConstructor<MessageGenerator>();
        return tid;
    }
    
    MessageGenerator()
        : m_nodeId(0),
          m_lastNode(0),
          m_minDensity(0.0f),
          m_maxDensity(100.0f),
          m_deltaDensity(10.0f),
          m_minInterval(Seconds(1.0)),
          m_maxInterval(Seconds(5.0)),
          m_currentDensity(0.0f),
          m_running(false)
    {
        m_randomInterval = CreateObject<UniformRandomVariable>();
        m_randomDensity = CreateObject<UniformRandomVariable>();
    }

    virtual ~MessageGenerator()
    {
        Stop();
    }

    void Setup(uint32_t nodeId,
               uint32_t lastNode,
               float minDensity,
               float maxDensity,
               float deltaDensity,
               Time minInterval,
               Time maxInterval)
    {
        NS_ASSERT(minDensity <= maxDensity);
        NS_ASSERT(deltaDensity > 0.0f);
        NS_ASSERT(minInterval <= maxInterval);

        m_nodeId = nodeId;
        m_lastNode = lastNode;
        m_minDensity = minDensity;
        m_maxDensity = maxDensity;
        m_deltaDensity = deltaDensity;
        m_minInterval = minInterval;
        m_maxInterval = maxInterval;

        // Initialize current density to a random value in range
        m_randomDensity->SetAttribute("Min", DoubleValue(minDensity));
        m_randomDensity->SetAttribute("Max", DoubleValue(maxDensity));
        m_currentDensity = m_randomDensity->GetValue();

        NS_LOG_INFO("MessageGenerator setup for node " << nodeId 
                    << " - Density range: [" << minDensity << ", " << maxDensity 
                    << "], Delta: " << deltaDensity
                    << ", Interval range: [" << minInterval.GetSeconds() 
                    << "s, " << maxInterval.GetSeconds() << "s]");
    }

    void Start()
    {
        if (m_running)
        {
            NS_LOG_WARN("MessageGenerator already running for node " << m_nodeId);
            return;
        }

        m_running = true;
        NS_LOG_INFO("MessageGenerator started for node " << m_nodeId);
        
        // Schedule first message with a small random delay to avoid collisions
        Time initialDelay = MilliSeconds(m_randomInterval->GetInteger(0, 500));
        m_generateEvent = Simulator::Schedule(initialDelay, &MessageGenerator::GenerateMessage, this);
    }

    void Stop()
    {
        if (!m_running){return;}

        m_running = false;
        
        if (m_generateEvent.IsPending())
        {
            Simulator::Cancel(m_generateEvent);
        }

        NS_LOG_INFO("MessageGenerator stopped for node " << m_nodeId);
    }

    void SetMessageGeneratedCallback(MessageGeneratedCallback callback)
    {
        m_messageCallback = callback;
    }

    DensityMessage GetLastMessage() const { return m_lastMessage; }

private:
    void GenerateMessage()
    {
        if (!m_running){return;}

        // Generate new density value
        float newDensity = GenerateNewDensity();
        m_currentDensity = newDensity;

        // Get current timestamp in milliseconds
        uint32_t timestamp = static_cast<uint32_t>(Simulator::Now().GetMilliSeconds());

        // Create message
        m_lastMessage = DensityMessage(m_nodeId, newDensity, timestamp);

        NS_LOG_INFO("Node " << m_nodeId << " generated message: density=" 
                    << newDensity << ", timestamp=" << timestamp << "ms");

        // Invoke callback if set
        if (!m_messageCallback.IsNull())
        {
            m_messageCallback(m_lastMessage);
        }

        // Schedule next message
        ScheduleNextMessage();
    }

    void ScheduleNextMessage()
    {
        if (!m_running)
        {
            return;
        }

        Time nextInterval = GetRandomInterval();
        m_generateEvent = Simulator::Schedule(nextInterval, &MessageGenerator::GenerateMessage, this);

        NS_LOG_DEBUG("Node " << m_nodeId << " scheduled next message in " 
                     << nextInterval.GetSeconds() << "s");
    }

    Time GetRandomInterval()
    {
        double minSeconds = m_minInterval.GetSeconds();
        double maxSeconds = m_maxInterval.GetSeconds();
        
        m_randomInterval->SetAttribute("Min", DoubleValue(minSeconds));
        m_randomInterval->SetAttribute("Max", DoubleValue(maxSeconds));
        
        double randomSeconds = m_randomInterval->GetValue();
        return Seconds(randomSeconds);
    }

    float GenerateNewDensity()
    {
        // Calculate the range for new density based on delta constraint
        float minAllowed = max(m_minDensity, m_currentDensity - m_deltaDensity);
        float maxAllowed = min(m_maxDensity, m_currentDensity + m_deltaDensity);

        // Generate random value within allowed range
        m_randomDensity->SetAttribute("Min", DoubleValue(minAllowed));
        m_randomDensity->SetAttribute("Max", DoubleValue(maxAllowed));
        
        float newDensity = m_randomDensity->GetValue();

        NS_LOG_DEBUG("Node " << m_nodeId << " density: " << m_currentDensity 
                     << " -> " << newDensity 
                     << " (allowed range: [" << minAllowed << ", " << maxAllowed << "])");

        return newDensity;
    }

    // Configuration
    uint32_t m_nodeId;
    uint32_t m_lastNode;
    float m_minDensity;
    float m_maxDensity;
    float m_deltaDensity;
    Time m_minInterval;
    Time m_maxInterval;

    // State
    float m_currentDensity;
    DensityMessage m_lastMessage;
    EventId m_generateEvent;
    bool m_running;

    // Callback
    MessageGeneratedCallback m_messageCallback;

    // Random number generation
    Ptr<UniformRandomVariable> m_randomInterval;
    Ptr<UniformRandomVariable> m_randomDensity;
};

NS_OBJECT_ENSURE_REGISTERED(MessageGenerator);

class GossipHeader : public Header
{
public:
    GossipHeader() : m_nodeId(0), m_lastNode(0), m_density(0.0f), m_timestamp(0) {}
    GossipHeader(uint32_t nodeId, uint32_t lastNode, float density, uint32_t timestamp)
        : m_nodeId(nodeId), m_lastNode(lastNode), m_density(density), m_timestamp(timestamp) {}

    void Set(uint32_t nodeId, uint32_t lastNode, float density, uint32_t timestamp)
    {
        m_nodeId = nodeId;
        m_lastNode = lastNode;
        m_density = density;
        m_timestamp = timestamp;
    }

    static TypeId GetTypeId(void)
    {
        static TypeId tid = TypeId("GossipHeader")
            .SetParent<Header>()
            .AddConstructor<GossipHeader>();
        return tid;
    }
    virtual TypeId GetInstanceTypeId(void) const override { return GetTypeId(); }

    virtual void Serialize(Buffer::Iterator start) const override
    {
        start.WriteU32(m_nodeId);
        start.WriteU32(m_lastNode);
        start.WriteHtonU32(m_timestamp);
        start.WriteU32(static_cast<uint32_t>(m_density * 10000)); // store as int
    }

    virtual uint32_t Deserialize(Buffer::Iterator start) override
    {
        m_nodeId = start.ReadU32();
        m_lastNode = start.ReadU32();
        m_timestamp = start.ReadNtohU32();
        m_density = static_cast<float>(start.ReadU32()) / 10000.0f;
        return GetSerializedSize();
    }

    virtual uint32_t GetSerializedSize(void) const override
    {
        return 4 + 4 + 4 + 4;
    }

    virtual void Print(ostream &os) const override
    {
        os << "NodeId=" << m_nodeId << ", LastNode=" << m_lastNode << ", Density=" << m_density << ", Timestamp=" << m_timestamp;
    }

    uint32_t GetNodeId() const { return m_nodeId; }
    uint32_t GetLastNode() const { return m_lastNode; }
    float GetDensity() const { return m_density; }
    uint32_t GetTimestamp() const { return m_timestamp; }

private:
    uint32_t m_nodeId;
    uint32_t m_lastNode;
    float m_density;
    uint32_t m_timestamp;
};
NS_OBJECT_ENSURE_REGISTERED(GossipHeader);

/**
 * GossipApp - Gossip protocol application with message generation
 */
class GossipApp : public Application
{
public:
    static TypeId GetTypeId(void);
    GossipApp();
    virtual ~GossipApp();

    void Setup(double range, uint16_t port, double percentage, uint32_t nCameraNodes);
    void SetupLogDir(string logDir);
    void SetNodeToIP(const map<uint32_t, Ipv4Address>& nodeToIP);
    
    // Message generator configuration
    void SetupMessageGenerator(float minDensity, float maxDensity, 
                               float deltaDensity, Time minInterval, Time maxInterval);

protected:
    virtual void StartApplication(void);
    virtual void StopApplication(void);

private:
    void SendMessage(uint32_t peerIndex, uint32_t nodeId, uint32_t lastNode, float density, uint32_t timestamp);
    void DiscoverPeers(void);
    void OnMessageGenerated(const DensityMessage& message);
    void ReceiveFromDevice(Ptr<Socket> socket);

    Ptr<Socket> m_socket;
    double m_commRange;
    uint16_t m_port;
    vector<uint32_t> m_peers;
    map<uint32_t, Ipv4Address> m_nodeToIP;
    ofstream m_logFile;
    double m_percentage;
    uint32_t m_nCameraNodes;
    
    // Message generator
    Ptr<MessageGenerator> m_messageGenerator;
    
    // Message generator parameters (stored until StartApplication)
    bool m_messageGeneratorConfigured;
    float m_minDensity;
    float m_maxDensity;
    float m_deltaDensity;
    Time m_minInterval;
    Time m_maxInterval;

    string m_logDir;

    // Density knowledge for gossip protocol
    map<uint32_t, deque<pair<float, uint32_t>>> m_densityKnowledge;  // nodeId -> deque of (density, timestamp) pairs, max 5
};

TypeId
GossipApp::GetTypeId(void)
{
    static TypeId tid = TypeId("GossipApp")
        .SetParent<Application>()
        .AddConstructor<GossipApp>();
    return tid;
}

GossipApp::GossipApp()
    : m_socket(0),
      m_commRange(50.0),
      m_port(9),
      m_messageGeneratorConfigured(false),
      m_minDensity(0.0f),
      m_maxDensity(100.0f),
      m_deltaDensity(10.0f),
      m_minInterval(Seconds(1.0)),
      m_maxInterval(Seconds(5.0))
{
    m_messageGenerator = CreateObject<MessageGenerator>();
}

GossipApp::~GossipApp()
{
    m_socket = 0;
}

void GossipApp::Setup(double range, uint16_t port, double percentage, uint32_t nCameraNodes)
{
    m_commRange = range;
    m_port = port;
    m_percentage = percentage;
    m_nCameraNodes = nCameraNodes;
}

void GossipApp::SetupLogDir(string logDir)
{
    m_logDir = logDir;
}

void GossipApp::SetNodeToIP(const map<uint32_t, Ipv4Address>& nodeToIP)
{
    m_nodeToIP = nodeToIP;
}

void
GossipApp::SetupMessageGenerator(float minDensity, float maxDensity, 
                                 float deltaDensity, Time minInterval, Time maxInterval)
{
    // Store parameters - actual setup happens in StartApplication when node is available
    m_messageGeneratorConfigured = true;
    m_minDensity = minDensity;
    m_maxDensity = maxDensity;
    m_deltaDensity = deltaDensity;
    m_minInterval = minInterval;
    m_maxInterval = maxInterval;
}

void GossipApp::DiscoverPeers(void)
{
    uint32_t myId = GetNode()->GetId();
    m_peers.clear();
    for (const auto& [key, device] : p2pDevices) {
        size_t dash = key.find('-');
        uint32_t n1 = stoul(key.substr(0, dash));
        uint32_t n2 = stoul(key.substr(dash + 1));
        if (n1 == myId) m_peers.push_back(n2);
        else if (n2 == myId) m_peers.push_back(n1);
    }
}

void
GossipApp::StartApplication(void)
{
    // Open log file for this node
    string fname = m_logDir + "gossip_log_node" + to_string(GetNode()->GetId()) + ".csv";
    m_logFile.open(fname, ios::out);
    m_logFile << "event,time,node_from,node_to,density,timestamp_ms" << endl;

    // Create UDP sockets
    m_socket = Socket::CreateSocket(GetNode(), UdpSocketFactory::GetTypeId());
    m_socket->Bind(InetSocketAddress(Ipv4Address::GetAny(), m_port));
    m_socket->SetRecvCallback(MakeCallback(&GossipApp::ReceiveFromDevice, this));

    // Discover peers
    DiscoverPeers();

    // Initialize density knowledge for all nodes
    for (uint32_t i = 0; i < m_nCameraNodes; ++i)
    {
        m_densityKnowledge[i] = deque<pair<float, uint32_t>>();  // Empty deque for each node
    }

    // Setup message generator now that node is available
    if (m_messageGeneratorConfigured)
    {
        uint32_t nodeId = GetNode()->GetId();
        uint32_t m_lastNode = nodeId; // Assuming lastNode is the same as nodeId for initial setup
        m_messageGenerator->Setup(nodeId,
                                    m_lastNode, 
                                    m_minDensity, 
                                    m_maxDensity, 
                                    m_deltaDensity,
                                    m_minInterval, 
                                    m_maxInterval);
        
        // Set callback for when messages are generated
        m_messageGenerator->SetMessageGeneratedCallback(
            MakeCallback(&GossipApp::OnMessageGenerated, this));
    }

    // Start message generator
    m_messageGenerator->Start();
}

void GossipApp::ReceiveFromDevice(Ptr<Socket> socket)
{
    Ptr<Packet> packet;
    Address from;
    while ((packet = socket->RecvFrom(from))) {
        uint32_t myId = GetNode()->GetId();
        double now = Simulator::Now().GetSeconds();
        Ptr<Packet> copy = packet->Copy();
        GossipHeader header;
        copy->RemoveHeader(header);
        uint32_t senderId = header.GetNodeId();
        uint32_t lastNode = header.GetLastNode();

        // Ignore messages from self
        if (senderId == myId) {
            NS_LOG_DEBUG("Node " << myId << " ignoring self-message");
            return;
        }

        // Log recv event to CSV
        m_logFile << "recv," << now << "," << senderId << "," << myId << ","
                  << header.GetDensity() << "," << header.GetTimestamp() << endl;

        NS_LOG_INFO("Node " << myId << " received message from node " << senderId
                    << " (density=" << header.GetDensity()
                    << ", timestamp=" << header.GetTimestamp() << " ms) at " << now << "s");

        // Update density knowledge for the sender node
        if (m_densityKnowledge.find(senderId) != m_densityKnowledge.end()) {
            // Check if the new timestamp is newer than existing ones
            uint32_t maxExistingTimestamp = 0;
            for (const auto& entry : m_densityKnowledge[senderId]) {
                if (entry.second > maxExistingTimestamp) {
                    maxExistingTimestamp = entry.second;
                }
            }

            if (header.GetTimestamp() > maxExistingTimestamp) {
                // Update knowledge
                m_densityKnowledge[senderId].push_back({header.GetDensity(), header.GetTimestamp()});
                // Keep only the last 5 entries
                if (m_densityKnowledge[senderId].size() > 5) {
                    m_densityKnowledge[senderId].pop_front();
                }
                NS_LOG_DEBUG("Node " << myId << " updated knowledge for node " << senderId 
                             << ": density=" << header.GetDensity() << ", timestamp=" << header.GetTimestamp());

                // Log update event to CSV
                m_logFile << "update," << now << "," << senderId << "," << myId << ","
                          << header.GetDensity() << "," << header.GetTimestamp() << endl;

                NS_LOG_INFO("Node " << myId << " got a message from " << lastNode);
                // Prepare list of peers excluding sender and self
                vector<uint32_t> shuffledPeers;
                shuffledPeers.reserve(m_peers.size());
                for (auto peer : m_peers) {
                    if (peer != lastNode && peer != myId) {
                        shuffledPeers.push_back(peer);
                    }
                }

                if (shuffledPeers.empty()) {
                    NS_LOG_DEBUG("Node " << myId << " has no peers to propagate to after excluding sender and self");
                    return;
                }

                // Propagate the update to a percentage of neighbors (excluding the sender, and self)
                double percentage = m_percentage;
                size_t nPeers = shuffledPeers.size();
                size_t nToSend = static_cast<size_t>(ceil(percentage * nPeers));
                if (nToSend == 0 && nPeers > 0) nToSend = 1;
                
                // Shuffle peers
                random_device rd;
                mt19937 g(rd());
                shuffle(shuffledPeers.begin(), shuffledPeers.end(), g);

                size_t limit = min(nToSend, shuffledPeers.size());
                size_t sentCount = 0;

                for (size_t idx = 0; idx < limit && sentCount < limit; ++idx) {
                    uint32_t peerIndex = shuffledPeers[idx];
                    SendMessage(peerIndex, header.GetNodeId(), myId, header.GetDensity(), header.GetTimestamp());
                }
            }
            else {
                NS_LOG_DEBUG("Node " << myId << " received outdated message from node " << senderId 
                            << " (timestamp=" << header.GetTimestamp() << " <= max=" << maxExistingTimestamp << "), dropping");

                // Log drop event to CSV
                m_logFile << "drop," << now << "," << senderId << "," << myId << ","
                        << header.GetDensity() << "," << header.GetTimestamp() << endl;
            }
        }
    }
}

void GossipApp::StopApplication(void)
{
    // Stop message generator
    m_messageGenerator->Stop();

    if (m_socket)
    {
        m_socket->Close();
    }

    // Write final knowledge to CSV with expanded tabular headers and columns
    ofstream knowledgeFile(m_logDir + "final_knowledge.csv", ios::app);
    if (knowledgeFile.is_open())
    {
        // Add header if file is empty (first write)
        knowledgeFile.seekp(0, ios::end);
        if (knowledgeFile.tellp() == 0)
        {
            knowledgeFile << "index";
            for (uint32_t node = 0; node < m_nCameraNodes; ++node)
            {
                for (size_t ent = 0; ent < 5; ++ent)  // Assuming max 5 entries per node
                {
                    knowledgeFile << ",node" << node << "density" << ent << ",node" << node << "time" << ent;
                }
            }
            knowledgeFile << endl;
        }

        // Write the node ID (index)
        knowledgeFile << GetNode()->GetId();

        // Write densities and timestamps for each node, padding with empty values
        for (uint32_t node = 0; node < m_nCameraNodes; ++node)
        {
            auto it = m_densityKnowledge.find(node);
            size_t entryCount = 0;
            if (it != m_densityKnowledge.end())
            {
                for (const auto& entry : it->second)
                {
                    knowledgeFile << "," << entry.first << "," << entry.second;
                    ++entryCount;
                }
            }
            // Pad with empty values for missing entries (up to 5)
            for (size_t ent = entryCount; ent < 5; ++ent)
            {
                knowledgeFile << ",,";
            }
        }
        knowledgeFile << endl;
        knowledgeFile.close();
    }
    else
    {
        NS_LOG_ERROR("Failed to open final_knowledge.csv for writing");
    }

    if (m_logFile.is_open())
    {
        m_logFile.close();
    }
}

void GossipApp::OnMessageGenerated(const DensityMessage& message)
{
    NS_LOG_INFO("Node " << GetNode()->GetId() 
                << " received new message from generator: density=" 
                << message.reading.density << ", timestamp=" 
                << message.reading.timestamp << "ms");

    double now = Simulator::Now().GetSeconds();

    // Log generated message to CSV
    m_logFile << "gen," << now << "," << GetNode()->GetId() << "," << GetNode()->GetId() << ","
              << message.reading.density << "," << message.reading.timestamp << endl;

    double percentage = m_percentage;
    size_t nPeers = m_peers.size();
    size_t nToSend = static_cast<size_t>(ceil(percentage * nPeers));
    if (nToSend == 0 && nPeers > 0) nToSend = 1;

    // Shuffle peers
    vector<uint32_t> shuffledPeers = m_peers;
    random_device rd;
    mt19937 g(rd());
    shuffle(shuffledPeers.begin(), shuffledPeers.end(), g);

    uint32_t lastNode = GetNode()->GetId();  // Define lastNode

    for (size_t idx = 0; idx < nToSend; ++idx)
    {
        uint32_t peerIndex = shuffledPeers[idx];
        SendMessage(peerIndex, message.nodeId, lastNode, message.reading.density, message.reading.timestamp);
    }
}

void GossipApp::SendMessage(uint32_t peerIndex, uint32_t nodeId, uint32_t lastNode, float density, uint32_t timestamp)
{
    if (m_peers.empty())
    {
        NS_LOG_WARN("Node " << GetNode()->GetId() << " has no peers to send to");
        return;
    }

    // Create packet and add custom header
    Ptr<Packet> packet = Create<Packet>();
    GossipHeader header(nodeId, lastNode, density, timestamp);
    packet->AddHeader(header);

    // Get the correct IP for this specific link
    uint32_t myId = GetNode()->GetId();
    uint32_t minId = min(myId, peerIndex);
    uint32_t maxId = max(myId, peerIndex);
    auto it = linkIPs.find({minId, maxId});
    if (it == linkIPs.end())
    {
        NS_LOG_ERROR("No IP found for link between " << myId << " and " << peerIndex << ", skipping send");
        return;
    }
    Ipv4Address peerIP = (myId == minId) ? it->second.second : it->second.first;

    m_socket->SendTo(packet, 0, InetSocketAddress(peerIP, m_port));

    // Log send event to CSV
    double now = Simulator::Now().GetSeconds();
    m_logFile << "send," << now << "," << GetNode()->GetId() << "," << peerIndex << ","
              << density << "," << timestamp << endl;

    NS_LOG_INFO("Packet sent to node " << peerIndex << " at IP " << peerIP);
}


/**
 * GossipAppHelper
 */
class GossipAppHelper
{
public:
    GossipAppHelper(double range = 50.0, uint16_t port = 9, string logDir = "", double percentage = 0.5, uint32_t nCameraNodes = 50);
    
    // Set message generator parameters
    void SetMessageGenerator(float minDensity, 
                                float maxDensity, 
                                float deltaDensity, 
                                Time minInterval, 
                                Time maxInterval);
    
    void SetNodeToIP(const map<uint32_t, Ipv4Address>& nodeToIP);
    
    ApplicationContainer Install(NodeContainer c);
    ApplicationContainer Install(Ptr<Node> node);

private:
    double m_percentage;
    double m_range;
    uint16_t m_port;
    string m_logDir;
    uint32_t m_nCameraNodes;
    map<uint32_t, Ipv4Address> m_nodeToIP;
    
    // Message generator parameters
    float m_minDensity;
    float m_maxDensity;
    float m_deltaDensity;
    Time m_minInterval;
    Time m_maxInterval;
    bool m_useMessageGenerator;

    void SendMessage(uint32_t peerIndex, uint32_t nodeId, uint32_t lastNode, float density, uint32_t timestamp);
};

GossipAppHelper::GossipAppHelper(double range, uint16_t port, string logDir, double percentage, uint32_t nCameraNodes)
    : m_percentage(percentage),
      m_range(range), 
      m_port(port),
      m_logDir(logDir),
      m_nCameraNodes(nCameraNodes),
      m_minDensity(0.0f),
      m_maxDensity(100.0f),
      m_deltaDensity(10.0f),
      m_minInterval(Seconds(1.0)),
      m_maxInterval(Seconds(5.0)),
      m_useMessageGenerator(true)
{
}

void GossipAppHelper::SetMessageGenerator(float minDensity, 
                                            float maxDensity, 
                                            float deltaDensity, 
                                            Time minInterval, 
                                            Time maxInterval)
{
    m_minDensity = minDensity;
    m_maxDensity = maxDensity;
    m_deltaDensity = deltaDensity;
    m_minInterval = minInterval;
    m_maxInterval = maxInterval;
    m_useMessageGenerator = true;
}

void GossipAppHelper::SetNodeToIP(const map<uint32_t, Ipv4Address>& nodeToIP)
{
    m_nodeToIP = nodeToIP;
}

ApplicationContainer
GossipAppHelper::Install(NodeContainer c)
{
    ApplicationContainer apps;
    for (NodeContainer::Iterator i = c.Begin(); i != c.End(); ++i)
    {
        apps.Add(Install(*i));
    }
    return apps;
}

ApplicationContainer
GossipAppHelper::Install(Ptr<Node> node)
{
    Ptr<GossipApp> app = CreateObject<GossipApp>();
    app->Setup(m_range, m_port, m_percentage, m_nCameraNodes);
    app->SetupLogDir(m_logDir);
    app->SetNodeToIP(m_nodeToIP);
    
    if (m_useMessageGenerator)
    {
        app->SetupMessageGenerator(m_minDensity, 
                                    m_maxDensity, 
                                    m_deltaDensity,
                                    m_minInterval, 
                                    m_maxInterval);
    }
    
    node->AddApplication(app);
    return ApplicationContainer(app);
}


NodeContainer CreateCameraNodes(uint32_t nNodes);
void SetCameraMobility(uint8_t node, NodeContainer& camera_nodes);
NodeContainer CreateMobileNodes();
void SetMobileNodeMobility(uint8_t node, NodeContainer& mobile_nodes);

NodeContainer CreateCameraNodes(uint32_t nNodes) {
    // Log the creation of camera nodes
    NS_LOG_INFO("Creating camera nodes");
    uint32_t node = nNodes;

    NodeContainer camera_nodes;
    camera_nodes.Create(node);

    NS_LOG_INFO("Setting camera node mobility");

    SetCameraMobility(node, camera_nodes);
    return camera_nodes;
}

void SetCameraMobility(uint8_t node, NodeContainer& camera_nodes) {
    NS_LOG_INFO("Setting camera node mobility");
    double radius = 200.0;         // radius of disc
    double wifiRange = 50.0;      // communication range
    double centerX = 0.0;
    double centerY = 0.0;

    NS_LOG_INFO("Placing camera nodes with coverage constraints");
    
    // Path to the directory
    const char* node_file = "scratch/placement.json";

    // Structure which would store the metadata
    struct stat sb;

    string placement;

    if (stat(node_file, &sb) == 0){

        // Read from the text file
        ifstream MyReadFile(node_file);

        // Use a while loop together with the getline() function to read the file line by line
        while (getline (MyReadFile, placement)) {
        // Output the text from the file
        cout << placement;
        }

        // Close the file
        MyReadFile.close();

    }

    else{
        string command = "python3 scratch/node_placement.py " + to_string(node) + " " + to_string(wifiRange) + " " + to_string(centerX) + " " + to_string(centerY) + " " + to_string(radius);
        // Read positions from file
        char buffer[128];
        FILE* pipe = popen(command.c_str(), "r");
        if (!pipe) {
            cerr << "popen() failed!" << endl;
        }

        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            placement += buffer;
        }

        pclose(pipe);
        cout << "Python returned: " << placement << endl;

    }

    auto nodes = nlohmann::json::parse(placement);

    // Create a vector of (nodeId, position) pairs and sort by nodeId numerically
    vector<pair<uint32_t, pair<double, double>>> sortedNodes;
    
    for (auto& [key, value] : nodes.items()) {
        uint32_t nodeId = stoi(key);
        double x = value[0];
        double y = value[1];
        sortedNodes.push_back({nodeId, {x, y}});
    }
    
    // Sort by node ID numerically (not alphabetically)
    sort(sortedNodes.begin(), sortedNodes.end(), 
              [](const auto& a, const auto& b) { return a.first < b.first; });

    // Now add positions in the correct order
    Ptr<ListPositionAllocator> listAlloc = CreateObject<ListPositionAllocator>();
    
    for (const auto& [nodeId, pos] : sortedNodes) {
        double x = pos.first;
        double y = pos.second;
        listAlloc->Add(Vector(x, y, 0.0));
        cout << "Loaded node " << nodeId << ": (" << x << ", " << y << ")\n";
    }

    MobilityHelper CameraMobility;
    CameraMobility.SetPositionAllocator(listAlloc);
    CameraMobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    CameraMobility.Install(camera_nodes);
}

NodeContainer CreateMobileNodes(){

    uint8_t node = 5;

    NodeContainer mobile_nodes;
    mobile_nodes.Create(node);

    SetMobileNodeMobility(node, mobile_nodes);
    NS_LOG_INFO("Setting mobile node mobility");
    return mobile_nodes;
}

void SetMobileNodeMobility(uint8_t node, NodeContainer& mobile_nodes){
    NS_LOG_INFO("Setting mobile node mobility");
    MobilityHelper MobileMobility;

    Ptr<ListPositionAllocator> initial = CreateObject<ListPositionAllocator>();
    initial->Add(Vector(0, 0, 0));
    MobileMobility.SetPositionAllocator(initial);


    MobileMobility.SetMobilityModel("ns3::RandomWalk2dOutdoorMobilityModel",
    "Bounds", RectangleValue(Rectangle(-200, 200, -200, 200)),
    "Speed", StringValue("ns3::UniformRandomVariable[Min=1.0|Max=20.0]"),
    "Direction", StringValue("ns3::UniformRandomVariable[Min=0|Max=6.2830]")
    );

    MobileMobility.Install(mobile_nodes);

}

NetDeviceContainer WifiStack(NodeContainer& nodes){
    // --- WiFi PHY ---
    YansWifiChannelHelper channel = YansWifiChannelHelper::Default();
    YansWifiPhyHelper phy;
    phy.SetChannel(channel.Create());
    phy.Set("TxPowerStart", DoubleValue(20.0));
    phy.Set("TxPowerEnd", DoubleValue(20.0));

    // --- WiFi MAC: Ad-hoc L2 ---
    WifiHelper wifi;
    wifi.SetStandard(WIFI_STANDARD_80211ax);

    WifiMacHelper mac;
    mac.SetType("ns3::AdhocWifiMac");
    NS_LOG_INFO("Installing WiFi devices on nodes");
    NetDeviceContainer devices = wifi.Install(phy, mac, nodes);
    return devices;
} 

map<string, NetDeviceContainer>& pointToPointEther(NodeContainer& nodes){
    // Find the neighbor within communication range from placement.json
    const char* file = "scratch/placement.json";
    double commRange = 50.0;

    // Call Python script to calculate neighbors
    string command = "python3 scratch/calculate_neighbors.py " + to_string((int)commRange) + " " + file + " neighbors.json";
    int result = system(command.c_str());
    
    if (result == 0) {
        NS_LOG_UNCOND("Successfully calculated neighbors and saved to neighbors.json");
    } else {
        NS_LOG_UNCOND("Error running calculate_neighbors.py");
    }
    
    // Read neighbors from neighbors.json
    ifstream inFile("neighbors.json");
    nlohmann::json neighborsJson;
    inFile >> neighborsJson;
    inFile.close();
    
    // Set up point-to-point links based on neighbors
    for (auto& [nodeId, neighborList] : neighborsJson.items()) {
        uint32_t nodeIndex = stoi(nodeId);
        Ptr<Node> node = nodes.Get(nodeIndex);
        for (auto& neighborId : neighborList) {
            uint32_t neighborIndex = neighborId.get<uint32_t>();
            
            if (neighborIndex == nodeIndex) continue;  // Skip self-links
            
            // Create unique key with smaller ID first
            uint32_t minId = min(nodeIndex, neighborIndex);
            uint32_t maxId = max(nodeIndex, neighborIndex);
            string key = to_string(minId) + "-" + to_string(maxId);
            
            // Skip if link already exists
            if (p2pDevices.find(key) != p2pDevices.end()) {
                continue;
            }
            
            Ptr<Node> neighborNode = nodes.Get(neighborIndex);
            
            // Set up point-to-point link between node and neighborNode
            PointToPointHelper p2p;
            p2p.SetDeviceAttribute("DataRate", StringValue("5Mbps"));
            p2p.SetChannelAttribute("Delay", StringValue("2ms"));
            NetDeviceContainer devices = p2p.Install(node, neighborNode);
            
            // Store devices in global map
            p2pDevices[key] = devices;
            
            NS_LOG_UNCOND("Established point-to-point link between Node " << nodeIndex << " and Node " << neighborIndex);
        }
    }
    
    // Print map
    for (const auto& [key, devices] : p2pDevices) {
        NS_LOG_UNCOND("Link: " << key << " with " << devices.GetN() << " devices");
    }
    NS_LOG_INFO("Point-to-point links established based on neighbors");
    
    return p2pDevices;
}


void animationSetup(AnimationInterface& anim, NodeContainer camera, NodeContainer mobile){
    for (uint8_t i = 0; i < mobile.GetN(); i++){
        anim.UpdateNodeColor(mobile.Get(i), 0, 0, 255);
    }

    for (uint8_t i = 0; i < camera.GetN(); i++){
        anim.UpdateNodeColor(camera.Get(i), 255, 0, 0);
    }

}

int main(int argc, char *argv[])
{
    // Configuration parameters
    uint32_t nNodes = 50;
    double commRange = 50.0;
    double simTime = 20;
    
    // Message generator parameters
    float minDensity = 0.0f;
    float maxDensity = 100.0f;
    float deltaDensity = 15.0f;
    double minInterval = 1.0;
    double maxInterval = 3.0;

    double percentage = 0.1;

    // Command line parsing
    CommandLine cmd;
    cmd.AddValue("percentage", "Percentage of neighbors to send messages to", percentage);
    cmd.AddValue("nNodes", "Number of nodes", nNodes);
    cmd.AddValue("commRange", "Communication range (meters)", commRange);
    cmd.AddValue("simTime", "Simulation time (seconds)", simTime);
    cmd.AddValue("minDensity", "Minimum density value", minDensity);
    cmd.AddValue("maxDensity", "Maximum density value", maxDensity);
    cmd.AddValue("deltaDensity", "Maximum density change", deltaDensity);
    cmd.AddValue("minInterval", "Minimum message interval (seconds)", minInterval);
    cmd.AddValue("maxInterval", "Maximum message interval (seconds)", maxInterval);
    cmd.Parse(argc, argv);
    
    LogComponentDisable("WifiPhy", LOG_LEVEL_INFO);
    LogComponentEnable("WifiMac", LOG_LEVEL_INFO);
    // Must be first (before any packets are created)
    PacketMetadata::Enable();

    LogComponentEnable("First", LOG_LEVEL_INFO);
    NS_LOG_UNCOND("Starting program");

    // Build network
    NodeContainer camera_nodes = CreateCameraNodes(nNodes);
    NodeContainer mobile_nodes = CreateMobileNodes();
    NetDeviceContainer devices = WifiStack(camera_nodes);

    pointToPointEther(camera_nodes);

    InternetStackHelper internet;
    internet.Install(camera_nodes);
    internet.Install(mobile_nodes);

    Ipv4AddressHelper address;
    map<uint32_t, Ipv4Address> nodeToIP;
    int subnetIndex = 0;
    for (auto& [key, devices] : p2pDevices) {
        ostringstream subnet;
        subnet << "10." << (subnetIndex / 256) << "." << (subnetIndex % 256) << ".0";
        address.SetBase(subnet.str().c_str(), "255.255.255.0");
        Ipv4InterfaceContainer interfaces = address.Assign(devices);
        subnetIndex++;
        
        // Store IPs
        try {
            size_t dash = key.find('-');
            uint32_t n1 = stoul(key.substr(0, dash));
            uint32_t n2 = stoul(key.substr(dash + 1));
            Ipv4Address ip1 = interfaces.GetAddress(0);
            Ipv4Address ip2 = interfaces.GetAddress(1);
            linkIPs[{n1, n2}] = {ip1, ip2};
            linkIPs[{n2, n1}] = {ip2, ip1};
        } catch (const exception& e) {
            NS_LOG_ERROR("Failed to parse key " << key << " or assign IPs: " << e.what());
        }
    }
    
    // Build IP mapping for gossip app
    for (uint32_t i = 0; i < camera_nodes.GetN(); ++i) {
        Ptr<Node> node = camera_nodes.Get(i);
        Ptr<Ipv4> ipv4 = node->GetObject<Ipv4>();
        // Use the first non-loopback interface address
        for (uint32_t j = 1; j < ipv4->GetNInterfaces(); ++j) {
            Ipv4Address addr = ipv4->GetAddress(j, 0).GetLocal();
            nodeToIP[i] = addr;
            break;
        }
    }

    // Get current date and time
    auto now = chrono::system_clock::now();
    auto time_t = chrono::system_clock::to_time_t(now);
    tm tm = *localtime(&time_t);
    stringstream ss;
    ss << put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    string timestamp = ss.str();

    // Create directory structure: test/YYYY-MM-DD_HH-MM-SS/
    string logDir = "test/" + timestamp + "/";
    filesystem::create_directories(logDir);

    // Install gossip application with message generator
    GossipAppHelper gossipApp(commRange, 9 , logDir, percentage, camera_nodes.GetN());
    gossipApp.SetMessageGenerator(minDensity, maxDensity, deltaDensity,
                                  Seconds(minInterval), Seconds(maxInterval));
    gossipApp.SetNodeToIP(nodeToIP);
    
    ApplicationContainer apps = gossipApp.Install(camera_nodes);
    apps.Start(Seconds(1.0));
    apps.Stop(Seconds(simTime - 1.0));


    AnimationInterface anim("test-netanim.xml");
    anim.SetMaxPktsPerTraceFile(5000000);
    anim.EnablePacketMetadata();
    animationSetup(anim, camera_nodes, mobile_nodes);

    // Optionally stop simulation after e.g. 2 seconds
    Simulator::Stop(Seconds(simTime));

    Simulator::Run();
    Simulator::Destroy();
    return 0;
}
