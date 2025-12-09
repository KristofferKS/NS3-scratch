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
        
        std::memcpy(buffer + offset, &nodeId, sizeof(nodeId));
        offset += sizeof(nodeId);

        std::memcpy(buffer + offset, &lastNode, sizeof(lastNode));
        offset += sizeof(lastNode);
        
        std::memcpy(buffer + offset, &reading.density, sizeof(reading.density));
        offset += sizeof(reading.density);
        
        std::memcpy(buffer + offset, &reading.timestamp, sizeof(reading.timestamp));
    }
    
    // Deserialize from buffer
    static DensityMessage Deserialize(const uint8_t* buffer)
    {
        DensityMessage msg;
        uint32_t offset = 0;
        
        std::memcpy(&msg.nodeId, buffer + offset, sizeof(msg.nodeId));
        offset += sizeof(msg.nodeId);

        std::memcpy(&msg.lastNode, buffer + offset,sizeof(msg.lastNode));
        offset += sizeof(msg.lastNode);

        std::memcpy(&msg.reading.density, buffer + offset, sizeof(msg.reading.density));
        offset += sizeof(msg.reading.density);
        
        std::memcpy(&msg.reading.timestamp, buffer + offset, sizeof(msg.reading.timestamp));
        
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
        float minAllowed = std::max(m_minDensity, m_currentDensity - m_deltaDensity);
        float maxAllowed = std::min(m_maxDensity, m_currentDensity + m_deltaDensity);

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

    virtual void Print(std::ostream &os) const override
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
    void SetupLogDir(std::string logDir);
    //void SetNodeToIP(const std::map<uint32_t, Ipv4Address>& nodeToIP);
    
    // Message generator configuration
    void SetupMessageGenerator(float minDensity, float maxDensity, 
                               float deltaDensity, Time minInterval, Time maxInterval);

protected:
    virtual void StartApplication(void);
    virtual void StopApplication(void);

private:
    void SendMessage(Ptr<NetDevice> dev, Mac48Address dst, const DensityMessage& message);
    void DiscoverPeers(void);
    void OnMessageGenerated(const DensityMessage& message);
    Mac48Address GetMacFromIpv4(const Ipv4Address& ipv4);
    bool ReceiveFromDevice(Ptr<NetDevice> device, Ptr<const Packet> packet, uint16_t protocol, const Address& src);

    Ptr<Socket> m_socket;
    double m_commRange;
    uint16_t m_port;
    std::vector<uint32_t> m_peers;
    std::ofstream m_logFile;
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

    std::string m_logDir;

    // Density knowledge for gossip protocol
    std::map<uint32_t, std::deque<std::pair<float, uint32_t>>> m_densityKnowledge;  // nodeId -> deque of (density, timestamp) pairs, max 5
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

void GossipApp::SetupLogDir(std::string logDir)
{
    m_logDir = logDir;
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

Mac48Address GossipApp::GetMacFromIpv4(const Ipv4Address& ip)
{
    for (uint32_t i = 0; i < m_nCameraNodes; ++i)
    {
        Ptr<Node> node = NodeList::GetNode(i);
        Ptr<Ipv4> ipv4 = node->GetObject<Ipv4>();
        if (ipv4)
        {
            for (uint32_t j = 0; j < ipv4->GetNInterfaces(); ++j)
            {
                for (uint32_t k = 0; k < ipv4->GetNAddresses(j); ++k)
                {
                    if (ipv4->GetAddress(j, k).GetLocal() == ip)
                    {
                        // Found the node and interface
                        Ptr<NetDevice> dev = node->GetDevice(j);
                        Mac48Address mac = Mac48Address::ConvertFrom(dev->GetAddress());
                        // Avoid returning all-zero MAC
                        if (mac != Mac48Address("00:00:00:00:00:00")) {
                            return mac;
                        }
                    }
                }
            }
        }
    }
    // Not found, return broadcast MAC as fallback
    return Mac48Address("ff:ff:ff:ff:ff:ff");
}

void GossipApp::DiscoverPeers(void)
{
    uint32_t myId = GetNode()->GetId();
    Ptr<MobilityModel> myMob = GetNode()->GetObject<MobilityModel>();
    Vector myPos = myMob->GetPosition();

    m_peers.clear();
    
    for (uint32_t i = 0; i < m_nCameraNodes; ++i)
    {
        if (i == myId) continue;

        Ptr<Node> peerNode = NodeList::GetNode(i);
        Ptr<MobilityModel> peerMob = peerNode->GetObject<MobilityModel>();
        Vector peerPos = peerMob->GetPosition();
        
        double dist = CalculateDistance(myPos, peerPos);
        
        if (dist <= m_commRange)
        {
            m_peers.push_back(i); // Store node index directly
        }
    }

    NS_LOG_INFO("Node " << myId << " at (" << myPos.x << "," << myPos.y 
                << ") found " << m_peers.size() << " peers within " 
                << m_commRange << "m range");
}

void
GossipApp::StartApplication(void)
{
    // Open log file for this node
    std::string fname = m_logDir + "gossip_log_node" + std::to_string(GetNode()->GetId()) + ".csv";
    m_logFile.open(fname, std::ios::out);
    m_logFile << "event,time,node_from,node_to,density,timestamp_ms" << std::endl;

    // Discover peers
    DiscoverPeers();

    Ptr<Node> node = GetNode();
    for (uint32_t i = 0; i < node->GetNDevices(); ++i)
    {
        Ptr<NetDevice> dev = node->GetDevice(i);
        if (DynamicCast<WifiNetDevice>(dev) != nullptr)  // Only WiFi devices
        {
            dev->SetReceiveCallback(MakeCallback(&GossipApp::ReceiveFromDevice, this));
        }
    }

    // Initialize density knowledge for all nodes
    for (uint32_t i = 0; i < m_nCameraNodes; ++i)
    {
        m_densityKnowledge[i] = std::deque<std::pair<float, uint32_t>>();  // Empty deque for each node
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

bool GossipApp::ReceiveFromDevice(Ptr<NetDevice> device, Ptr<const Packet> packet, uint16_t protocol, const Address& src)
{
    uint32_t myId = GetNode()->GetId();
    double now = Simulator::Now().GetSeconds();

    // Copy the packet and remove the header
    Ptr<Packet> copy = packet->Copy();
    GossipHeader header;
    copy->RemoveHeader(header);
    
    uint32_t senderId = header.GetNodeId();
    uint32_t lastNode = header.GetLastNode();

    // Ignore messages from self
    if (senderId == myId)
    {
        NS_LOG_DEBUG("Node " << myId << " ignoring self-message");
        return true;
    }

    // Log recv event to CSV
    m_logFile << "recv," << now << "," << senderId << "," << myId << ","
              << header.GetDensity() << "," << header.GetTimestamp() << std::endl;

    NS_LOG_INFO("Node " << myId << " received message from node " << senderId
                << " (density=" << header.GetDensity()
                << ", timestamp=" << header.GetTimestamp() << " ms) at " << now << "s");

    // Update density knowledge for the sender node
    if (m_densityKnowledge.find(senderId) != m_densityKnowledge.end())
    {
        // Check if the new timestamp is newer than existing ones
        uint32_t maxExistingTimestamp = 0;
        for (const auto& entry : m_densityKnowledge[senderId])
        {
            if (entry.second > maxExistingTimestamp)
            {
                maxExistingTimestamp = entry.second;
            }
        }

        if (header.GetTimestamp() > maxExistingTimestamp)
        {
            // Update knowledge
            m_densityKnowledge[senderId].push_back({header.GetDensity(), header.GetTimestamp()});
            // Keep only the last 5 entries
            if (m_densityKnowledge[senderId].size() > 5)
            {
                m_densityKnowledge[senderId].pop_front();
            }
            NS_LOG_DEBUG("Node " << myId << " updated knowledge for node " << senderId 
                         << ": density=" << header.GetDensity() << ", timestamp=" << header.GetTimestamp());

            // Log update event to CSV
            m_logFile << "update," << now << "," << senderId << "," << myId << ","
                      << header.GetDensity() << "," << header.GetTimestamp() << std::endl;

            NS_LOG_INFO("Node " << myId << " got a message from " << lastNode);
            // Prepare list of peers excluding sender and self
            std::vector<uint32_t> shuffledPeers;
            shuffledPeers.reserve(m_peers.size());
            for (auto peer : m_peers) {
                if (peer != lastNode && peer != myId) {
                    shuffledPeers.push_back(peer);
                }
            }

            if (shuffledPeers.empty()) {
                NS_LOG_DEBUG("Node " << myId << " has no peers to propagate to after excluding sender and self");
                return true;
            }

            // Propagate the update to a percentage of neighbors (excluding the sender, and self)
            double percentage = m_percentage;
            size_t nPeers = shuffledPeers.size();
            size_t nToSend = static_cast<size_t>(std::ceil(percentage * nPeers));
            if (nToSend == 0 && nPeers > 0) nToSend = 1;
            
            // Shuffle peers
            std::random_device rd;
            std::mt19937 g(rd());
            std::shuffle(shuffledPeers.begin(), shuffledPeers.end(), g);

            size_t limit = std::min(nToSend, shuffledPeers.size());
            size_t sentCount = 0;

            Ptr<Node> node = GetNode();
            for (uint32_t i = 0; i < node->GetNDevices(); ++i)
            {
                Ptr<NetDevice> dev = node->GetDevice(i);
                // Only use WiFiNetDevice
                if (DynamicCast<WifiNetDevice>(dev) == nullptr)
                    continue;

                for (size_t idx = 0; idx < limit && sentCount < limit; ++idx)
                {
                    uint32_t peerIndex = shuffledPeers[idx];

                    Ptr<Node> peerNode = NodeList::GetNode(peerIndex);
                    Ptr<NetDevice> peerDev;
                    for (uint32_t j = 0; j < peerNode->GetNDevices(); ++j) {
                        if (DynamicCast<WifiNetDevice>(peerNode->GetDevice(j)) != nullptr) {
                            peerDev = peerNode->GetDevice(j);
                            break;
                        }
                    }
                    if (!peerDev) continue; // No WiFi device found, skip

                    Mac48Address dstMac = Mac48Address::ConvertFrom(peerDev->GetAddress());
                    Mac48Address myMac = Mac48Address::ConvertFrom(dev->GetAddress());
                    if (dstMac == myMac) continue; // Skip self

                    // Create DensityMessage from header for sending
                    DensityMessage messageToSend(header.GetNodeId(), header.GetDensity(), header.GetTimestamp());
                    SendMessage(dev, dstMac, messageToSend);
                }
            }
        }
        else
        {
            NS_LOG_DEBUG("Node " << myId << " received outdated message from node " << senderId 
                        << " (timestamp=" << header.GetTimestamp() << " <= max=" << maxExistingTimestamp << "), dropping");

            // Log drop event to CSV
            m_logFile << "drop," << now << "," << senderId << "," << myId << ","
                    << header.GetDensity() << "," << header.GetTimestamp() << std::endl;
        }
    }

    return true; // Indicate packet was handled
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
    std::ofstream knowledgeFile(m_logDir + "final_knowledge.csv", std::ios::app);
    if (knowledgeFile.is_open())
    {
        // Add header if file is empty (first write)
        knowledgeFile.seekp(0, std::ios::end);
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
            knowledgeFile << std::endl;
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
        knowledgeFile << std::endl;
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
              << message.reading.density << "," << message.reading.timestamp << std::endl;

    double percentage = m_percentage;
    size_t nPeers = m_peers.size();
    size_t nToSend = static_cast<size_t>(std::ceil(percentage * nPeers));
    if (nToSend == 0 && nPeers > 0) nToSend = 1;

    // Shuffle peers
    std::vector<uint32_t> shuffledPeers = m_peers;
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(shuffledPeers.begin(), shuffledPeers.end(), g);

    Ptr<Node> node = GetNode();
    for (uint32_t i = 0; i < node->GetNDevices(); ++i)
    {
        Ptr<NetDevice> dev = node->GetDevice(i);
        // Only use WiFiNetDevice
        if (DynamicCast<WifiNetDevice>(dev) == nullptr)
            continue;

        for (size_t idx = 0; idx < nToSend; ++idx)
        {
            uint32_t peerIndex = shuffledPeers[idx];
            Ptr<Node> peerNode = NodeList::GetNode(peerIndex);
            Ptr<NetDevice> peerDev;
            for (uint32_t j = 0; j < peerNode->GetNDevices(); ++j) {
                if (DynamicCast<WifiNetDevice>(peerNode->GetDevice(j)) != nullptr) {
                    peerDev = peerNode->GetDevice(j);
                    break;
                }
            }
            if (!peerDev) continue; // No WiFi device found, skip

            Mac48Address dstMac = Mac48Address::ConvertFrom(peerDev->GetAddress());
            Mac48Address myMac = Mac48Address::ConvertFrom(dev->GetAddress());
            if (dstMac == myMac) continue;
            SendMessage(dev, dstMac, message);
        }
    }
}

void GossipApp::SendMessage(Ptr<NetDevice> dev, Mac48Address dst, const DensityMessage& message)
{
    if (m_peers.empty())
    {
        NS_LOG_WARN("Node " << GetNode()->GetId() << " has no peers to send to");
        return;
    }

    uint32_t lastNode = GetNode()->GetId();

    // Create packet and add custom header
    Ptr<Packet> packet = Create<Packet>();
    GossipHeader header(message.nodeId, lastNode, message.reading.density, message.reading.timestamp);
    packet->AddHeader(header);

    uint16_t ethertype = 0x88B5;

    Mac48Address macAddr = Mac48Address::ConvertFrom(dev->GetAddress());

    // Log send event to CSV
    double now = Simulator::Now().GetSeconds();
    m_logFile << "send," << now << "," << GetNode()->GetId() << "," << dst << ","
              << message.reading.density << "," << message.reading.timestamp << std::endl;

    dev->Send(packet, dst, ethertype);

    NS_LOG_INFO("Packet sent from " << macAddr << " to " << dst);
}


/**
 * GossipAppHelper
 */
class GossipAppHelper
{
public:
    GossipAppHelper(double range = 50.0, uint16_t port = 9, std::string logDir = "", double percentage = 0.5, uint32_t nCameraNodes = 50);
    
    // Set message generator parameters
    void SetMessageGenerator(float minDensity, 
                                float maxDensity, 
                                float deltaDensity, 
                                Time minInterval, 
                                Time maxInterval);
    
    ApplicationContainer Install(NodeContainer c);
    ApplicationContainer Install(Ptr<Node> node);

private:
    double m_percentage;
    double m_range;
    uint16_t m_port;
    std::string m_logDir;
    uint32_t m_nCameraNodes;
    
    // Message generator parameters
    float m_minDensity;
    float m_maxDensity;
    float m_deltaDensity;
    Time m_minInterval;
    Time m_maxInterval;
    bool m_useMessageGenerator;
};

GossipAppHelper::GossipAppHelper(double range, uint16_t port, std::string logDir, double percentage, uint32_t nCameraNodes)
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


NodeContainer CreateCameraNodes();
void SetCameraMobility(uint8_t node, NodeContainer& camera_nodes);
NodeContainer CreateMobileNodes();
void SetMobileNodeMobility(uint8_t node, NodeContainer& mobile_nodes);

NodeContainer CreateCameraNodes() {
    // Log the creation of camera nodes
    NS_LOG_INFO("Creating camera nodes");
    uint8_t node = 50;

    NodeContainer camera_nodes;
    camera_nodes.Create(node);

    NS_LOG_INFO("Setting camera node mobility");

    SetCameraMobility(node, camera_nodes);
    return camera_nodes;
}

void SetCameraMobility(uint8_t node, NodeContainer& camera_nodes) {
    NS_LOG_INFO("Setting camera node mobility");
    double radius = 200.0;         // radius of disc
    double wifiRange = 50.0;       // Wi-Fi range
    double centerX = 0.0;
    double centerY = 0.0;

    NS_LOG_INFO("Placing camera nodes with coverage constraints");
    
    // Path to the directory
    const char* node_file = "scratch/placement.json";

    // Structure which would store the metadata
    struct stat sb;

    std::string placement;

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
        std::string command = "python3 scratch/node_placement.py " + std::to_string(node) + " " + std::to_string(wifiRange) + " " + std::to_string(centerX) + " " + std::to_string(centerY) + " " + std::to_string(radius);
        // Read positions from file
        char buffer[128];
        FILE* pipe = popen(command.c_str(), "r");
        if (!pipe) {
            std::cerr << "popen() failed!" << std::endl;
        }

        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            placement += buffer;
        }

        pclose(pipe);
        std::cout << "Python returned: " << placement << std::endl;

    }

    auto nodes = nlohmann::json::parse(placement);

    // Create a vector of (nodeId, position) pairs and sort by nodeId numerically
    std::vector<std::pair<uint32_t, std::pair<double, double>>> sortedNodes;
    
    for (auto& [key, value] : nodes.items()) {
        uint32_t nodeId = std::stoi(key);
        double x = value[0];
        double y = value[1];
        sortedNodes.push_back({nodeId, {x, y}});
    }
    
    // Sort by node ID numerically (not alphabetically)
    std::sort(sortedNodes.begin(), sortedNodes.end(), 
              [](const auto& a, const auto& b) { return a.first < b.first; });

    // Now add positions in the correct order
    Ptr<ListPositionAllocator> listAlloc = CreateObject<ListPositionAllocator>();
    
    for (const auto& [nodeId, pos] : sortedNodes) {
        double x = pos.first;
        double y = pos.second;
        listAlloc->Add(Vector(x, y, 0.0));
        std::cout << "Loaded node " << nodeId << ": (" << x << ", " << y << ")\n";
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


void TxCallback(Ptr<const Packet> packet)
{
    NS_LOG_UNCOND("Packet transmitted at " << Simulator::Now().GetSeconds() << "s: " << packet->GetSize() << " bytes");
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
    uint32_t nNodes = 4;
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
    NodeContainer camera_nodes = CreateCameraNodes();
    NodeContainer mobile_nodes = CreateMobileNodes();
    NetDeviceContainer devices = WifiStack(camera_nodes);

    // Get current date and time
    auto now = std::chrono::system_clock::now();
    auto time_t = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&time_t);
    std::stringstream ss;
    ss << std::put_time(&tm, "%Y-%m-%d_%H-%M-%S");
    std::string timestamp = ss.str();

    // Create directory structure: test/YYYY-MM-DD_HH-MM-SS/
    std::string logDir = "test/" + timestamp + "/";
    std::filesystem::create_directories(logDir);

    // Install gossip application with message generator
    GossipAppHelper gossipApp(commRange, 9 , logDir, percentage, camera_nodes.GetN());
    gossipApp.SetMessageGenerator(minDensity, maxDensity, deltaDensity,
                                  Seconds(minInterval), Seconds(maxInterval));
    
    ApplicationContainer apps = gossipApp.Install(camera_nodes);
    apps.Start(Seconds(1.0));
    apps.Stop(Seconds(simTime - 1.0));


    //AnimationInterface anim("test-netanim.xml");
    //anim.SetMaxPktsPerTraceFile(5000000);
    //anim.EnablePacketMetadata();
    //animationSetup(anim, camera_nodes, mobile_nodes);

    // Optionally stop simulation after e.g. 2 seconds
    Simulator::Stop(Seconds(simTime));

    Simulator::Run();
    Simulator::Destroy();
    return 0;
}
