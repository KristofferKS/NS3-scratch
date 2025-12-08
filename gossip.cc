/**
 * Gossip Protocol Simulation with MessageGenerator
 * Single-file version for NS-3 scratch directory
 *
 * Usage: ./ns3 run scratch/gossip-with-messages
 */

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/netanim-module.h"
#include <vector>
#include <random>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("GossipWithMessages");

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
    DensityReading reading;

    DensityMessage() : nodeId(0) {}
    DensityMessage(uint32_t id, float density, uint32_t timestamp)
        : nodeId(id), reading(density, timestamp) {}
    
    // Get serialized size
    uint32_t GetSerializedSize() const
    {
        return sizeof(nodeId) + sizeof(reading.density) + sizeof(reading.timestamp);
    }
    
    // Serialize to buffer
    void Serialize(uint8_t* buffer) const
    {
        uint32_t offset = 0;
        
        std::memcpy(buffer + offset, &nodeId, sizeof(nodeId));
        offset += sizeof(nodeId);
        
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
        if (!m_running)
        {
            return;
        }

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
        if (!m_running)
        {
            return;
        }

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

/**
 * GossipApp - Gossip protocol application with message generation
 */
class GossipApp : public Application
{
public:
    static TypeId GetTypeId(void);
    GossipApp();
    virtual ~GossipApp();

    void Setup(double range, uint16_t port);
    
    // Message generator configuration
    void SetupMessageGenerator(float minDensity, float maxDensity, 
                               float deltaDensity, Time minInterval, Time maxInterval);

protected:
    virtual void StartApplication(void);
    virtual void StopApplication(void);

private:
    void SendMessage(const DensityMessage& message);
    void ReceivePacket(Ptr<Socket> socket);
    void DiscoverPeers(void);
    void OnMessageGenerated(const DensityMessage& message);

    Ptr<Socket> m_socket;
    double m_commRange;
    uint16_t m_port;
    std::vector<Ipv4Address> m_peers;
    std::ofstream m_logFile;
    
    // Message generator
    Ptr<MessageGenerator> m_messageGenerator;
    
    // Message generator parameters (stored until StartApplication)
    bool m_messageGeneratorConfigured;
    float m_minDensity;
    float m_maxDensity;
    float m_deltaDensity;
    Time m_minInterval;
    Time m_maxInterval;
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

void
GossipApp::Setup(double range, uint16_t port)
{
    m_commRange = range;
    m_port = port;
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

void
GossipApp::DiscoverPeers(void)
{
    uint32_t myId = GetNode()->GetId();
    Ptr<MobilityModel> myMob = GetNode()->GetObject<MobilityModel>();
    Vector myPos = myMob->GetPosition();

    m_peers.clear();
    
    for (uint32_t i = 0; i < NodeList::GetNNodes(); ++i)
    {
        if (i == myId) continue;

        Ptr<Node> peerNode = NodeList::GetNode(i);
        Ptr<MobilityModel> peerMob = peerNode->GetObject<MobilityModel>();
        Vector peerPos = peerMob->GetPosition();
        
        double dist = CalculateDistance(myPos, peerPos);
        
        if (dist <= m_commRange)
        {
            Ptr<Ipv4> ipv4 = peerNode->GetObject<Ipv4>();
            if (ipv4 && ipv4->GetNInterfaces() > 1)
            {
                Ipv4Address peerAddr = ipv4->GetAddress(1, 0).GetLocal();
                m_peers.push_back(peerAddr);
            }
        }
    }

    NS_LOG_INFO("Node " << myId << " at (" << myPos.x << "," << myPos.y 
                << ") found " << m_peers.size() << " peers within " 
                << m_commRange << "m range");
}

void
GossipApp::StartApplication(void)
{
    // Create and bind UDP socket
    if (!m_socket)
    {
        m_socket = Socket::CreateSocket(GetNode(), UdpSocketFactory::GetTypeId());
        m_socket->Bind(InetSocketAddress(Ipv4Address::GetAny(), m_port));
        m_socket->SetRecvCallback(MakeCallback(&GossipApp::ReceivePacket, this));
    }

    // Open log file for this node
    std::string fname = "gossip_log_node" + std::to_string(GetNode()->GetId()) + ".csv";
    m_logFile.open(fname, std::ios::out);
    m_logFile << "event,time,node_from,node_to,density,timestamp_ms" << std::endl;

    // Discover peers
    DiscoverPeers();

    // Setup message generator now that node is available
    if (m_messageGeneratorConfigured)
    {
        uint32_t nodeId = GetNode()->GetId();
        m_messageGenerator->Setup(nodeId, m_minDensity, m_maxDensity, m_deltaDensity,
                                 m_minInterval, m_maxInterval);
        
        // Set callback for when messages are generated
        m_messageGenerator->SetMessageGeneratedCallback(
            MakeCallback(&GossipApp::OnMessageGenerated, this));
    }

    // Start message generator
    m_messageGenerator->Start();
}

void
GossipApp::StopApplication(void)
{
    // Stop message generator
    m_messageGenerator->Stop();

    if (m_socket)
    {
        m_socket->Close();
    }
    if (m_logFile.is_open())
    {
        m_logFile.close();
    }
}

void
GossipApp::OnMessageGenerated(const DensityMessage& message)
{
    NS_LOG_INFO("Node " << GetNode()->GetId() 
                << " received new message from generator: density=" 
                << message.reading.density << ", timestamp=" 
                << message.reading.timestamp << "ms");
    
    // Send the generated message to all peers
    SendMessage(message);
}

void
GossipApp::SendMessage(const DensityMessage& message)
{
    if (m_peers.empty())
    {
        NS_LOG_WARN("Node " << GetNode()->GetId() << " has no peers to send to");
        return;
    }

    // Serialize message
    uint32_t msgSize = message.GetSerializedSize();
    uint8_t* buffer = new uint8_t[msgSize];
    message.Serialize(buffer);
    
    Ptr<Packet> packet = Create<Packet>(buffer, msgSize);
    delete[] buffer;

    uint32_t myId = GetNode()->GetId();
    double now = Simulator::Now().GetSeconds();

    // Send to all peers within range
    for (const auto& peerAddr : m_peers)
    {
        m_socket->SendTo(packet->Copy(), 0, InetSocketAddress(peerAddr, m_port));
        NS_LOG_INFO("Node " << myId << " sent density message to " << peerAddr 
                    << " (density=" << message.reading.density << ")");
        
        m_logFile << "send," << now << "," << myId << "," << peerAddr << "," 
                  << message.reading.density << "," << message.reading.timestamp << std::endl;
    }
}

void
GossipApp::ReceivePacket(Ptr<Socket> socket)
{
    Ptr<Packet> packet;
    Address from;

    while ((packet = socket->RecvFrom(from)))
    {
        uint32_t myId = GetNode()->GetId();
        InetSocketAddress fromAddr = InetSocketAddress::ConvertFrom(from);
        double now = Simulator::Now().GetSeconds();
        
        // Deserialize message
        uint32_t size = packet->GetSize();
        uint8_t* buffer = new uint8_t[size];
        packet->CopyData(buffer, size);
        
        DensityMessage message = DensityMessage::Deserialize(buffer);
        delete[] buffer;
        
        NS_LOG_INFO("Node " << myId << " received density message from " 
                    << fromAddr.GetIpv4() << " (source node=" << message.nodeId 
                    << ", density=" << message.reading.density 
                    << ", timestamp=" << message.reading.timestamp << "ms)");
        
        m_logFile << "recv," << now << "," << message.nodeId << "," << myId << "," 
                  << message.reading.density << "," << message.reading.timestamp << std::endl;
    }
}

/**
 * GossipAppHelper
 */
class GossipAppHelper
{
public:
    GossipAppHelper(double range = 50.0, uint16_t port = 9);
    
    // Set message generator parameters
    void SetMessageGenerator(float minDensity, float maxDensity, 
                            float deltaDensity, Time minInterval, Time maxInterval);
    
    ApplicationContainer Install(NodeContainer c);
    ApplicationContainer Install(Ptr<Node> node);

private:
    double m_range;
    uint16_t m_port;
    
    // Message generator parameters
    float m_minDensity;
    float m_maxDensity;
    float m_deltaDensity;
    Time m_minInterval;
    Time m_maxInterval;
    bool m_useMessageGenerator;
};

GossipAppHelper::GossipAppHelper(double range, uint16_t port)
    : m_range(range), 
      m_port(port),
      m_minDensity(0.0f),
      m_maxDensity(100.0f),
      m_deltaDensity(10.0f),
      m_minInterval(Seconds(1.0)),
      m_maxInterval(Seconds(5.0)),
      m_useMessageGenerator(true)
{
}

void
GossipAppHelper::SetMessageGenerator(float minDensity, float maxDensity, 
                                     float deltaDensity, Time minInterval, Time maxInterval)
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
    app->Setup(m_range, m_port);
    
    if (m_useMessageGenerator)
    {
        app->SetupMessageGenerator(m_minDensity, m_maxDensity, m_deltaDensity,
                                  m_minInterval, m_maxInterval);
    }
    
    node->AddApplication(app);
    return ApplicationContainer(app);
}

/**
 * Main simulation function
 */
int main(int argc, char *argv[])
{
    // Configuration parameters
    uint32_t nNodes = 4;
    double commRange = 120.0;
    double simTime = 20.0;
    
    // Message generator parameters
    float minDensity = 0.0f;
    float maxDensity = 100.0f;
    float deltaDensity = 15.0f;
    double minInterval = 1.0;
    double maxInterval = 3.0;

    // Command line parsing
    CommandLine cmd;
    cmd.AddValue("nNodes", "Number of nodes", nNodes);
    cmd.AddValue("commRange", "Communication range (meters)", commRange);
    cmd.AddValue("simTime", "Simulation time (seconds)", simTime);
    cmd.AddValue("minDensity", "Minimum density value", minDensity);
    cmd.AddValue("maxDensity", "Maximum density value", maxDensity);
    cmd.AddValue("deltaDensity", "Maximum density change", deltaDensity);
    cmd.AddValue("minInterval", "Minimum message interval (seconds)", minInterval);
    cmd.AddValue("maxInterval", "Maximum message interval (seconds)", maxInterval);
    cmd.Parse(argc, argv);

    // Enable logging
    LogComponentEnable("GossipWithMessages", LOG_LEVEL_INFO);
    
    NS_LOG_INFO("Simulation Configuration:");
    NS_LOG_INFO("  Nodes: " << nNodes);
    NS_LOG_INFO("  Communication Range: " << commRange << "m");
    NS_LOG_INFO("  Simulation Time: " << simTime << "s");
    NS_LOG_INFO("  Density Range: [" << minDensity << ", " << maxDensity << "]");
    NS_LOG_INFO("  Max Density Change: " << deltaDensity);
    NS_LOG_INFO("  Message Interval: [" << minInterval << "s, " << maxInterval << "s]");

    // Create nodes
    NodeContainer nodes;
    nodes.Create(nNodes);

    // Setup WiFi
    WifiHelper wifi;
    wifi.SetStandard(WIFI_STANDARD_80211be);

    YansWifiChannelHelper channel;
    channel.SetPropagationDelay("ns3::ConstantSpeedPropagationDelayModel");
    channel.AddPropagationLoss("ns3::RangePropagationLossModel",
                               "MaxRange", DoubleValue(commRange + 0.5));
    
    YansWifiPhyHelper phy;
    phy.SetChannel(channel.Create());
    
    double txPower = 16.0206;
    if (commRange < 50)
        txPower = 10.0;
    else if (commRange < 150)
        txPower = 20.0;
    else
        txPower = 30.0;
    
    phy.Set("TxPowerStart", DoubleValue(txPower));
    phy.Set("TxPowerEnd", DoubleValue(txPower));

    WifiMacHelper mac;
    mac.SetType("ns3::AdhocWifiMac");

    NetDeviceContainer devices = wifi.Install(phy, mac, nodes);

    // Install Internet stack
    InternetStackHelper stack;
    stack.Install(nodes);

    // Assign IP addresses
    Ipv4AddressHelper address;
    address.SetBase("10.0.0.0", "255.255.255.0");
    address.Assign(devices);

    // Setup mobility - grid layout
    MobilityHelper mobility;
    Ptr<ListPositionAllocator> positionAlloc = CreateObject<ListPositionAllocator>();
    positionAlloc->Add(Vector(0, 0, 0));
    positionAlloc->Add(Vector(0, 100, 0));
    positionAlloc->Add(Vector(100, 0, 0));
    positionAlloc->Add(Vector(100, 100, 0));
    positionAlloc->Add(Vector(50, 50, 0));

    mobility.SetPositionAllocator(positionAlloc);
    mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    mobility.Install(nodes);

    // Install gossip application with message generator
    GossipAppHelper gossipApp(commRange);
    gossipApp.SetMessageGenerator(minDensity, maxDensity, deltaDensity,
                                  Seconds(minInterval), Seconds(maxInterval));
    
    ApplicationContainer apps = gossipApp.Install(nodes);
    apps.Start(Seconds(1.0));
    apps.Stop(Seconds(simTime - 1.0));

    // Setup animation
    AnimationInterface anim("gossip-animation.xml");
    anim.SetMaxPktsPerTraceFile(500000);

    for (uint32_t i = 0; i < nodes.GetN(); i++)
    {
        anim.UpdateNodeDescription(nodes.Get(i), "Node " + std::to_string(i));
        anim.UpdateNodeColor(nodes.Get(i), 255, 0, 0);
    }

    // Run simulation
    Simulator::Stop(Seconds(simTime));
    Simulator::Run();
    Simulator::Destroy();

    return 0;
}