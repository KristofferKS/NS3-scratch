#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/applications-module.h"
#include "ns3/netanim-module.h"
#include "ns3/buildings-module.h"
#include <cstdint>  // for uint32_t
#include <iostream>
#include <sys/stat.h>
#include <fstream>
#include <nlohmann/json.hpp>


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
    void SendMessage(Ptr<NetDevice> dev, Mac48Address dst, const DensityMessage& message);
    void ReceivePacket(Ptr<Socket> socket);
    void DiscoverPeers(void);
    void OnMessageGenerated(const DensityMessage& message);
    Mac48Address GetMacFromIpv4(const Ipv4Address& ipv4);
    bool ReceiveFromDevice(Ptr<NetDevice> device, Ptr<const Packet> packet, uint16_t protocol, const Address& src);

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

Mac48Address GossipApp::GetMacFromIpv4(const Ipv4Address& ip)
{
    for (uint32_t i = 0; i < NodeList::GetNNodes(); ++i)
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
                        return Mac48Address::ConvertFrom(dev->GetAddress());
                    }
                }
            }
        }
    }
    // Not found, return broadcast MAC as fallback
    return Mac48Address("ff:ff:ff:ff:ff:ff");
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
    /*
    // Create and bind UDP socket
    if (!m_socket)
    {
        m_socket = Socket::CreateSocket(GetNode(), UdpSocketFactory::GetTypeId());
        m_socket->Bind(InetSocketAddress(Ipv4Address::GetAny(), m_port));
        m_socket->SetRecvCallback(MakeCallback(&GossipApp::ReceivePacket, this));
    }*/

    // Open log file for this node
    std::string fname = "gossip_log_node" + std::to_string(GetNode()->GetId()) + ".csv";
    m_logFile.open(fname, std::ios::out);
    m_logFile << "event,time,node_from,node_to,density,timestamp_ms" << std::endl;

    // Discover peers
    DiscoverPeers();

    Ptr<Node> node = GetNode();
    for (uint32_t i = 0; i < node->GetNDevices(); ++i)
    {
        Ptr<NetDevice> dev = node->GetDevice(i);
        dev->SetReceiveCallback(MakeCallback(&GossipApp::ReceiveFromDevice, this));
    }

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

bool GossipApp::ReceiveFromDevice(Ptr<NetDevice> device, Ptr<const Packet> packet, uint16_t protocol, const Address& src)
{
    uint32_t myId = GetNode()->GetId();
    double now = Simulator::Now().GetSeconds();

    // Deserialize message
    uint32_t size = packet->GetSize();
    uint8_t* buffer = new uint8_t[size];
    packet->CopyData(buffer, size);

    DensityMessage message = DensityMessage::Deserialize(buffer);
    delete[] buffer;
/*
    NS_LOG_INFO("Node " << myId << " received MAC packet from "
                << Mac48Address::ConvertFrom(src) << " (source node=" << message.nodeId
                << ", density=" << message.reading.density
                << ", timestamp=" << message.reading.timestamp << "ms)");*/

    m_logFile << "recv," << now << "," << message.nodeId << "," << myId << ","
              << message.reading.density << "," << message.reading.timestamp << std::endl;

    return true; // Indicate packet was handled
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

    // Replace SendMessage with SendPacket for all peers
    Ptr<Node> node = GetNode();
    for (uint32_t i = 0; i < node->GetNDevices(); ++i)
    {
        Ptr<NetDevice> dev = node->GetDevice(i);
        //Mac48Address srcMac = Mac48Address::ConvertFrom(dev->GetAddress());
        for (const auto& peerAddr : m_peers)
        {
            Mac48Address dstMac = GetMacFromIpv4(peerAddr);
            SendMessage(dev, dstMac, message);
        }
    }
}

void
GossipApp::SendMessage(Ptr<NetDevice> dev, Mac48Address dst, const DensityMessage& message)
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

    uint16_t ethertype = 0x88B5;

    WifiMacHeader hdr;
    hdr.SetType(WIFI_MAC_DATA);
    Mac48Address macAddr = Mac48Address::ConvertFrom(dev->GetAddress());
    hdr.SetAddr1(dst);
    hdr.SetAddr2(macAddr);
    hdr.SetAddr3(macAddr);
    hdr.SetDsNotFrom();
    hdr.SetDsNotTo();

    dev->Send(packet, dst, ethertype);

    NS_LOG_INFO("Packet sent from " << macAddr << " to " << dst);
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

struct NodeDate{
    int nodeid;
    float density;
    int64_t timestamp;
};

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
        std::string command = "python3 /node_placement.py " + std::to_string(node) + " " + std::to_string(wifiRange) + " " + std::to_string(centerX) + " " + std::to_string(centerY) + " " + std::to_string(radius);
        //system(command.c_str());
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

    // Set up camera node mobility
    Ptr<ListPositionAllocator> listAlloc = CreateObject<ListPositionAllocator>();
        // Parse JSON → add positions
    for (auto& [key, value] : nodes.items()) {
        double x = value[0];
        double y = value[1];
        listAlloc->Add(Vector(x, y, 0.0));
        std::cout << "Loaded node " << key << ": (" << x << ", " << y << ")\n";
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

void LoraStack(){
    ;
}

void Cellular(){
    ;
}

void CostumHeader(){
    ;
}

void GossipProtocol(){
    ;
}

void neighborDiscovery(){
    ;
}

void TxCallback(Ptr<const Packet> packet)
{
    NS_LOG_UNCOND("Packet transmitted at " << Simulator::Now().GetSeconds() << "s: " << packet->GetSize() << " bytes");
}

/*
void SendPacket(Ptr<NetDevice> dev, Mac48Address dst)
{
    NS_LOG_UNCOND("SendPacket() called at " << Simulator::Now().GetSeconds() << "s");

    Ptr<Packet> packet = Create<Packet>(100); // 100-byte payload
    uint16_t ethertype = 0x88B5;

    WifiMacHeader hdr;
    hdr.SetType(WIFI_MAC_DATA);
    Mac48Address macAddr = Mac48Address::ConvertFrom(dev->GetAddress());
    hdr.SetAddr1(dst);
    hdr.SetAddr2(macAddr);
    hdr.SetAddr3(macAddr);
    hdr.SetDsNotFrom();
    hdr.SetDsNotTo();

    dev->Send(packet, dst, ethertype);
}
*/

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
    double simTime = 5;
    
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

    //Ptr<NetDevice> dev1 = devices.Get(1);
    //Mac48Address addr1 = Mac48Address::ConvertFrom(dev1->GetAddress());
    //NS_LOG_INFO("Device 1 MAC Address: " << addr1);

    // Prepare device
    //Ptr<NetDevice> dev0 = devices.Get(0);
    InternetStackHelper internet;
    internet.Install(camera_nodes);

    Ipv4AddressHelper ipv4;
    ipv4.SetBase("10.1.1.0", "255.255.255.0");
    ipv4.Assign(devices);

    // Install gossip application with message generator
    GossipAppHelper gossipApp(commRange);
    gossipApp.SetMessageGenerator(minDensity, maxDensity, deltaDensity,
                                  Seconds(minInterval), Seconds(maxInterval));
    
    ApplicationContainer apps = gossipApp.Install(camera_nodes);
    apps.Start(Seconds(1.0));
    apps.Stop(Seconds(simTime - 1.0));


    AnimationInterface anim("test-netanim.xml");
    anim.EnablePacketMetadata();
    animationSetup(anim, camera_nodes, mobile_nodes);

    // Schedule packet send at t = 0.1s
    //Simulator::Schedule(Seconds(1), &SendPacket, dev0, addr1);

    // Optionally stop simulation after e.g. 2 seconds
    Simulator::Stop(Seconds(simTime));

    Simulator::Run();
    Simulator::Destroy();
    return 0;
}
