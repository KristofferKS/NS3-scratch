#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/applications-module.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("CustomApp");

// Custom application class
class MyApp : public Application
{
public:
    MyApp();
    virtual ~MyApp();
    
    void Setup(Ptr<Socket> socket, Address address, uint32_t packetSize,
               uint32_t nPackets, DataRate dataRate);

private:
    virtual void StartApplication(void);
    virtual void StopApplication(void);
    
    void ScheduleTx(void);
    void SendPacket(void);
    
    Ptr<Socket>     m_socket;
    Address         m_peer;
    uint32_t        m_packetSize;
    uint32_t        m_nPackets;
    DataRate        m_dataRate;
    EventId         m_sendEvent;
    bool            m_running;
    uint32_t        m_packetsSent;
};

MyApp::MyApp()
    : m_socket(0),
      m_peer(),
      m_packetSize(0),
      m_nPackets(0),
      m_dataRate(0),
      m_sendEvent(),
      m_running(false),
      m_packetsSent(0)
{
}

MyApp::~MyApp()
{
    m_socket = 0;
}

void MyApp::Setup(Ptr<Socket> socket, Address address, uint32_t packetSize,
                  uint32_t nPackets, DataRate dataRate)
{
    m_socket = socket;
    m_peer = address;
    m_packetSize = packetSize;
    m_nPackets = nPackets;
    m_dataRate = dataRate;
}

void MyApp::StartApplication(void)
{
    m_running = true;
    m_packetsSent = 0;
    m_socket->Bind();
    m_socket->Connect(m_peer);
    SendPacket();
}

void MyApp::StopApplication(void)
{
    m_running = false;
    
    if (m_sendEvent.IsPending())
    {
        Simulator::Cancel(m_sendEvent);
    }
    
    if (m_socket)
    {
        m_socket->Close();
    }
}

void MyApp::SendPacket(void)
{
    Ptr<Packet> packet = Create<Packet>(m_packetSize);
    m_socket->Send(packet);
    
    NS_LOG_INFO("Sent packet " << m_packetsSent << " of size " << m_packetSize);
    
    if (++m_packetsSent < m_nPackets)
    {
        ScheduleTx();
    }
}

void MyApp::ScheduleTx(void)
{
    if (m_running)
    {
        Time tNext(Seconds(m_packetSize * 8 / static_cast<double>(m_dataRate.GetBitRate())));
        m_sendEvent = Simulator::Schedule(tNext, &MyApp::SendPacket, this);
    }
}

// Receive callback
static void RxCallback(Ptr<Socket> socket)
{
    Ptr<Packet> packet;
    Address from;
    while ((packet = socket->RecvFrom(from)))
    {
        if (packet->GetSize() > 0)
        {
            InetSocketAddress iaddr = InetSocketAddress::ConvertFrom(from);
            NS_LOG_INFO("Received " << packet->GetSize() << " bytes from " 
                        << iaddr.GetIpv4());
        }
    }
}

int main(int argc, char *argv[])
{
    LogComponentEnable("CustomApp", LOG_LEVEL_INFO);
    
    // Create nodes
    NodeContainer nodes;
    nodes.Create(2);
    
    // WiFi setup
    YansWifiChannelHelper channel = YansWifiChannelHelper::Default();
    YansWifiPhyHelper phy;
    phy.SetChannel(channel.Create());
    
    WifiHelper wifi;
    wifi.SetStandard(WIFI_STANDARD_80211ax);
    wifi.SetRemoteStationManager("ns3::IdealWifiManager");
    
    WifiMacHelper mac;
    Ssid ssid = Ssid("test-network");
    mac.SetType("ns3::AdhocWifiMac", "Ssid", SsidValue(ssid));
    
    NetDeviceContainer devices = wifi.Install(phy, mac, nodes);
    
    // Mobility
    MobilityHelper mobility;
    mobility.SetPositionAllocator("ns3::GridPositionAllocator",
                                  "MinX", DoubleValue(0.0),
                                  "MinY", DoubleValue(0.0),
                                  "DeltaX", DoubleValue(10.0),
                                  "DeltaY", DoubleValue(0.0),
                                  "GridWidth", UintegerValue(2),
                                  "LayoutType", StringValue("RowFirst"));
    mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    mobility.Install(nodes);
    
    // Internet stack
    InternetStackHelper stack;
    stack.Install(nodes);
    
    Ipv4AddressHelper address;
    address.SetBase("10.1.1.0", "255.255.255.0");
    Ipv4InterfaceContainer interfaces = address.Assign(devices);
    
    // Create sockets
    TypeId tid = TypeId::LookupByName("ns3::UdpSocketFactory");
    
    // Receiver socket (node 1)
    Ptr<Socket> recvSink = Socket::CreateSocket(nodes.Get(1), tid);
    InetSocketAddress local = InetSocketAddress(Ipv4Address::GetAny(), 9);
    recvSink->Bind(local);
    recvSink->SetRecvCallback(MakeCallback(&RxCallback));
    
    // Sender socket (node 0)
    Ptr<Socket> source = Socket::CreateSocket(nodes.Get(0), tid);
    InetSocketAddress remote = InetSocketAddress(interfaces.GetAddress(1), 9);
    
    // Create and install custom app
    Ptr<MyApp> app = CreateObject<MyApp>();
    app->Setup(source, remote, 512, 20, DataRate("1Mbps"));
    nodes.Get(0)->AddApplication(app);
    app->SetStartTime(Seconds(1.0));
    app->SetStopTime(Seconds(10.0));
    
    Simulator::Stop(Seconds(11.0));
    Simulator::Run();
    Simulator::Destroy();
    return 0;
}