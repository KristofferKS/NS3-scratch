#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/applications-module.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("MeshNetwork");

void GenerateTraffic(Ptr<Socket> socket, uint32_t pktSize, uint32_t pktCount, Time pktInterval)
{
    if (pktCount > 0)
    {
        socket->Send(Create<Packet>(pktSize));
        Simulator::Schedule(pktInterval, &GenerateTraffic, socket, pktSize, pktCount - 1, pktInterval);
    }
    else
    {
        socket->Close();
    }
}

int main(int argc, char *argv[])
{
    uint32_t nNodes = 5;

    LogComponentEnable("MeshNetwork", LOG_LEVEL_INFO);

    // Create nodes
    NodeContainer nodes;
    nodes.Create(nNodes);

    // Wifi
    WifiHelper wifi;
    wifi.SetStandard(WIFI_STANDARD_80211ax);
    wifi.SetRemoteStationManager("ns3::IdealWifiManager");

    WifiMacHelper mac;
    mac.SetType("ns3::AdhocWifiMac");

    // WiFi Adhoc
    YansWifiChannelHelper channel = YansWifiChannelHelper::Default();
    YansWifiPhyHelper phy;
    phy.SetChannel(channel.Create());

    NetDeviceContainer devices = wifi.Install(phy, mac, nodes);

    // Position nodes in  Circle
    MobilityHelper mobility;
    mobility.SetPositionAllocator("ns3::UniformDiscPositionAllocator",
                                  "X", DoubleValue(0.0),
                                  "Y", DoubleValue(0.0),
                                  "rho", DoubleValue(20.0));
    mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    mobility.Install(nodes);
    
    InternetStackHelper stack;
    stack.Install(nodes);

    Ipv4AddressHelper address;
    address.SetBase("10.1.1.0", "255.255.255.0");
    Ipv4InterfaceContainer interfaces = address.Assign(devices);

    // Print node info
    for (uint32_t i = 0; i < nodes.GetN(); i++)
    {
        Ptr<MobilityModel> mob = nodes.Get(i)->GetObject<MobilityModel>();
        Vector pos = mob->GetPosition();
        NS_LOG_INFO("Node " << i << ": IP=" << interfaces.GetAddress(i) << " Position=(" << pos.x << ", " << pos.y << ")");
    }

    // Create communication: Node 0 sends to Node 2, Node 1 sends to Node 3
    TypeId tid = TypeId::LookupByName("ns3::UdpSocketFactory");

    // Node 0 -> Node 2
    Ptr<Socket> recvSink2 = Socket::CreateSocket(nodes.Get(2), tid);
    InetSocketAddress local2 = InetSocketAddress(Ipv4Address::GetAny(), 9);
    recvSink2->Bind(local2);
    
    Ptr<Socket> source0 = Socket::CreateSocket(nodes.Get(0), tid);
    InetSocketAddress remote2 = InetSocketAddress(interfaces.GetAddress(2), 9);
    source0->Connect(remote2);
    
    // Node 1 -> Node 3
    Ptr<Socket> recvSink3 = Socket::CreateSocket(nodes.Get(3), tid);
    InetSocketAddress local3 = InetSocketAddress(Ipv4Address::GetAny(), 10);
    recvSink3->Bind(local3);
    
    Ptr<Socket> source1 = Socket::CreateSocket(nodes.Get(1), tid);
    InetSocketAddress remote3 = InetSocketAddress(interfaces.GetAddress(3), 10);
    source1->Connect(remote3);

    // Start traffic generation
    Simulator::Schedule(Seconds(2.0), &GenerateTraffic, source0, 512, 5, Seconds(0.5));
    Simulator::Schedule(Seconds(2.5), &GenerateTraffic, source1, 512, 5, Seconds(0.5));
    
    // Enable routing protocol
    Ipv4GlobalRoutingHelper::PopulateRoutingTables();
    
    Simulator::Stop(Seconds(10.0));
    Simulator::Run();
    Simulator::Destroy();
    return 0;

}