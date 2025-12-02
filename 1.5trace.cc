#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("TraceExample");

// Callback function for packet transmission
void TxCallback(Ptr<const Packet> packet)
{
    NS_LOG_INFO("Packet transmitted: " << packet->GetSize() << " bytes");
}

// Callback function for packet reception
void RxCallback(Ptr<const Packet> packet, const Address &from, const Address &local)
{
    NS_LOG_INFO("Packet received: " << packet->GetSize() << " bytes from " 
                << InetSocketAddress::ConvertFrom(from).GetIpv4());
}

void DropCallback(Ptr<const Packet> packet)
{
    NS_LOG_INFO("Packet dropped: " << packet->GetSize() << " bytes");
}

int main(int argc, char *argv[])
{
    LogComponentEnable("TraceExample", LOG_LEVEL_INFO);
    
    // Create nodes
    NodeContainer nodes;
    nodes.Create(2);
    
    // Point-to-point link
    PointToPointHelper pointToPoint;
    pointToPoint.SetDeviceAttribute("DataRate", StringValue("5Mbps"));
    pointToPoint.SetChannelAttribute("Delay", StringValue("2ms"));
    
    NetDeviceContainer devices;
    devices = pointToPoint.Install(nodes);
    
    // Install Internet stack
    InternetStackHelper stack;
    stack.Install(nodes);
    
    // Assign IP addresses
    Ipv4AddressHelper address;
    address.SetBase("10.1.1.0", "255.255.255.0");
    Ipv4InterfaceContainer interfaces = address.Assign(devices);
    
    // Create server
    UdpEchoServerHelper echoServer(9);
    ApplicationContainer serverApps = echoServer.Install(nodes.Get(1));
    serverApps.Start(Seconds(1.0));
    serverApps.Stop(Seconds(10.0));
    
    // Create client
    UdpEchoClientHelper echoClient(interfaces.GetAddress(1), 9);
    echoClient.SetAttribute("MaxPackets", UintegerValue(5));
    echoClient.SetAttribute("Interval", TimeValue(Seconds(1.0)));
    echoClient.SetAttribute("PacketSize", UintegerValue(1024));
    
    ApplicationContainer clientApps = echoClient.Install(nodes.Get(0));
    clientApps.Start(Seconds(2.0));
    clientApps.Stop(Seconds(10.0));
    
    // Connect trace callbacks
    devices.Get(0)->TraceConnectWithoutContext("PhyTxEnd", MakeCallback(&TxCallback));
    serverApps.Get(0)->TraceConnectWithoutContext("RxWithAddresses", MakeCallback(&RxCallback));
    devices.Get(1)->TraceConnectWithoutContext("PhyRxDrop", MakeCallback(&DropCallback));
    devices.Get(0)->TraceConnectWithoutContext("PhyTxDrop", MakeCallback(&DropCallback));
    
    // Enable pcap tracing
    pointToPoint.EnablePcapAll("trace-example");

    Ptr<Packet> p = Create<Packet>(1024); // or whatever size you want
    Simulator::Schedule(Seconds(2.0), &DropCallback, p);

    
    Simulator::Run();
    Simulator::Destroy();
    return 0;
}