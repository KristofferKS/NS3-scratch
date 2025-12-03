# include "ns3/core-module.h"
# include "ns3/network-module.h"
# include "ns3/internet-module.h"
# include "ns3/point-to-point-module.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("P2PSimple");

int main(int argc, char* argv[])
{
    LogComponentEnable("P2PSimple", LOG_LEVEL_INFO);

    // Create two nodes
    NodeContainer nodes;
    nodes.Create(3);

    // Create point to point link
    PointToPointHelper pointToPoint;
    pointToPoint.SetDeviceAttribute("DataRate", StringValue("10Mbps"));
    pointToPoint.SetChannelAttribute("Delay", StringValue("5ms"));

    // link 0-1
    NodeContainer link1(nodes.Get(0), nodes.Get(1));
    NetDeviceContainer devices1;
    devices1 = pointToPoint.Install(link1);

    // link 1-2
    NodeContainer link2(nodes.Get(1), nodes.Get(2));
    NetDeviceContainer devices2;
    devices2 = pointToPoint.Install(link2);

    // Install Internet stack on nodes
    InternetStackHelper stack;
    stack.Install(nodes);

    // Assign IP addresses to devices
    Ipv4AddressHelper address;

    // addresses for link 1
    address.SetBase("192.168.1.0", "255.255.255.0");
    Ipv4InterfaceContainer interfaces1 = address.Assign(devices1);

    // addresses for link 2
    address.SetBase("192.168.2.0", "255.255.255.0");
    Ipv4InterfaceContainer interfaces2 = address.Assign(devices2);

    NS_LOG_INFO("Node 0 IP: " << interfaces1.GetAddress(0));
    NS_LOG_INFO("Node 1 (in subnet 1) IP: " << interfaces1.GetAddress(1));
    NS_LOG_INFO("Node 1 (in subnet 2) IP: " << interfaces2.GetAddress(0));
    NS_LOG_INFO("Node 2 IP: " << interfaces2.GetAddress(1));

    Simulator::Run();
    Simulator::Destroy();
    return 0;

}