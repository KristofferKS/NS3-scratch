#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/applications-module.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("WifiAP");

int main(int argc, char *argv[])
{
    Time::SetResolution(Time::NS);
    LogComponentEnable("UdpEchoClientApplication", LOG_LEVEL_INFO);
    LogComponentEnable("UdpEchoServerApplication", LOG_LEVEL_INFO);
    LogComponentEnable("WifiAP", LOG_LEVEL_INFO);

    uint32_t nStations = 5;

    // Crete AP Node and Station nodes
    NodeContainer apNode;
    apNode.Create(1);

    NodeContainer staNodes;
    staNodes.Create(nStations);

    // WiFi Steup
    YansWifiChannelHelper channel = YansWifiChannelHelper::Default();
    YansWifiPhyHelper phy;
    phy.SetChannel(channel.Create());

    WifiHelper wifi;
    wifi.SetStandard(WIFI_STANDARD_80211ax);
    wifi.SetRemoteStationManager("ns3::IdealWifiManager");

    // Configure stations
    WifiMacHelper mac;
    Ssid ssid = Ssid("My-network");
    mac.SetType("ns3::StaWifiMac",
                "Ssid", SsidValue(ssid),
                "ActiveProbing", BooleanValue(false));

    NetDeviceContainer staDevices;
    staDevices = wifi.Install(phy, mac, staNodes);

    // Configure AP
    mac.SetType("ns3::ApWifiMac",
                "Ssid", SsidValue(ssid));

    NetDeviceContainer apDevices;
    apDevices = wifi.Install(phy, mac, apNode);

    // Mobility - Ap at center
    MobilityHelper mobility;

    // AP Position
    mobility.SetPositionAllocator("ns3::GridPositionAllocator",
                                  "MinX", DoubleValue(0.0),
                                  "MinY", DoubleValue(0.0),
                                  "DeltaX", DoubleValue(1.0),
                                  "DeltaY", DoubleValue(1.0),
                                  "GridWidth", UintegerValue(1),
                                  "LayoutType", StringValue("RowFirst"));
    mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    mobility.Install(apNode);

    // Stations in  Circle around AP
    mobility.SetPositionAllocator("ns3::UniformDiscPositionAllocator",
                                  "X", DoubleValue(0.0),
                                  "Y", DoubleValue(0.0),
                                  "rho", DoubleValue(20.0));
    mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    mobility.Install(staNodes);

    // Install internet stack
    InternetStackHelper stack;
    stack.Install(apNode);
    stack.Install(staNodes);

    // Assign IP Addresses
    Ipv4AddressHelper address;
    address.SetBase("10.1.1.0", "255.255.255.0");

    Ipv4InterfaceContainer apInterface;
    apInterface = address.Assign(apDevices);

    Ipv4InterfaceContainer staInterfaces;
    staInterfaces = address.Assign(staDevices);

    // Create an UDP Echo Server on node 0
    UdpEchoServerHelper echoServer(9);       // on port 9
    ApplicationContainer serverApps = echoServer.Install(apNode.Get(0));
    serverApps.Start(Seconds(1.0));
    serverApps.Stop(Seconds(10.0));

    // Create an UDP Echo Client on node 1
    UdpEchoClientHelper echoClient(apInterface.GetAddress(0), 9);       // client connect to server (ip, port)
    echoClient.SetAttribute("MaxPackets", UintegerValue(10));
    echoClient.SetAttribute("Interval", TimeValue(Seconds(0.5)));
    echoClient.SetAttribute("PacketSize", UintegerValue(512));
    
    for (uint32_t i = 0; i < staNodes.GetN(); i++)
    {
        ApplicationContainer clientApps = echoClient.Install(staNodes.Get(i));
        clientApps.Start(Seconds(3.0));
        clientApps.Stop(Seconds(10.0));
    }
    
    // Print COnfiguration
    NS_LOG_INFO("Access Point IP" << apInterface.GetAddress(0));
    for (uint32_t i = 0; i < staNodes.GetN(); i++)
    {
        Ptr<MobilityModel> mob = staNodes.Get(i)->GetObject<MobilityModel>();
        Vector pos = mob->GetPosition();
        NS_LOG_INFO("Station " << i << " IP: " << staInterfaces.GetAddress(i)
                    << " Position: (" << pos.x << ", " << pos.y << ")");
    }

    Simulator::Stop(Seconds(10.0));
    Simulator::Run();
    Simulator::Destroy();
    return 0;
}
