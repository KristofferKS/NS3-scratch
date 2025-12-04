#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/applications-module.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("WifiBasic");

int main(int argc, char *argv[])
{
    LogComponentEnable("WifiBasic", LOG_LEVEL_INFO);
    
    // Create 4 WiFi nodes
    NodeContainer wifiNodes;
    wifiNodes.Create(6);
    
    // Create WiFi channel
    YansWifiChannelHelper channel = YansWifiChannelHelper::Default();
    YansWifiPhyHelper phy;
    phy.SetChannel(channel.Create());
    
    // Configure WiFi
    WifiHelper wifi;
    wifi.SetStandard(WIFI_STANDARD_80211ax);
    wifi.SetRemoteStationManager("ns3::IdealWifiManager");
    
    WifiMacHelper mac;
    Ssid ssid = Ssid("ns3-wifi");
    mac.SetType("ns3::StaWifiMac",
                "Ssid", SsidValue(ssid),
                "ActiveProbing", BooleanValue(false));
    
    NetDeviceContainer staDevices;
    staDevices = wifi.Install(phy, mac, wifiNodes);
    
    // Set up positions (static for now)
    MobilityHelper mobility;
    mobility.SetPositionAllocator("ns3::GridPositionAllocator",
                                  "MinX", DoubleValue(0.0),
                                  "MinY", DoubleValue(0.0),
                                  "DeltaX", DoubleValue(10.0),
                                  "DeltaY", DoubleValue(20.0),
                                  "GridWidth", UintegerValue(3),
                                  "LayoutType", StringValue("RowFirst"));
    
    mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    mobility.Install(wifiNodes);
    
    // Install Internet stack
    InternetStackHelper stack;
    stack.Install(wifiNodes);
    
    // Assign IP addresses
    Ipv4AddressHelper address;
    address.SetBase("10.1.1.0", "255.255.255.0");
    Ipv4InterfaceContainer interfaces = address.Assign(staDevices);

    Ptr<Node> node = wifiNodes.Get(0);
    Ptr<MobilityModel> mobility0 = node->GetObject<MobilityModel>();
    
    // Print node positions and IPs
    for (uint32_t i = 0; i < wifiNodes.GetN(); i++)
    {
        Ptr<Node> node = wifiNodes.Get(i);
        Ptr<MobilityModel> mobility = node->GetObject<MobilityModel>();
        Vector pos = mobility->GetPosition();
        double distance = mobility0->GetDistanceFrom(mobility);
        
        NS_LOG_INFO("Node " << i << ": IP=" << interfaces.GetAddress(i) 
                    << " Position=(" << pos.x << ", " << pos.y << ", " << pos.z << ") and distance to Node0: " << distance);
    }
    
    Simulator::Stop(Seconds(10.0));
    Simulator::Run();
    Simulator::Destroy();
    return 0;
}