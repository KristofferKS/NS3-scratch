#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/netanim-module.h"
#include <map>
#include <set>
#include <fstream>
#include <iomanip>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("MixedMobility");

// Global variables for tracking connections
std::map<uint32_t, std::map<double, std::set<uint32_t>>> connectionHistory;
AnimationInterface* globalAnim = nullptr;
NodeContainer globalStaticNodes;
NodeContainer globalMobileNodes;
Ipv4InterfaceContainer globalInterfaces;
double commRange = 50.0;
double simTime = 300.0;
uint32_t nMobileNodes = 5;
std::ofstream connectionFile;

void CheckConnections()
{
    double currentTime = Simulator::Now().GetSeconds();
    
    for (uint32_t m = 0; m < globalMobileNodes.GetN(); m++)
    {
        Ptr<MobilityModel> mobileMob = globalMobileNodes.Get(m)->GetObject<MobilityModel>();
        Vector mobilePos = mobileMob->GetPosition();
        
        bool inRange = false;
        std::set<uint32_t> connectedStatic;
        
        // Find closest static node within range
        double minDistance = commRange + 1;
        
        for (uint32_t s = 0; s < globalStaticNodes.GetN(); s++)
        {
            Ptr<MobilityModel> staticMob = globalStaticNodes.Get(s)->GetObject<MobilityModel>();
            Vector staticPos = staticMob->GetPosition();
            
            double distance = CalculateDistance(mobilePos, staticPos);
            
            if (distance <= commRange)
            {
                inRange = true;
                connectedStatic.insert(s);
                
                // Track the closest one
                if (distance < minDistance)
                {
                    minDistance = distance;
                }
            }
        }
        
        if (!connectedStatic.empty())
        {
            connectionHistory[m][currentTime] = connectedStatic;
        }
        
        // Update mobile node color
        if (inRange)
        {
            globalAnim->UpdateNodeColor(globalMobileNodes.Get(m), 0, 255, 0); // Green
        }
        else
        {
            globalAnim->UpdateNodeColor(globalMobileNodes.Get(m), 255, 0, 0); // Red
        }
    }
    
    if (Simulator::Now().GetSeconds() < simTime - 1.0)
    {
        Simulator::Schedule(Seconds(1.0), &CheckConnections);
    }
}

void PrintNodeStatus(NodeContainer staticNodes, NodeContainer mobileNodes)
{
    for (uint32_t i = 0; i < staticNodes.GetN(); i++)
    {
        Ptr<MobilityModel> mob = staticNodes.Get(i)->GetObject<MobilityModel>();
    }

    for (uint32_t i = 0; i < mobileNodes.GetN(); i++)
    {
        Ptr<MobilityModel> mob = mobileNodes.Get(i)->GetObject<MobilityModel>();
        Ptr<MobilityModel> mobileMob = mobileNodes.Get(i)->GetObject<MobilityModel>();
        Vector mobilePos = mobileMob->GetPosition();
        std::vector<std::pair<uint32_t, double>> connected;
        
        for (uint32_t s = 0; s < staticNodes.GetN(); s++)
        {
            Ptr<MobilityModel> staticMob = staticNodes.Get(s)->GetObject<MobilityModel>();
            Vector staticPos = staticMob->GetPosition();
            double distance = CalculateDistance(mobilePos, staticPos);
            if (distance <= commRange)
            {
                connected.push_back(std::make_pair(s, distance));
            }
        }
    }

    if (Simulator::Now().GetSeconds() < 29)
    {
        Simulator::Schedule(Seconds(5), &PrintNodeStatus, staticNodes, mobileNodes);
    }
}

void PrintFinalSummary()
{
    std::cout << "\n========================================" << std::endl;
    std::cout << "SIMULATION SUMMARY" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Communication Range: " << commRange << " meters" << std::endl;
    std::cout << "Simulation Duration: " << Simulator::Now().GetSeconds() << " seconds" << std::endl;
    std::cout << "\n";

    double totalConnectivityPercentage = 0.0;

    for (uint32_t m = 0; m < globalMobileNodes.GetN(); m++)
    {
        std::cout << "Mobile Node " << m << ":" << std::endl;
        
        if (connectionHistory[m].empty())
        {
            std::cout << "  - Never connected to any static node" << std::endl;
        }
        else
        {
            std::set<uint32_t> uniqueConnections;
            
            for (auto& timePair : connectionHistory[m])
            {
                for (uint32_t staticId : timePair.second)
                {
                    uniqueConnections.insert(staticId);
                }
            }
            
            std::cout << "  - Connected to " << uniqueConnections.size() << " unique static node(s): ";
            std::cout << std::endl;
            std::cout << "  - Number of connection events: " << connectionHistory[m].size() << std::endl;
            double connectivityPercentage = (connectionHistory[m].size() * 1.0 / simTime) * 100.0;
            std::cout << "  - Approximate time in range: " << std::fixed << std::setprecision(1)
                      << connectivityPercentage << "%" << std::endl;
            totalConnectivityPercentage += connectivityPercentage;
        }
        std::cout << std::endl;
    }

    // Calculate and display average connectivity
    double averageConnectivity = totalConnectivityPercentage / globalMobileNodes.GetN();
    std::cout << "========================================" << std::endl;
    std::cout << "AVERAGE CONNECTION TIME (ALL MOBILE NODES): " << std::fixed << std::setprecision(1)
              << averageConnectivity << "%" << std::endl;
    std::cout << "========================================\n" << std::endl;

    std::cout << "Static Node Connection Statistics:" << std::endl;
    for (uint32_t s = 0; s < globalStaticNodes.GetN(); s++)
    {
        std::set<uint32_t> mobileNodesConnected;
        int connectionCount = 0;
        
        for (auto& mobilePair : connectionHistory)
        {
            for (auto& timePair : mobilePair.second)
            {
                if (timePair.second.find(s) != timePair.second.end())
                {
                    mobileNodesConnected.insert(mobilePair.first);
                    connectionCount++;
                }
            }
        }
        
        Ptr<MobilityModel> mob = globalStaticNodes.Get(s)->GetObject<MobilityModel>();
        Vector pos = mob->GetPosition();
        
        std::cout << "Static Node " << s << " at (" << std::fixed << std::setprecision(2) 
                  << pos.x << ", " << pos.y << "):" << std::endl;
        std::cout << "  - Total connection events: " << connectionCount << std::endl;
        std::cout << std::endl;
    }
    
    std::cout << "========================================" << std::endl;
    std::cout << "\nFiles generated:" << std::endl;
    std::cout << "  - connections.log (detailed connection log)" << std::endl;
    std::cout << "  - mixed-mobility.xml (NetAnim file)" << std::endl;
    std::cout << "\nVisualization:" << std::endl;
    std::cout << "  - Blue (large) = Static nodes" << std::endl;
    std::cout << "  - Green (small) = Mobile nodes in range" << std::endl;
    std::cout << "  - Red (small) = Mobile nodes out of range" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

int main(int argc, char *argv[])
{
    LogComponentEnable("MixedMobility", LOG_LEVEL_INFO);

    uint32_t nNodes = 5;
    
    CommandLine cmd;
    // cmd.AddValue("nNodes", "Number of static nodes", nNodes);
    cmd.AddValue("range", "Communication range in meters", commRange);
    cmd.AddValue("simTime", "Simulation time in seconds", simTime);
    cmd.AddValue("nMobileNodes", "Number of mobile nodes", nMobileNodes);
    cmd.Parse(argc, argv);

    connectionFile.open("connections.log");
    connectionFile << "Connection Log - Communication Range: " << commRange << " meters" << std::endl;
    connectionFile << "==========================================================" << std::endl << std::endl;

    NodeContainer staticNodes;
    staticNodes.Create(nNodes);

    NodeContainer mobileNodes;
    mobileNodes.Create(nMobileNodes);

    globalStaticNodes = staticNodes;
    globalMobileNodes = mobileNodes;

    NodeContainer allNodes;
    allNodes.Add(staticNodes);
    allNodes.Add(mobileNodes);

    // WiFi setup
    YansWifiChannelHelper channel = YansWifiChannelHelper::Default();
    YansWifiPhyHelper phy;
    phy.SetChannel(channel.Create());

    WifiHelper wifi;
    wifi.SetStandard(WIFI_STANDARD_80211g);
    wifi.SetRemoteStationManager("ns3::AarfWifiManager");

    WifiMacHelper mac;
    mac.SetType("ns3::AdhocWifiMac");

    NetDeviceContainer devices = wifi.Install(phy, mac, allNodes);

    // Static node mobility
    MobilityHelper staticMobility;
    Ptr<ListPositionAllocator> staticPos = CreateObject<ListPositionAllocator>();
    staticPos->Add(Vector(30,19.125,0));
    staticPos->Add(Vector(70,19.125,0));
    staticPos->Add(Vector(17.64,57.165,0));
    staticPos->Add(Vector(50,80.875,0));
    staticPos->Add(Vector(82.36,57.165,0));

    staticMobility.SetPositionAllocator(staticPos);
    staticMobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
    staticMobility.Install(staticNodes);

    // Mobile node mobility
    MobilityHelper mobileMobility;
    Ptr<RandomRectanglePositionAllocator> waypointAllocator = CreateObject<RandomRectanglePositionAllocator>();
    waypointAllocator->SetAttribute("X", StringValue("ns3::UniformRandomVariable[Min=0|Max=100.0]"));
    waypointAllocator->SetAttribute("Y", StringValue("ns3::UniformRandomVariable[Min=0|Max=100.0]"));
    
    mobileMobility.SetMobilityModel("ns3::RandomWaypointMobilityModel",
                              "Speed", StringValue("ns3::UniformRandomVariable[Min=1.0|Max=5.0]"),
                              "Pause", StringValue("ns3::UniformRandomVariable[Min=1.0|Max=5.0]"),
                              "PositionAllocator", PointerValue(waypointAllocator));
    
    mobileMobility.SetPositionAllocator("ns3::RandomRectanglePositionAllocator",
                                  "X", StringValue("ns3::UniformRandomVariable[Min=0|Max=100.0]"),
                                  "Y", StringValue("ns3::UniformRandomVariable[Min=0|Max=100.0]"));
    
    mobileMobility.Install(mobileNodes);

    // Internet stack
    InternetStackHelper stack;
    stack.Install(allNodes);

    Ipv4AddressHelper address;
    address.SetBase("10.0.0.0", "255.255.0.0");
    globalInterfaces = address.Assign(devices);

    // Animation setup
    AnimationInterface anim("mixed-mobility.xml");
    globalAnim = &anim;

    // Configure node appearance
    for (uint32_t i = 0; i < staticNodes.GetN(); i++)
    {
        anim.UpdateNodeDescription(staticNodes.Get(i), "Static-" + std::to_string(i));
        anim.UpdateNodeColor(staticNodes.Get(i), 0, 0, 255); // Blue
        anim.UpdateNodeSize(staticNodes.Get(i)->GetId(), 5, 5);
    }

    for (uint32_t i = 0; i < mobileNodes.GetN(); i++)
    {
        anim.UpdateNodeDescription(mobileNodes.Get(i), "Mobile-" + std::to_string(i));
        anim.UpdateNodeColor(mobileNodes.Get(i), 255, 0, 0); // Red initially
        anim.UpdateNodeSize(mobileNodes.Get(i)->GetId(), 1, 1);
    }

    // Schedule functions
    Simulator::Schedule(Seconds(0), &PrintNodeStatus, staticNodes, mobileNodes);
    Simulator::Schedule(Seconds(0.5), &CheckConnections);

    Simulator::Stop(Seconds(simTime));
    Simulator::Run();
    
    connectionFile.close();
    PrintFinalSummary();
    
    Simulator::Destroy();
    return 0;
}