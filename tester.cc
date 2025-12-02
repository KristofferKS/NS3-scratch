/*
 * SPDX-License-Identifier: GPL-2.0-only

 * This is the first example script to learn ns-3
 *     first.cc       from examples/tutorial/first.cc
 */

#include "ns3/applications-module.h"
#include "ns3/core-module.h"
#include "ns3/internet-module.h"
#include "ns3/network-module.h"
#include "ns3/point-to-point-module.h"

// Default Network Topology
//
//       10.1.1.0
// n0 -------------- n1
//    point-to-point
//

// using the ns3 namespace to avoid prefixing all classes with ns3::
using namespace ns3;

// Define logging component
NS_LOG_COMPONENT_DEFINE("FirstScriptExample");

// main function - execution starts here
int
main(int argc, char* argv[])
{
    // Enable command line arguments - much like argparse in python
    CommandLine cmd(__FILE__);
    cmd.Parse(argc, argv);

    // Set time resolution to nanoseconds
    Time::SetResolution(Time::NS);

    // Enable logging for UdpEchoClient and UdpEchoServer applications
    LogComponentEnable("UdpEchoClientApplication", LOG_LEVEL_INFO);
    LogComponentEnable("UdpEchoServerApplication", LOG_LEVEL_INFO);

    // Create two nodes, NodeContainer is a helper class, nodes are an object - we create two
    // we can later call nodes.Get(0) and nodes.Get(1) to get the two nodes
    // we can later install applications, netdevices, channels etc on these nodes
    NodeContainer nodes;
    nodes.Create(2);

    // making pointToPoint object from PointToPointHelper class
    // this is the channel that will connect the two nodes
    // setting attributes for the point to point connection
    PointToPointHelper pointToPoint;
    pointToPoint.SetDeviceAttribute("DataRate", StringValue("5Mbps"));
    pointToPoint.SetChannelAttribute("Delay", StringValue("2ms"));

    // installing the point to point channel on the nodes
    NetDeviceContainer devices;                     // creating a NetDeviceContainer object: devices, to hold the netdevices
    devices = pointToPoint.Install(nodes);
    
    // installing the internet stack protocols (tcp, udp, ip etc) on the nodes
    // this is a large piece to install if not all is needed - but for simplicity we install all
    InternetStackHelper stack;
    stack.Install(nodes);

    // assigning IP addresses to the netdevices, 
    // the setbase takes the network address base and the network mask, 
    // we set the address base, and it will then increment automatically for each node
    Ipv4AddressHelper address;
    address.SetBase("10.1.1.0", "255.255.255.0");

    // assigning the addresses to the devices
    Ipv4InterfaceContainer interfaces = address.Assign(devices);

    // creating a UdpEchoServerHelper object to make a server application
    UdpEchoServerHelper echoServer(9);

    // installing the server application on node 1 (the second node)
    ApplicationContainer serverApps = echoServer.Install(nodes.Get(1));
    serverApps.Start(Seconds(1));       // starting at 1 second
    serverApps.Stop(Seconds(10));       // stopping at 10 seconds

    // creating a UdpEchoClientHelper object to make a client application and setting the server address and port (port 9)
    UdpEchoClientHelper echoClient(interfaces.GetAddress(1), 9);
    echoClient.SetAttribute("MaxPackets", UintegerValue(1));        // setting max packets to send
    echoClient.SetAttribute("Interval", TimeValue(Seconds(1)));     // setting interval between packets
    echoClient.SetAttribute("PacketSize", UintegerValue(1024));     // setting packet size

    // installing the client application on node 0 (the first node)
    ApplicationContainer clientApps = echoClient.Install(nodes.Get(0));
    clientApps.Start(Seconds(2));       // starting at 2 seconds
    clientApps.Stop(Seconds(10));       // stopping at 10 seconds

    
    Simulator::Run();           // running the simulation
    Simulator::Destroy();       // destroying the simulation
    return 0;   // return success
}
