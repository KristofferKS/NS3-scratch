
# include "ns3/core-module.h"
# include "ns3/network-module.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("CreateNodes");


int main(int argc, char* argv[])
{
    LogComponentEnable("CreateNodes", LOG_LEVEL_INFO);

    // Create a Nodecontainer to hold 3 Nodes
    NodeContainer nodes;
    NodeContainer computeNodes;
    nodes.Create(5);
    computeNodes.Create(5);

    NS_LOG_INFO("Created " << nodes.GetN() << " nodes");
    NS_LOG_INFO("Created " << computeNodes.GetN() << " compute nodes");

    // Access individual nodes in for loop and print their IDs
    for (int i = 0; i < 5; i++) {
        NS_LOG_INFO("Node " << i << " ID: " << nodes.Get(i)->GetId());
        NS_LOG_INFO("ComputeNode " << i << " ID: " << computeNodes.Get(i)->GetId());
    }

    Simulator::Run();
    Simulator::Destroy();
    return 0;

}