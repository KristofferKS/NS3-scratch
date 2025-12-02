
/*
simple simulation exercises to learn ns3

hello ns3 script
*/

#include "ns3/core-module.h"

using namespace ns3;
NS_LOG_COMPONENT_DEFINE("HelloNS3");

void scheduledFunction()
{
    NS_LOG_INFO("Halfway Done!");
}

int main(int arc, char* argv[])
{
    // print hello world to console

    // enable logging
    
    LogComponentEnable("HelloNS3", LOG_LEVEL_INFO);

    // create simulator instance
    Simulator::Stop(Seconds(20.0));
    Simulator::Schedule(Seconds(2.5), &scheduledFunction);
    Simulator::Run();

    
    

    NS_LOG_INFO("Simulation Complete! Hello, NS-3!");
    Simulator::Destroy();
    return 0;
}