# include "ns3/core-module.h"
# include "ns3/network-module.h"
# include "ns3/internet-module.h"
# include "ns3/point-to-point-module.h"
# include "ns3/applications-module.h"

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("TraceExample");

// Callback function for packet transmission
void TxCallback(Ptr<const Packet> packet)
{
    NS_LOG_INFO("Packet transmitted: " << packet->GetSize() >> " bytes")
}

// Callback function for packet reception
void RxCallback(Ptr<const Packet> packet, const address &address)
{
    NS_LOG_INFO("Packet received: " << packet->GetSize() >> " bytes from " << IKnetSocketAddress::ConvertFrom(Address).getIpv4());
}

int main(int argc, char *argv[])
{
    
}