#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <typeinfo>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/traffic-control-module.h"

using namespace ns3;
// using Random = effolkronium::random_static;

// ---------------------------------------------------- //
// --- BEGIN OF SIMULATION CONFIGURATION PARAMETERS --- //
// ---------------------------------------------------- //

// ---[TRACES]---
// Store Traces
bool storeTraces = true;

// Location to store traces
std::string tracesPath = "traces/";

// ---[SIMULATION PARAMETERS]---
// Time to initiate the UDP connection
Time startTimeUDP = Seconds(0.1);

// Time to stop the UDP connection
Time stopTimeUDP = Seconds(5);

// Simulation stop time
Time stopTimeSimulation = Seconds(60);

// UDP flow interval
Time intervalUDP = Seconds(1);

// ---[UDP PARAMETERS] ---
uint32_t segmentSize = 512;
uint32_t MTU_bytes = segmentSize + 54;
uint8_t delAckCount = 2;
uint8_t initialCwnd = 4;
std::string delAckTimeout = "200ms";
std::string socketFactory = "ns3::TcpSocketFactory";
std::string qdiscTypeId = "ns3::PFifoFastQueueDisc";
double minRTO = 1.0;

// ---[TOPOLOGY PARAMETERS]---
std::string bandwidth_bottleneck = "10Mbps";
std::string bandwidth_access = "20Mbps";
std::string delay_access = "5ms";
std::string delay_serialization = "1.9ms";
std::string interframe_gap = "0ms";

// ---[POINTER TO THE DEVICE THAT WILL IMPLEMENT PACKET DROPING]
NetDeviceContainer * netDeviceToDropPacket = NULL;


// -------------------------------------------------- //
// --- END OF SIMULATION CONFIGURATION PARAMETERS --- //
// -------------------------------------------------- //

// Global Variables
Ptr<PacketSink> sinker;

int packetsDroppedInQueue = 0;
int64_t lastTotalRx = 0;
uint32_t bytes_to_send = 1000;

uint32_t cnt_packets = 0;

void TraceCwnd(uint32_t node, uint32_t cwndWindow, Callback <void, uint32_t, uint32_t> CwndTrace){
    Config::ConnectWithoutContext("/NodeList/" + std::to_string(node) + "/$ns3::TcpL4Protocol/SocketList/" + std::to_string(cwndWindow) + "/CongestionWindow", CwndTrace);
}

static void CwndChange(uint32_t oldCwnd, uint32_t newCwnd){
  std::ofstream fPlotQueue(tracesPath + "cwnd.txt", std::ios::out | std::ios::app);
  fPlotQueue << Simulator::Now().GetSeconds() << " " << newCwnd / segmentSize << " " << newCwnd << std::endl;
  fPlotQueue.close();
}

void InstallBulkSend(Ptr<Node> node, Ipv4Address address, uint16_t port, std::string socketFactory,
  uint32_t nodeId, uint32_t cwndWindow, Callback <void, uint32_t, uint32_t> CwndTrace, uint32_t maxBytesToSend, Time startTime){
    BulkSendHelper source(socketFactory, InetSocketAddress(address, port));
    source.SetAttribute("MaxBytes", UintegerValue(maxBytesToSend));
    ApplicationContainer sourceApps = source.Install(node);
    sourceApps.Start(startTime);
    sourceApps.Stop (stopTimeUDP);
        /* Trace the CWND */
    if(storeTraces == false){
        Simulator::Schedule(startTime + Seconds(0.001), &TraceCwnd, nodeId, cwndWindow, CwndTrace);
    }
}

void InstallPacketSink(Ptr<Node> node, uint16_t port, std::string socketFactory, Time startTime, Time stopTime){
    PacketSinkHelper sink (socketFactory, InetSocketAddress(Ipv4Address::GetAny(), port));
    ApplicationContainer sinkApps = sink.Install(node);
    sinker = StaticCast<PacketSink>(sinkApps.Get(0));
    sinkApps.Start(startTime);
    sinkApps.Stop(stopTime);
}

void RxDrop(Ptr<const Packet> packet) {
    // Print Source and Destination IP Addresses
    /* 
     * Need to copy packet since headers need to be removed
     * to be inspected. Alternatively, remove the headers,
     * and add them back.
     */
    Ptr<Packet> copy = packet->Copy();

    // Headers must be removed in the order they're present.
    PppHeader pppHeader;
    copy->RemoveHeader(pppHeader);
    Ipv4Header ipHeader;
    copy->RemoveHeader(ipHeader);

    std::cout << "Source IP: ";
    ipHeader.GetSource().Print(std::cout);
    std::cout << std::endl;

    std::cout << "Destination IP: ";
    ipHeader.GetDestination().Print(std::cout);
    std::cout << std::endl;
}

void ExaminePacket(Ptr< const Packet > packet){
    // Todo: repurpose for UDP packets
    // Extract UDP Header from the packet
    TcpHeader * tcpHeader = (TcpHeader *) packet->extractTcpHeader();
    uint32_t payloadSize = packet->GetPayloadSize();
    int dPort = int(tcpHeader->GetDestinationPort());

    // Extract the SEQ and ACK numbers
    uint32_t seq = tcpHeader->GetSequenceNumber().GetValue();
    uint32_t ack = tcpHeader->GetAckNumber().GetValue();

    std::cout<<"Time: "<< Simulator::Now() << "[UDP PACKET] [SEQ: " << seq << "] [ACK: " << ack << "] [Payload Length: " << payloadSize << "] [Dest port:" << dPort << "] " << std::endl;
}

int main (int argc, char *argv[]){
    // Command line arguments
    CommandLine cmd;
    cmd.AddValue("bytesToSend", "Number of bytes to send", bytes_to_send);
    cmd.AddValue("delaySerialization", "Delay for serialization for all links", delay_serialization);
    cmd.AddValue("delayAccess", "Delay time for all links", delay_access);
    cmd.Parse(argc, argv);

    // Initialise nodes
    NodeContainer leftNodes, rightNodes, routers;
    routers.Create(1);
    leftNodes.Create(1);
    rightNodes.Create(2);

    Names::Add("Router1", routers.Get(0));
    Names::Add("Client", leftNodes.Get(0));
    Names::Add("Server1", rightNodes.Get(0));
    Names::Add("Server2", rightNodes.Get(1));

    // Create the point-to-point link helpers and connect two router nodes
    PointToPointHelper pointToPointRouter;
    pointToPointRouter.SetDeviceAttribute("DataRate", StringValue(bandwidth_access));
    pointToPointRouter.SetChannelAttribute("Delay", StringValue(delay_access));

    NetDeviceContainer leftToRouter = pointToPointRouter.Install(leftNodes.Get(0), routers.Get(0));
    NetDeviceContainer routerToRight1 = pointToPointRouter.Install(routers.Get(0), rightNodes.Get(0));
    NetDeviceContainer routerToRight2 = pointToPointRouter.Install(routers.Get(0), rightNodes.Get(1));
    Names::Add("IntCl->R1", leftToRouter.Get(0));
    Names::Add("IntR1->Cl", leftToRouter.Get(1));
    Names::Add("IntR1->Se1", routerToRight1.Get(0));
    Names::Add("IntSe1->R1", routerToRight1.Get(1));
    Names::Add("IntR1->Se2", routerToRight2.Get(0));
    Names::Add("IntSe2->R1", routerToRight2.Get(1));

    InternetStackHelper internetStack;
    internetStack.Install(leftNodes);
    internetStack.Install(rightNodes);
    internetStack.Install(routers);

    Ipv4AddressHelper ipAddresses("10.0.0.0", "255.255.255.0");
    Ipv4InterfaceContainer leftToRouterIPAddress = ipAddresses.Assign(leftToRouter);
    ipAddresses.NewNetwork();
    Ipv4InterfaceContainer routerToRight1IPAddress = ipAddresses.Assign(routerToRight1);
    ipAddresses.NewNetwork();
    Ipv4InterfaceContainer routerToRight2IPAddress = ipAddresses.Assign(routerToRight2);
    ipAddresses.NewNetwork();

    Ipv4GlobalRoutingHelper::PopulateRoutingTables();

    /* Packet Printing is mandatory for the Packet Drop Test */
    Packet::EnablePrinting();

    uint16_t server_port = 50000;
    /* Install packet sink at receiver side */
    InstallPacketSink(rightNodes.Get(0), server_port, socketFactory, Seconds(0.01), stopTimeSimulation);
    InstallPacketSink(rightNodes.Get(1), server_port, socketFactory, Seconds(0.01), stopTimeSimulation);

    /* Install BulkSend application */
    InstallBulkSend(leftNodes.Get(0), routerToRight1IPAddress.GetAddress(1), server_port, socketFactory, 2, 0, MakeCallback(&CwndChange), bytes_to_send, Seconds(0.2));
    InstallBulkSend(leftNodes.Get(0), routerToRight2IPAddress.GetAddress(1), server_port, socketFactory, 2, 0, MakeCallback(&CwndChange), bytes_to_send, Seconds(0.2));
    
    /* Trace ASCII */
    pointToPointRouter.EnableAsciiAll(tracesPath);

    Simulator::Stop(stopTimeSimulation);
    Simulator::Run();
    Simulator::Destroy();
    return 0;
}