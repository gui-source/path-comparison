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

#include "../libs/random/random.h"

using namespace ns3;
using Random = effolkronium::random_static;

// ---------------------------------------------------- //
// --- BEGIN OF SIMULATION CONFIGURATION PARAMETERS --- //
// ---------------------------------------------------- //

// Random Seed
uint32_t seed = 1;

// ---[TRACES]---
// Store Traces
bool storeTraces = true;

// Location to store traces
std::string tracesPath = "traces/";

// ---[SIMULATION PARAMETERS]---
// Time to initiate the TCP connection
Time startTimeTCP = Seconds(0.1);

// Time to stop the TCP connection
Time stopTimeTCP = Seconds(100);

// Simulation stop time
Time stopTimeSimulation = Seconds(101);

// ---[TCP PARAMETERS] ---
Time intervalTCP = Seconds(1);
uint32_t segmentSize = 512;
uint32_t MTU_bytes = segmentSize + 54;
uint8_t delAckCount = 2;
uint8_t initialCwnd = 4;
std::string delAckTimeout = "200ms";
std::string socketFactory = "ns3::TcpSocketFactory";
std::string qdiscTypeId = "ns3::PFifoFastQueueDisc";
std::string tcpRecovery = "ns3::TcpClassicRecovery";
std::string tcpVariantId = "ns3::TcpNewReno";
bool enableSack = false;
double minRTO = 1.0;

// ---[TOPOLOGY PARAMETERS]---
std::string flowDataRateStr = "10Gbps";
std::string delay_serialization = "0ms";

std::string bandwidth_access = "10Gbps";
std::string delay_access = "5ms";

std::string long_delay_access = "200ms";

double lambda_1 = 1000.00;
double lambda_2 = 1000.00;

// ---[POINTER TO THE DEVICE THAT WILL IMPLEMENT PACKET DROPING]
NetDeviceContainer * netDeviceToDropPacket = NULL;


// -------------------------------------------------- //
// --- END OF SIMULATION CONFIGURATION PARAMETERS --- //
// -------------------------------------------------- //

// Global Variables
Ptr<PacketSink> sinker;
int packetsDroppedInQueue = 0;
int64_t lastTotalRx = 0;
uint32_t bytes_to_send = 50000;

void InstallPoissonSend(Ptr<Node> src, Ipv4Address destAddress, uint16_t destPort, std::string socketFactory, 
  double lambda, uint32_t packetSize, Time startTime, Time stopTime){
    // Set up On time (just time for 1 packet)
    Ptr<ConstantRandomVariable> onTime = CreateObject<ConstantRandomVariable> ();
    DataRate dataRate(flowDataRateStr);
    double onRate = ((double) packetSize) * 8.0 / dataRate.GetBitRate();
    onTime->SetAttribute("Constant", DoubleValue(onRate));

    std::cout << "On time: " << onRate << std::endl;
    std::cout << "Mean inter-arrival time: " << 1/lambda - onRate << std::endl;

    // Set up Poisson traffic
    Ptr<ExponentialRandomVariable> interArrival = CreateObject<ExponentialRandomVariable> ();
    interArrival->SetAttribute ("Mean", DoubleValue (1/lambda - onRate)); // Set the mean inter-arrival time (lambda)

    OnOffHelper onoff (socketFactory, Address (InetSocketAddress (destAddress, destPort)));
    onoff.SetAttribute ("DataRate", StringValue (flowDataRateStr));
    onoff.SetAttribute ("OnTime", PointerValue (onTime));
    onoff.SetAttribute ("OffTime", PointerValue (interArrival));
    onoff.SetAttribute("PacketSize", UintegerValue(packetSize));

    // Install the Poisson traffic source on src
    ApplicationContainer apps = onoff.Install (src);
    apps.Start (startTime);
    apps.Stop (stopTime);
}

void InstallPacketSink(Ptr<Node> dest, uint16_t port, std::string socketFactory, Time startTime, Time stopTime){
    PacketSinkHelper sink (socketFactory, InetSocketAddress(Ipv4Address::GetAny(), port));
    ApplicationContainer sinkApps = sink.Install(dest);
    sinker = StaticCast<PacketSink>(sinkApps.Get(0));
    sinkApps.Start(startTime);
    sinkApps.Stop(stopTime);
}


void enableLinkLoss(NetDeviceContainer ndToEnableLinkLoss, double percentage){
    // Enable link loss on initial direction of link
    Ptr<PointToPointNetDevice> device = ndToEnableLinkLoss.Get(0)->GetObject<PointToPointNetDevice> ();
    device->EnableLinkLoss(percentage, true, tracesPath);

    // Enable link loss for other direction as well
    device = ndToEnableLinkLoss.Get(1)->GetObject<PointToPointNetDevice> ();
    device->EnableLinkLoss(percentage, true, tracesPath);
}

int main (int argc, char *argv[]){
    // Command line arguments
    CommandLine cmd;
    cmd.AddValue("tcpVariantId", "TCP variant", tcpVariantId);
    // cmd.AddValue("enableSack", "Enable/disable sack in TCP", enableSack);
    cmd.AddValue("seed", "The random seed", seed);
    cmd.AddValue("simStopTime", "The simulation stop time", stopTimeSimulation);
    cmd.AddValue("initialCwnd", "Initial CWND window", initialCwnd);
    cmd.Parse(argc, argv);

    // Set Random Seed
    Random::seed(seed);

    // TCP Recovery Algorithm
    Config::SetDefault("ns3::TcpL4Protocol::RecoveryType", TypeIdValue(TypeId::LookupByName(tcpRecovery)));

    // Set Congestion Control Algorithm
    Config::SetDefault("ns3::TcpL4Protocol::SocketType", StringValue(tcpVariantId));

    Config::SetDefault("ns3::TcpSocket::DelAckTimeout", TimeValue(Time(delAckTimeout)));
    Config::SetDefault("ns3::TcpSocket::SndBufSize", UintegerValue(1073741824));
    Config::SetDefault("ns3::TcpSocket::RcvBufSize", UintegerValue(1073741824));

    // Set default initial congestion window
    Config::SetDefault("ns3::TcpSocket::InitialCwnd", UintegerValue(initialCwnd));

    // Set default delayed ack count to a specified value
    Config::SetDefault("ns3::TcpSocket::DelAckCount", UintegerValue(delAckCount));

    // Set default segment size of TCP packet to a specified value
    Config::SetDefault("ns3::TcpSocket::SegmentSize", UintegerValue(segmentSize));

    // Enable/Disable SACK in TCP
    Config::SetDefault("ns3::TcpSocketBase::Sack", BooleanValue(enableSack));

    Config::SetDefault("ns3::TcpSocketBase::MinRto", TimeValue(Seconds(minRTO)));

    // Initialise nodes
    NodeContainer leftNodes, rightNodes, routers;
    routers.Create(5);
    leftNodes.Create(2);
    rightNodes.Create(1);

    Names::Add("Router1", routers.Get(0));
    Names::Add("Router2", routers.Get(1));
    Names::Add("Router3", routers.Get(2));
    Names::Add("Router4", routers.Get(3));
    Names::Add("Router5", routers.Get(4));
    Names::Add("Client1", leftNodes.Get(0));
    Names::Add("Client2", leftNodes.Get(1));
    Names::Add("Server", rightNodes.Get(0));

    DataRate b_access(bandwidth_access);
    Time d_access(delay_access);
    Time d_serialization(delay_serialization);

    // Create the point-to-point link helpers and connect two router nodes
    PointToPointHelper pointToPointRouter;
    pointToPointRouter.SetDeviceAttribute("DataRate", StringValue(bandwidth_access));
    pointToPointRouter.SetChannelAttribute("Delay", StringValue(delay_access));


    PointToPointHelper pointToPointRouterWithDelay;
    pointToPointRouterWithDelay.SetDeviceAttribute("DataRate", StringValue(bandwidth_access));
    pointToPointRouterWithDelay.SetChannelAttribute("Delay", StringValue(long_delay_access ));

    NetDeviceContainer r1r2ND = pointToPointRouter.Install(routers.Get(0), routers.Get(1));
    Names::Add("IntR1->R2", r1r2ND.Get(0));
    Names::Add("IntR2->R1", r1r2ND.Get(1));

    NetDeviceContainer r2r3ND = pointToPointRouter.Install(routers.Get(1), routers.Get(2));
    Names::Add("IntR2->R3", r2r3ND.Get(0));
    Names::Add("IntR3->R2", r2r3ND.Get(1));
    NetDeviceContainer r2r4ND = pointToPointRouter.Install(routers.Get(1), routers.Get(3));
    Names::Add("IntR2->R4", r2r4ND.Get(0));
    Names::Add("IntR4->R2", r2r4ND.Get(1));

    NetDeviceContainer r3r5ND = pointToPointRouter.Install(routers.Get(2), routers.Get(4));
    Names::Add("IntR3->R5", r3r5ND.Get(0));
    Names::Add("IntR5->R3", r3r5ND.Get(1));
    NetDeviceContainer r4r5ND = pointToPointRouter.Install(routers.Get(3), routers.Get(4));
    Names::Add("IntR4->R5", r4r5ND.Get(0));
    Names::Add("IntR5->R4", r4r5ND.Get(1));

    // Enable link loss 
    enableLinkLoss(r3r5ND, 0.02);
    enableLinkLoss(r4r5ND, 0.01);

    NetDeviceContainer src1r1ND = pointToPointRouter.Install(leftNodes.Get(0), routers.Get(0));
    Names::Add("IntCl1->R1", src1r1ND.Get(0));
    Names::Add("IntR1->Cl1", src1r1ND.Get(1));
    NetDeviceContainer src2r1ND = pointToPointRouter.Install(leftNodes.Get(1), routers.Get(0));
    Names::Add("IntCl2->R1", src2r1ND.Get(0));
    Names::Add("IntR1->Cl2", src2r1ND.Get(1));
    NetDeviceContainer r5destND = pointToPointRouter.Install(routers.Get(4), rightNodes.Get(0));
    Names::Add("IntR5->Se", r5destND.Get(0));
    Names::Add("IntSe->R5", r5destND.Get(1));

    InternetStackHelper internetStack;
    internetStack.Install(leftNodes);
    internetStack.Install(rightNodes);
    internetStack.Install(routers);

    Ipv4AddressHelper ipAddresses("10.0.0.0", "255.255.255.0");
    Ipv4InterfaceContainer r1r2IPAddress = ipAddresses.Assign(r1r2ND);
    ipAddresses.NewNetwork();
    Ipv4InterfaceContainer r2r3IPAddress = ipAddresses.Assign(r2r3ND);
    ipAddresses.NewNetwork();
    Ipv4InterfaceContainer r2r4IPAddress = ipAddresses.Assign(r2r4ND);
    ipAddresses.NewNetwork();
    Ipv4InterfaceContainer r3r5IPAddress = ipAddresses.Assign(r3r5ND);
    ipAddresses.NewNetwork();
    Ipv4InterfaceContainer r4r5IPAddress = ipAddresses.Assign(r4r5ND);
    ipAddresses.NewNetwork();
    Ipv4InterfaceContainer src1r1IPAddress = ipAddresses.Assign(src1r1ND);
    ipAddresses.NewNetwork();
    Ipv4InterfaceContainer src2r1IPAddress = ipAddresses.Assign(src2r1ND);
    ipAddresses.NewNetwork();
    Ipv4InterfaceContainer r5destIPAddress = ipAddresses.Assign(r5destND);
    ipAddresses.NewNetwork(); 

    std::cout << "From src 1 : " << src1r1IPAddress.GetAddress(0) << " to dest " << r5destIPAddress.GetAddress(1) << std::endl;

    std::cout << "From src 2 : " << src2r1IPAddress.GetAddress(0) << " to dest " << r5destIPAddress.GetAddress(1) << std::endl;

    Ipv4GlobalRoutingHelper::PopulateRoutingTables();

    uint16_t server_port = 50000;

    /* Install packet sink at receiver side */  
    InstallPacketSink(rightNodes.Get(0), server_port, socketFactory, Seconds(0.01), stopTimeSimulation);

    /* Install PoissonSend application */
    InstallPoissonSend(leftNodes.Get(0), r5destIPAddress.GetAddress(1), server_port, socketFactory, lambda_1, segmentSize, Seconds(0.2), stopTimeTCP);
    InstallPoissonSend(leftNodes.Get(1), r5destIPAddress.GetAddress(1), server_port, socketFactory, lambda_1, segmentSize, Seconds(0.2), stopTimeTCP);

    pointToPointRouter.EnableAsciiAll(tracesPath);

    Simulator::Stop(stopTimeSimulation);
    Simulator::Run();
    Simulator::Destroy();
    return 0;
}
