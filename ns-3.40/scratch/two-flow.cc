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
#include "ns3/rng-seed-manager.h"

using namespace ns3;
// using Random = effolkronium::random_static;

// Random Seed
uint32_t seed = 1;

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
Time stopTimeUDP = Seconds(100);

// Simulation stop time
Time stopTimeSimulation = Seconds(101);

// ---[UDP PARAMETERS] ---
std::string flowDataRateStr = "10Gbps";
uint32_t maxPacketSize = 1024;

// ---[TOPOLOGY PARAMETERS]---
std::string bandwidth_access = "10Gbps";
std::string delay_access = "5ms";
std::string delay_serialization = "0ms";
std::string interframe_gap = "0ms";
double lambda_1 = 1000.00;
double lambda_2 = 1000.00;

// ---[POINTER TO THE DEVICE THAT WILL IMPLEMENT PACKET DROPING]
NetDeviceContainer * netDeviceToDropPacket = NULL;


// -------------------------------------------------- //
// --- END OF SIMULATION CONFIGURATION PARAMETERS --- //
// -------------------------------------------------- //

/*
    lambda should be the expected number of packets per second in the Poisson process
*/
void InstallSourceSink(Ptr<Node> src, Ipv4Address destAddress, double lambda){
    // Set up On time (just time for 1 packet)
    Ptr<ConstantRandomVariable> onTime = CreateObject<ConstantRandomVariable> ();
    DataRate dataRate(flowDataRateStr);
    double onRate = ((double) maxPacketSize) * 8.0 / dataRate.GetBitRate();
    onTime->SetAttribute("Constant", DoubleValue(onRate));

    std::cout << "On time: " << onRate << std::endl;
    std::cout << "Mean inter-arrival time: " << 1/lambda - onRate << std::endl;

    // Set up Poisson traffic
    Ptr<ExponentialRandomVariable> interArrival = CreateObject<ExponentialRandomVariable> ();
    interArrival->SetAttribute ("Mean", DoubleValue (1/lambda - onRate)); // Set the mean inter-arrival time (lambda)

    OnOffHelper onoff ("ns3::UdpSocketFactory", Address (InetSocketAddress (destAddress, 9)));
    onoff.SetAttribute ("DataRate", StringValue (flowDataRateStr));
    onoff.SetAttribute ("OnTime", PointerValue (onTime));
    onoff.SetAttribute ("OffTime", PointerValue (interArrival));
    onoff.SetAttribute("PacketSize", UintegerValue(maxPacketSize));

    // Install the Poisson traffic source on src
    ApplicationContainer apps = onoff.Install (src);
    apps.Start (Seconds (0.2));
    apps.Stop (stopTimeUDP);
}

int main (int argc, char *argv[]){
    // Command line arguments
    CommandLine cmd;
    cmd.AddValue("seed", "The random seed", seed);
    cmd.AddValue("delaySerialization", "Delay for serialization for all links", delay_serialization);
    cmd.AddValue("delayAccess", "Delay time for all links", delay_access);
    cmd.AddValue("lambda1", "Parameter for Poisson process traffic, average number of packets a second", lambda_1);
    cmd.AddValue("lambda2", "Parameter for Poisson process traffic, average number of packets a second", lambda_2);
    cmd.Parse(argc, argv);

    RngSeedManager::SetSeed(seed);

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

    Ipv4GlobalRoutingHelper::PopulateRoutingTables();

    InstallSourceSink(leftNodes.Get(0), routerToRight1IPAddress.GetAddress(1), lambda_1);
    InstallSourceSink(leftNodes.Get(0), routerToRight2IPAddress.GetAddress(1), lambda_1);

    /* Trace ASCII */
    pointToPointRouter.EnableAsciiAll(tracesPath);

    Simulator::Stop(stopTimeSimulation);
    Simulator::Run();
    Simulator::Destroy();
    return 0;
}