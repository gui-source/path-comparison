#include <iostream>
#include <string>
#include <fstream>
#include <sys/stat.h>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <typeinfo>
#include <map>
#include <list>
#include <stdexcept>

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/flow-monitor-module.h"
#include "ns3/point-to-point-module.h"
#include "ns3/applications-module.h"
#include "ns3/traffic-control-module.h"

#include "../libs/random/random.h"

using namespace ns3;
using Random = effolkronium::random_static;

// ---------------------------------------------------- //
// --- BEGIN OF SIMULATION CONFIGURATION PARAMETERS --- //
// ---------------------------------------------------- //

// Random loss
float lossDuration = 0.0;

// Time-bound
bool timeBound = true;
bool targetZeAchieved = false;
double targetZE;
double targetSeconds;

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
uint32_t stopTimeTCPSeconds = 30;

// Time to stop the TCP connection
Time stopTimeTCP = Seconds(stopTimeTCPSeconds);

// Simulation stop time
Time stopTimeSimulation = Seconds(stopTimeTCPSeconds + 1);

// ---[TCP PARAMETERS] ---
Time intervalTCP = Seconds(1);
uint32_t segmentSize = 512;
uint32_t MTU_bytes = segmentSize + 54;
uint8_t delAckCount = 2;
uint8_t initialCwnd = 10;
std::string delAckTimeout = "200ms";
std::string socketFactory = "ns3::TcpSocketFactory";
std::string tcpRecovery = "ns3::TcpClassicRecovery";
std::string tcpVariantId = "ns3::TcpNewReno";
double minRTO = 1.0;
uint32_t tcpFlows = 2;

// ---[TOPOLOGY PARAMETERS]---
std::string srcThroughput1 = "10Mbps";
std::string srcThroughput2 = "20Mbps";
std::string throughput = "40Mbps";
std::string delay = "3ms";

double linkLoss1 = 0.001;
double linkLoss2 = 0.002;

// --- [DATA STRUCTS FOR SYNC] ---
struct FlowSeqId{
    uint32_t flow;
    uint32_t seq;

    bool operator==(const FlowSeqId &o) const {
        return flow == o.flow && seq == o.seq;
    }

    bool operator<(const FlowSeqId &o) const {
        return flow < o.flow || (flow == o.flow && seq < o.seq);
    }
};

struct SyncTrack{
    FlowSeqId id;
    int flow;
    Time syncExpire;
    Time measureExpire;
    bool active;
};

/// --- [GLOBAL VARIABLES FOR SYNC] --- 
// to keep track of possible synchronisation
SyncTrack * syncTrackerPtr = new SyncTrack{{}, 0, Seconds(0), Seconds(0), false};

// to keep track of what packets we should check rtx for 
Ipv4FlowClassifier classifier;
// map from flow seq id of packets to lookout for rtxs
std::map<FlowSeqId, Time> lookoutRTX1;
std::map<FlowSeqId, Time> lookoutRTX2;
// number of measurements for each flow - number of synced packets
uint32_t syncedPkts = 0;
// number of rtx for each flow
uint32_t rtxFlow1 = 0;
uint32_t rtxFlow2 = 0;
// number of pkts for each flow
uint32_t totalFlow1 = 0;
uint32_t totalFlow2 = 0;

// maximum time between packets that may be synced
double delta = 0.0001;
Time t_delta;

// maximum time to stop looking for rtx
Time stopMeasuring = Seconds(2);

// -------------------------------------------------- //
// --- END OF SIMULATION CONFIGURATION PARAMETERS --- //
// -------------------------------------------------- //

void updateStopTime() {
    stopTimeTCP = Seconds(stopTimeTCPSeconds);
    stopTimeSimulation = Seconds(stopTimeTCPSeconds + 1);
}

void InstallBulkSend(Ptr<Node> node, Ipv4Address address, uint16_t port, std::string socketFactory, Time startTime){
    BulkSendHelper source(socketFactory, InetSocketAddress(address, port));
    ApplicationContainer sourceApps = source.Install(node);
    // // Add random start time between 0 to 2 seconds
    Time randomStartTime = Seconds(Random::get<double>(0.f, 2.f)) + startTime ;
    // std::cout << "src start time " << ((double)randomStartTime.GetMilliSeconds()/1000) << std::endl;

    sourceApps.Start(startTime);
    sourceApps.Stop(stopTimeTCP);
}

void InstallPacketSink(Ptr<Node> dest, uint16_t port, std::string socketFactory, Time startTime, Time stopTime){
    PacketSinkHelper sink (socketFactory, InetSocketAddress(Ipv4Address::GetAny(), port));
    ApplicationContainer sinkApps = sink.Install(dest);
    StaticCast<PacketSink>(sinkApps.Get(0));
    sinkApps.Start(startTime);
    sinkApps.Stop(stopTime);
}

void EnableLinkLoss(NetDeviceContainer ndToEnableLinkLoss, double percentage){
    // Enable link loss on initial direction of link
    Ptr<PointToPointNetDevice> device = ndToEnableLinkLoss.Get(0)->GetObject<PointToPointNetDevice> ();
    device->EnableLinkLoss(percentage, true, tracesPath);

    // // Enable link loss for other direction as well
    // device = ndToEnableLinkLoss.Get(1)->GetObject<PointToPointNetDevice> ();
    // device->EnableLinkLoss(percentage, true, tracesPath);
}

void UpdateSyncTracker(FlowSeqId id, int flow, Time syncExpire, Time measureExpire, bool active){
    syncTrackerPtr->id = id;
    syncTrackerPtr->flow = flow;
    syncTrackerPtr->syncExpire = syncExpire;
    syncTrackerPtr->measureExpire = measureExpire;
    syncTrackerPtr->active = active;
}

FlowSeqId GetFlowSeqId(Ptr<const Packet> pkt) {
    Ipv4Header * ipHeader = (Ipv4Header *) pkt->extractIpHeader();
    TcpHeader * tcpHeader = (TcpHeader *) pkt->extractTcpHeader();
    uint32_t seq = tcpHeader->GetSequenceNumber().GetValue();
    uint32_t flowId;
    uint32_t packetId;

    Ptr<Packet> pkt2 = Copy(pkt);

    if (!classifier.CustomClassify(*ipHeader, *tcpHeader, &flowId, &packetId) ){
        std::cout << "WARNING: CLASSIFYING FAIL" << std::endl;
    }
    return FlowSeqId{flowId, seq};
}

void UpdateTargetZe(){
    if (syncedPkts > 0 && rtxFlow1 > 0 && rtxFlow2 > 0 && rtxFlow1 != syncedPkts && rtxFlow2 != syncedPkts){
        double phat1 = (double) rtxFlow1 / (double) syncedPkts;
        double phat2 = (double) rtxFlow2 / (double) syncedPkts;

        // std::cout << "syncedPkts " << syncedPkts << " rtx flow 1 " << rtxFlow1 << " 2 " << rtxFlow2 << std::endl;
        // std::cout << "phat 1 " << phat1 << " phat 2 " << phat2 << std::endl;
        double currentZE1 = std::sqrt(syncedPkts / (phat1 * (1 - phat1)));
        double currentZE2 = std::sqrt(syncedPkts / (phat2 * (1 - phat2)));

        // std::cout << "currentZE1 " << currentZE1 << " 2 " << currentZE2 << std::endl;
        targetZeAchieved = (bool) (currentZE1 > targetZE && currentZE2 > targetZE) ;
    }
}

void SyncGeneral(Ptr< const Packet > pkt, int flow) {
    if (!timeBound ) {
        if (targetZeAchieved) {
            if (targetSeconds == 0) {
                std::cout << "Target ze achieved, syncedPkts : " << syncedPkts << std::endl;
                targetSeconds = Simulator::Now().GetMilliSeconds() / 1000.0 ;
            }
            return;
        }
        UpdateTargetZe();
    }
    Time now = Simulator::Now();
    Time syncExpire = now + t_delta;
    Time measureExpire = now + stopMeasuring;

    std::map<FlowSeqId, Time> *lookoutRTXPtr;
    std::map<FlowSeqId, Time> *otherlookoutRTXPtr;
    uint32_t *rtxCounterPtr;

    FlowSeqId id = GetFlowSeqId(pkt);

    if (flow == 1) {
        lookoutRTXPtr = &lookoutRTX1;
        otherlookoutRTXPtr = &lookoutRTX2;
        rtxCounterPtr = &rtxFlow1;
    }
    else {
        lookoutRTXPtr = &lookoutRTX2;
        otherlookoutRTXPtr = &lookoutRTX1;
        rtxCounterPtr = &rtxFlow2;
    }

    // check for seq in syncmap and time < expire
    if ((*lookoutRTXPtr).find(id) != (*lookoutRTXPtr).end() && now <= (*lookoutRTXPtr)[id]){
        (*rtxCounterPtr) ++;
        (*lookoutRTXPtr).erase(id);
    }

    // all conditions meet to synchronise
    if (syncTrackerPtr->active && syncTrackerPtr->flow != flow && syncTrackerPtr->syncExpire >= now) {
        syncedPkts ++;
        
        (*lookoutRTXPtr)[id] = measureExpire;
        (*otherlookoutRTXPtr)[syncTrackerPtr->id] = syncTrackerPtr->measureExpire;

        syncTrackerPtr->active = false;
    }
    else{
        UpdateSyncTracker(id, flow, syncExpire, measureExpire, true);
    }
}

void CleanMaps() {
    Time now = Simulator::Now();
    auto it = lookoutRTX1.begin();
    while (it != lookoutRTX1.end()) {
        if (it->second > now) {
            it = lookoutRTX1.erase(it);
        } else {
            ++it;
        }
    }

    it = lookoutRTX2.begin();
    while (it != lookoutRTX2.end()) {
        if (it->second > now) {
            it = lookoutRTX2.erase(it);
        } else {
            ++it;
        }
    }
}

void SyncPacketFlow1(Ptr< const Packet > pkt) {
    
    SyncGeneral(pkt, 1);
    totalFlow1 ++;
}

void SyncPacketFlow2(Ptr< const Packet > pkt) {
    SyncGeneral(pkt, 2);
    totalFlow2 ++;
}


int main (int argc, char *argv[]){
    // Command line arguments
    CommandLine cmd;
    cmd.AddValue("seed", "The random seed", seed);
    cmd.AddValue("stopTime", "The stop time for TCP", stopTimeTCPSeconds);
    cmd.AddValue("linkLoss1", "link1 packet loss rate", linkLoss1);
    cmd.AddValue("linkLoss2", "link2 packet loss rate", linkLoss2);
    cmd.AddValue("delta", "Maximum time between packets we can still sync", delta);
    cmd.AddValue("throughput", "throughtput for all other links - at least srcThroughput 1 + 2", throughput);
    cmd.AddValue("srcThroughput1", "throughput at src 1", srcThroughput1);
    cmd.AddValue("srcThroughput2", "throughput at src 2", srcThroughput2);
    cmd.AddValue("tcpFlows", "Number of tcp flows on each path", tcpFlows);
    cmd.AddValue("lossDuration", "Duration for total loss starting at t=10s", lossDuration);
    cmd.AddValue("targetZE", "Turn on time-unbound mode", targetZE);
    cmd.Parse(argc, argv);

    if (targetZE != 0){
        timeBound = false;
    }
    t_delta = Seconds(delta);
    updateStopTime();

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

    // Create the point-to-point link helpers and connect two router nodes
    PointToPointHelper pointToPointRouter;
    pointToPointRouter.SetDeviceAttribute("DataRate", StringValue(throughput));
    pointToPointRouter.SetChannelAttribute("Delay", StringValue(delay));

    PointToPointHelper leftToRouter1;
    leftToRouter1.SetDeviceAttribute("DataRate", StringValue(srcThroughput1));
    leftToRouter1.SetChannelAttribute("Delay", StringValue(delay));
    PointToPointHelper leftToRouter2;
    leftToRouter2.SetDeviceAttribute("DataRate", StringValue(srcThroughput2));
    leftToRouter2.SetChannelAttribute("Delay", StringValue(delay));

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
    EnableLinkLoss(r3r5ND, linkLoss1);
    EnableLinkLoss(r4r5ND, linkLoss2);

    NetDeviceContainer src1r1ND = leftToRouter1.Install(leftNodes.Get(0), routers.Get(0));
    Names::Add("IntCl1->R1", src1r1ND.Get(0));
    Names::Add("IntR1->Cl1", src1r1ND.Get(1));
    NetDeviceContainer src2r1ND = leftToRouter2.Install(leftNodes.Get(1), routers.Get(0));
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

    std::cout << "From src 1 : " << src1r1IPAddress.GetAddress(0) << " > " << r5destIPAddress.GetAddress(1) << std::endl;
    std::cout << "From src 2 : " << src2r1IPAddress.GetAddress(0) << " > " << r5destIPAddress.GetAddress(1) << std::endl;

    Ipv4GlobalRoutingHelper::PopulateRoutingTables();

    uint16_t server_port = 50000;

    for (int i = 0; i < tcpFlows; i++) {
        /* Install packet sink at receiver side */  
        InstallPacketSink(rightNodes.Get(0), server_port + i, socketFactory, Seconds(0.01), stopTimeSimulation);
        
        /* Install Send application */
        InstallBulkSend(leftNodes.Get(0), r5destIPAddress.GetAddress(1), server_port + i, socketFactory, Seconds(0.1));
        InstallBulkSend(leftNodes.Get(1), r5destIPAddress.GetAddress(1), server_port + i, socketFactory, Seconds(0.1));
    }

    // Install at router 2
    r2r3ND.Get(0)->TraceConnectWithoutContext("MacTx",  MakeCallback(&SyncPacketFlow1));
    r2r4ND.Get(0)->TraceConnectWithoutContext("MacTx",  MakeCallback(&SyncPacketFlow2));

    pointToPointRouter.EnableAsciiAll(tracesPath);

    // schedule link loss
    Simulator::Schedule(Seconds(10), &EnableLinkLoss, r5destND, 1);
    Simulator::Schedule(Seconds(10 + lossDuration), &EnableLinkLoss, r5destND, 0);

    // periodically clean maps
    for (int i = 5; i < stopTimeTCPSeconds; i += 5) {
        Simulator::Schedule(Seconds(i), &CleanMaps);
    }

    Simulator::Stop(stopTimeSimulation);
    Simulator::Run();
    Simulator::Destroy();
    delete syncTrackerPtr;
    if (!timeBound && !targetZeAchieved) {
        std::cout << "Not enough time for time-unbound mode" << std::endl;
        double phat1 = (double) rtxFlow1 / (double) syncedPkts;
        double phat2 = (double) rtxFlow2 / (double) syncedPkts;

        double currentZE1 = std::sqrt(syncedPkts / (phat1 * (1 - phat1)));
        double currentZE2 = std::sqrt(syncedPkts / (phat2 * (1 - phat2)));

        std::cout << "currentZE1 " << currentZE1 << " 2 " << currentZE2 << std::endl;
        return 1;
    }

    std::cout << "Finished simulation!" << std::endl;
    if (!timeBound && targetZeAchieved) std::cout << "Time [" << targetSeconds << "]" << std::endl;
    std::cout << "Ran for [" << stopTimeTCP.GetSeconds() << "] seconds, throughput ["<< throughput << "], srcThroughput1 [" << srcThroughput1 << "], srcThroughput2 [" << srcThroughput2 << "], num of flows [" << tcpFlows << "] linkLoss1 [" << linkLoss1 << "], linkLoss2 [" << linkLoss2 << "]." << std::endl;

    std:: cout << "Synced ["<< syncedPkts << "] total packets (n) , and ["<< ((double) syncedPkts) / stopTimeTCP.GetSeconds() << "] packets per second" << std::endl;
    std::cout << "Total size of rtxs 1 [" << rtxFlow1 << "] 2 [" << rtxFlow2 << "]" << std::endl;
    std::cout << "Total number of sent pkts 1 [" << totalFlow1 << "] 2 [" << totalFlow2 << "]" << std::endl;
    std::cout << "Lambda 1 [" << ((double) totalFlow1 / stopTimeTCPSeconds) << "] 2 [" << ((double) totalFlow2 / stopTimeTCPSeconds)  << "]" << std::endl;
    
    return 0;
}
