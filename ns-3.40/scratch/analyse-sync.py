# %%
import math
import pickle
import random
import re
import time
import os
from enum import Enum

import numpy as np
import matplotlib.pyplot as plt

import scipy.stats as st

# cwd
os.chdir(r'/Users/guichengtong/Desktop/ns-allinone-3.40/ns-3.40')

print(os.getcwd())


# %%
# [MUST RUN] helper functions for running, parsing info, etc.
class EventType(Enum):
    ENQUEUE = 1
    SEND_OR_DEQUEUE = 2
    RECEIVE = 3

def run_test(throughput1, throughput2, link_loss1, link_loss2, delta, seconds, tcp_flows, random_loss, targetZE=None):
    start = time.time()
    os.system("rm -rf traces")
    os.system("mkdir traces")
    seed = random.randrange(1, 2**31)
    args = "-srcThroughput1={}Mbps -srcThroughput2={}Mbps -throughput={}Mbps -linkLoss1={} -linkLoss2={} -stopTime={} -delta={} ".format(
        throughput1,
        throughput2, 
        throughput1+throughput2+10, 
        link_loss1,
        link_loss2,
        seconds,
        delta
    )

    args2 = "-seed={} -tcpFlows={} -randomLoss={}".format(
        seed,
        tcp_flows,
        random_loss
    )

    logs_path = "data/t1{}Mbps-t2{}Mbps-t{}Mbps-l1{}-l2{}-time{}-delta{}-flows{}-randLoss{}-seed{}".format(
        throughput1,
        throughput2, 
        throughput1+throughput2+10, 
        link_loss1,
        link_loss2,
        seconds,
        delta,
        tcp_flows,
        random_loss,
        seed
    )

    if targetZE is not None:
        args2 +=  " -targetZE={}".format(targetZE)
        logs_path += "-ze{}".format(targetZE)
    logs_path += ".txt"
    os.system("./ns3 run scratch/sync-impl.cc -- {} {}> {} 2>&1".format(args, args2, logs_path))

    print("args: {} {}".format(args, args2))

    end = time.time()
    print("ran for ", end - start)
    return logs_path


def get_span(matcher, txt):
    return txt[matcher.span()[0]:matcher.span()[1]]

def parse_info(line, link):
    event = None
    if line[0] == '+':
        event = EventType.ENQUEUE
    elif line[0] == '-':
        event = EventType.SEND_OR_DEQUEUE
    elif line[0] == 'r':
        event = EventType.RECEIVE

    time = float(line[1:].split()[0])
    length_matcher = re.search("(size=\d*)", line)
    if length_matcher is None:
        length = 0
    else:
        length = int(line[length_matcher.span()[0]+5:length_matcher.span()[1]])
    src_dest_matcher = re.search("\d*\.\d*\.\d*\.\d*\s[\>\<]\s\d*\.\d*\.\d*\.\d", line)
    src_dest_str = get_span(src_dest_matcher, line)
    seq_matcher = re.search("Seq=\d*", line)
    seq_no = int(line[seq_matcher.span()[0]+4:seq_matcher.span()[1]])
    id_matcher = re.search("id \d*", line)
    id_no = int(line[id_matcher.span()[0]+3:id_matcher.span()[1]])
    port_matcher = re.search("[0-9]{5} > [0-9]{5}", line)
    port = get_span(port_matcher, line)
    
    return {"link": link,
            "event": event, 
            "time": time, 
            "len":length,
            "src_dest":src_dest_str,
            "seq":seq_no,
            "id":id_no,
            "port":port}

def parse_link(host_node, int1, int2, events, save_data=True):
    # returns the link string
    filename = "traces/-{}-Int{}->{}.tr".format(host_node, int1, int2)

    f = open(filename, "r")
    link = "{}->{}".format(str(int1), str(int2))
    while True:
        line = f.readline()
        if len(line) == 0:
            break
        info = parse_info(line, link)
        events.append(info)

    return link

def event_str(event):
    if event == EventType.ENQUEUE:
        return "Enqueue"
    elif event == EventType.SEND_OR_DEQUEUE:
        return "SEND"
    elif event == EventType.RECEIVE:
        return "RCV"
    else:
        return "ERROR"

def read_logfile(logfile_name):
    f = open(logfile_name, "r")
    txt = ""
    while True:
        line = f.readline()
        if len(line) == 0:
            break
        txt += line
    
    return parse_results_sync(txt)

def parse_results_sync(txt):
    n_matcher = re.search("Synced \[(\d*)\]", txt)
    n = int(n_matcher.group(1))

    rtx_matcher = re.search("Total size of rtxs 1 \[(\d*)\] 2 \[(\d*)\]", txt)
    rtx1 = int(rtx_matcher.group(1))
    rtx2 = int(rtx_matcher.group(2))

    total_pkts_matcher = re.search("Total number of sent pkts 1 \[(\d*)\] 2 \[(\d*)\]", txt)
    total_pkts1 = int(total_pkts_matcher.group(1))
    total_pkts2 = int(total_pkts_matcher.group(2))

    lambda_matcher = re.search("Lambda 1 \[(\d*.?\d*)\] 2 \[(\d*.?\d*)\]", txt)
    lambda1 = float(lambda_matcher.group(1))
    lambda2 = float(lambda_matcher.group(2))

    time_matcher = re.search("Time \[(\d*.?\d*)\]", txt)

    phat1 = rtx1 / n
    phat2 = rtx2 / n
    
    dict = {"n":n, 
            "pkts1": total_pkts1,
            "pkts2": total_pkts2,
            "rtx1": rtx1, 
            "rtx2": rtx2, 
            "lambda1" : lambda1,
            "lambda2" : lambda2,
            "phat1": phat1, 
            "phat2": phat2}
    
    if time_matcher is not None:
        print("not none")
        time = float(time_matcher.group(1))
        dict["time"] = time
    return dict

# %%
# parse link on Se1, plot inter-arrival times between packets
def plot_inter_arrival_times():
    # only use saved data
    events = []
    parse_link("Router1", "R1", "R2", events)

    s_e =  sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, events), key=lambda x: x["time"])

    first_p_time =  s_e[0]["time"]
    prev_time = first_p_time
    last_p_time = s_e[-1]["time"]


    times = []

    for i in range(1, len(s_e)):
        times.append(s_e[i]["time"] - prev_time)
        prev_time = s_e[i]["time"]
    
    plt.figure()
    plt.title("Frequency of packet inter-arrival times")
    plt.xlabel("Time ($s$)")
    plt.ylabel("Frequency")
    plt.yscale("log")

    plt.hist(list(filter(lambda x: x < 0.005, times)), bins=100)

    plt.show()

plot_inter_arrival_times()


# %%
# proves that our process is poisson - graph frequency of packets per second
def plot_distribution():
    events = []
    parse_link("Router1", "R1", "R2", events)

    s_e =  sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, events), key=lambda x: x["time"])

    first_p_time =  s_e[0]["time"]
    last_p_time = s_e[-1]["time"]

    print(first_p_time, last_p_time)

    # number of bins = number of seconds
    num_bins = math.ceil(last_p_time - first_p_time)  + 1
    print("Number of bins: {}".format(num_bins))
    counts = [0 for _ in range(num_bins)]

    # need to adjust times to be time - first_p_time

    for e in s_e:
        # count number of packets in each second
        counts[math.floor(e["time"])] += 1

    print(counts)
        
    # counts = list(filter(lambda x : x > 22084 and x < 22086, counts))

    print("total of {} packets".format(sum(counts)))
    
    plt.figure()
    plt.title("Frequency Distribution")
    plt.xlabel("Number of packets per second")
    plt.ylabel("Frequency")
    # plt.hist(counts, bins=50)
    plt.plot(range(num_bins), counts)

    plt.show()

plot_distribution()


# %%
# verify all packets from src1 10.0.5.1 are r2->r3
def verify_f1_r2r4():
    # only use saved data
    events = []
    parse_link("Router2", "R2", "R3", events)

    # look for events with src from src2
    src2_dest_str = "10.0.6.1 > 10.0.7.2"
    s_e =  sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE,filter(lambda x: x["src_dest"] == src2_dest_str, events)), key=lambda x: x["time"])

    # should be 0
    assert(len(s_e) == 0)

# verify all packets from src2 10.0.6.1 are r2->r4
def verify_f2_r2r3():
    # only use saved data
    events = []
    parse_link("Router2", "R2", "R4", events)

    # look for events with src from src1
    src1_dest_str = "10.0.5.1 > 10.0.7.2"
    s_e =  sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE,filter(lambda x: x["src_dest"] == src1_dest_str, events)), key=lambda x: x["time"])

    # should be 0
    assert(len(s_e) == 0)

verify_f1_r2r4()
verify_f2_r2r3()

# %%
# print total number of packets sent on r1 -> r2 and r5 -> se, group by flow 1 src and flow 2 src

def pkt_loss_src1():
    events = []
    parse_link("Client1", "Cl1", "R1", events)

    src_s_e =  sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE,events), key=lambda x: x["time"])

    events = []
    parse_link("Router5", "R5", "Se", events)
    
    src1_dest_str = "10.0.5.1 > 10.0.7.2"
    post_loss_s_e = sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE,filter(lambda x: x["src_dest"] == src1_dest_str, events)), key=lambda x: x["time"])

    print("pkt loss of src 1 :", 1-(len(post_loss_s_e) / len(src_s_e)))


def pkt_loss_src2():
    events = []
    parse_link("Client2", "Cl2", "R1", events)

    src_s_e =  sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE,events), key=lambda x: x["time"])

    events = []
    parse_link("Router5", "R5", "Se", events)
    
    src1_dest_str = "10.0.6.1 > 10.0.7.2"
    post_loss_s_e = sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE,filter(lambda x: x["src_dest"] == src1_dest_str, events)), key=lambda x: x["time"])

    print("pkt loss of src 2 :", 1-(len(post_loss_s_e) / len(src_s_e)))

pkt_loss_src1()
pkt_loss_src2()

# %%
# [MUST RUN] helper functions to count synchronised packets

class DefaultStrategy:
    # memory of 1, only sync if packets are within delta
    def __init__(self, sorted_events):
        self.events = sorted_events

    def count_sync(self, delta):
        sync_count = 0
        prev_e = None
        for e in self.events:
            if prev_e is None:
                prev_e = e
            
            else:
                # sync condition - e is of different stream within prev_e[time] + delta
                if e["src_dest"] != prev_e["src_dest"] and e["time"] <= prev_e["time"] + delta:
                    sync_count += 1
                    # print(sync_count)
                    prev_e = None
                # otherwise update prev event
                else:
                    prev_e = e
        
        return sync_count
    
    def count_rtx(self, delta, measure_expire):
        sync_count = 0
        rtx_count = {}
        prev_e = None

        dest_to_lookout_rtx = {}

        for e in self.events:
            if e["src_dest"] not in dest_to_lookout_rtx.keys():
                dest_to_lookout_rtx[e["src_dest"]] = {}
                rtx_count[e["src_dest"]] = 0
                print("added lookout rtx dict for ", e["src_dest"])

        for e in self.events:
            lookout_rtx = dest_to_lookout_rtx[e["src_dest"]]
            other_lookout_rtx = dest_to_lookout_rtx[[k for k in dest_to_lookout_rtx.keys() if k != e["src_dest"]][0]]
            port_seq = e["port"] +"|" +  str(e["seq"])
            if port_seq in lookout_rtx.keys() and e["time"] <= lookout_rtx[port_seq]:
                    rtx_count[e["src_dest"]] += 1

            if prev_e is None:
                prev_e = e
            
            else:
                # rtx condition 
                if e["src_dest"] != prev_e["src_dest"] and e["time"] <= prev_e["time"] + delta:
                    sync_count += 1
                    lookout_rtx[port_seq] = e["time"] + measure_expire
                    other_lookout_rtx[prev_e["port"] + "|" + str(prev_e["seq"])] = prev_e["time"] + measure_expire
                    prev_e = None
                # otherwise update prev event
                else:
                    prev_e = e
        return rtx_count
    
class NoSyncStrategy:
    def __init__(self, sorted_events):
        self.events = sorted_events
        self.flow_1_src_dest = "10.0.5.1 > 10.0.7.2"
        self.flow_2_src_dest = "10.0.6.1 > 10.0.7.2"


    def count_rtx(self, foo, bar):
        rtx_count = {}

        dest_to_lookout_rtx = {}

        for e in self.events:
            if e["src_dest"] not in dest_to_lookout_rtx.keys():
                dest_to_lookout_rtx[e["src_dest"]] = set()
                rtx_count[e["src_dest"]] = 0
                print("added lookout rtx dict for ", e["src_dest"])

        for e in self.events:
            lookout_rtx = dest_to_lookout_rtx[e["src_dest"]]
            port_seq = e["port"] +"|" +  str(e["seq"])
            if port_seq in lookout_rtx:
                    rtx_count[e["src_dest"]] += 1

            lookout_rtx.add(port_seq)
            
        return rtx_count

        
    
# return total number of packets that can be synced 
def count_sync(delta, strategy=DefaultStrategy):
    r3r2_events = []
    parse_link("Router2", "R2", "R3", r3r2_events)
    # print(len(sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, r3r2_events), key=lambda x: x["time"])))
        
    r4r2_events = []
    parse_link("Router2", "R2", "R4", r4r2_events)
    # print(len(sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, r4r2_events), key=lambda x: x["time"])))

    events =  sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, r3r2_events + r4r2_events), key=lambda x: x["time"])
    counter = strategy(events).count_sync(delta)

    return counter

# return total number of packets sent
def count_total_pkts(flow):
    events = []
    if flow == 1:
        parse_link("Router2", "R2", "R3", events)
    elif flow == 2:
        parse_link("Router2", "R2", "R4", events)

    events = list(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, events))

    return len(events)

# return total number of rtxs
def count_rtxs(delta, measure_expire, strategy=DefaultStrategy):
    r3r2_events = []
    parse_link("Router2", "R2", "R3", r3r2_events)
    # print(len(sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, r3r2_events), key=lambda x: x["time"])))
        
    r4r2_events = []
    parse_link("Router2", "R2", "R4", r4r2_events)
    # print(len(sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, r4r2_events), key=lambda x: x["time"])))

    events =  sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, r3r2_events + r4r2_events), key=lambda x: x["time"])
    counter_dict = strategy(events).count_rtx(delta, measure_expire)

    return counter_dict

def get_expected(l_1, l_2, Delta):
    exp_time = (l_1 + l_2) / (l_1 * l_2)

    if Delta <= exp_time/2:
        return_val = math.floor((1 / Delta) * (1- math.e ** (-Delta * l_1)) * (1 - math.e ** (-Delta * l_2))) 
    else:
        return_val = math.floor(1 / exp_time)

    return return_val

def zscore_to_alpha(z):
    one_tail_alpha = st.norm.cdf(z)
    p = 1 - one_tail_alpha

    return 1 - (2 * p)

def alpha_to_zscore(alpha):
    p = 1 - alpha
    one_tail_alpha = 1 - (0.5 * p)

    return st.norm.ppf(one_tail_alpha)



# %%
# plot n against delta 

lambdas=[]

def run_and_plot_n(link_loss1, link_loss2):
    seconds = 20
    run_test("10Mbps", "20Mbps", link_loss1, link_loss2, 0.0001, seconds)
    print("lambda1 {}, lambda2 {}".format(count_total_pkts(1)/seconds, count_total_pkts(2)/seconds))
    lambdas.append((count_total_pkts(1)/seconds, count_total_pkts(2)/seconds))

    # plot the number of packets synced for different values of delta
    deltas = [i/10000 for i in range(1, 100)]
    counts = [count_sync(d) for d in deltas]
    plt.plot(deltas, counts, label="ll1:{} ll2:{}".format(link_loss1, link_loss2))

plt.figure()
plt.title("Number of packets that can be synced $n$ against $\delta$")
plt.xlabel("Delta, $s$")
plt.ylabel("Total number of packets synced")
run_and_plot_n(0.001, 0.002)
run_and_plot_n(0.002, 0.001)
run_and_plot_n(0.0001, 0.0002)
run_and_plot_n(0.0002, 0.0001)
run_and_plot_n(0, 0)
plt.legend()

# %%
#  helper functions for plotting 

def plot_expected_sync(lambda1, lambda2, seconds):
    # plot the number of packets synced for different values of delta
    deltas = [i/10000 for i in range(1, 100)]
    expecteds = [get_expected(lambda1, lambda2, d) * seconds for d in deltas]

    plt.plot(deltas, expecteds, label="$\lambda$ = {}, {}".format(lambda1, lambda2))

    # # to show the point of delta = 1/Lambda
    # plt.scatter([1/Lambda], [get_expected(Lambda, Lambda, 1/Lambda)], c='black', marker='o')


# %%
# plot E(n) against delta
fig = plt.figure()
ax = plt.subplot(111)
plt.title("Expected number of packets synced against delta")
# same as 
# plt.ylim(top=10200, bottom=-200)

plot_expected_sync(12561.7, 9638.18, 100)

# for (l1,l2) in lambdas:
#     l_min = min(l1, l2)
#     plot_expected_sync(l_min, l_min, 20)


plt.xlabel("Delta, $s$")
plt.ylabel("Total number of packets synced")
plt.legend()
ax.legend(bbox_to_anchor=(1.45, 1))
plt.show()

# %%
# [MUST RUN] define advice functions

def advice_min_sync_needed(alpha, E, max_phat=0.5):
    # returns advice on the minimum number of packets needed to sync to achieve \alpha, E
    z = alpha_to_zscore(alpha)

    n = max_phat * (1 - max_phat) * pow(z, 2) / pow(E, 2)
    
    return n

def advice_run_test_time(alpha, E, l1, l2, d, max_phat=0.5):
    # returns time in seconds for minimum time to achieve \alpha, E given we have l1 = lambda_1, l2 = lambda_2, d
    exp_sync_per_sec = get_expected(l1, l2, d)
    min_sync_total = advice_min_sync_needed(alpha, E, max_phat)

    time_needed = min_sync_total / exp_sync_per_sec

    return time_needed

def advice_alpha(E, n, max_phat = 0.5):
    # returns advice on minimum confidence level for E, n

    min_zscore = math.sqrt(pow(E, 2) * n / (max_phat * (1 - max_phat)))

    return zscore_to_alpha(min_zscore)

def advice_precision(alpha, n):
    # returns maximum precision we can achieve for \alpha, n
    z = alpha_to_zscore(alpha)
    phat = 0.5

    E = z * math.sqrt((phat * (1 - phat) / n))

    return E

print(advice_min_sync_needed(.98, 0.0001, max_phat=0.01))
print(advice_run_test_time(.98, 0.0001, 10000, 10000, 0.001, 0.001))

# %%
# [MUST RUN] - confidence, precision, classic mode, time-unbound mode

def get_confidence(n, E, phat):
    # returns total packets on both paths
    z = math.sqrt(((E**2) * n) / (phat * (1 - phat)))

    # since we are two-tailed, we want to remove the -z score
    confidence = st.norm.cdf(z) - st.norm.cdf(-z)

    return confidence

def change_confidence(old_E, new_E, conf):
    z = alpha_to_zscore(conf)
    z = z * new_E / old_E

    new_conf = st.norm.cdf(z) - st.norm.cdf(-z)

    return new_conf

    # returns maximum precision we can achieve for \alpha, n
def get_precision(n, alpha, phat):
    z = alpha_to_zscore(alpha)

    E = z * math.sqrt((phat * (1 - phat) / n))

    return E

# classic mode - run test for fixed time and E or alpha
# returns [flow1, flow2]
def run_classic_mode (throughput1, throughput2, seconds, delta, ll1, ll2, E=None, alpha=None, random_loss=0, tcpFlows=20):
    logspath = run_test(throughput1, throughput2, ll1, ll2, delta, seconds, tcpFlows, random_loss)

    # flow_1_src_dest = "10.0.5.1 > 10.0.7.2"
    # flow_2_src_dest = "10.0.6.1 > 10.0.7.2"

    return_dict = read_logfile(logspath)

    n = return_dict["n"]
    rtx1 = return_dict["rtx1"]
    rtx2 = return_dict["rtx2"]
    phat1 = rtx1 / n
    phat2 = rtx2 / n


    if (E is not None) and (alpha is not None):
        print("INVALID ARGS: E AND alpha PROVIDED")
        raise Exception
    if (E is not None):
        return_dict["alpha1"] = get_confidence(n, E, phat1)
        return_dict["alpha2"] = get_confidence(n, E, phat2)
    if (alpha is not None):
        return_dict["E1"] = get_precision(n, alpha, phat1)
        return_dict["E2"] = get_precision(n, alpha, phat2)

    return return_dict

def time_unbound_mode (throughput1, throughput2, delta, ll1, ll2, E, alpha, max_iter=10, random_loss=0, tcpFlows=20, seconds=100):
    success = False
    iters = 0
    ze = alpha_to_zscore(alpha) / E
    while success is False or iters < max_iter:
        iters += 1
        logspath = run_test(throughput1, throughput2, ll1, ll2, delta, seconds, tcpFlows, random_loss, targetZE=ze)
        return_dict = None
        
        try:
            return_dict = read_logfile(logspath)
            success = True
            print("SUCCESS!")
        except Exception as e:
            seconds *= 2
            print("except, retrying with ", seconds)

            continue

        n = return_dict["n"]
        rtx1 = return_dict["rtx1"]
        rtx2 = return_dict["rtx2"]
        phat1 = rtx1 / n
        phat2 = rtx2 / n


        return_dict["alpha1"] = get_confidence(n, E, phat1)
        return_dict["alpha2"] = get_confidence(n, E, phat2)

        return_dict["E1"] = get_precision(n, alpha, phat1)
        return_dict["E2"] = get_precision(n, alpha, phat2)
    return return_dict

# %%
# experiment for sync vs no sync - vary throughput, 

throughputs = [i for i in range(10, 101, 10)]

acc1_plot = []
acc2_plot = []

nosync_acc1_plot = []
nosync_acc2_plot = []

# flow 1 - flow 2
sync_diff_plot=[]
nosync_diff_plot=[]

for t in throughputs:
    ll1 = 0.001
    ll2 = 0.002
        
    acc1s = []
    acc2s = []
    nosync_acc1s = []
    nosync_acc2s = []

    sync_diffs = []
    nosync_diffs = []

    for i in range(5):
        return_dict = run_classic_mode(10, t, 20, 0.01, ll1, ll2, E=0.0001, random_loss=1)

        acc1s.append(abs(return_dict["phat1"] - ll1) / ll1)
        acc2s.append(abs(return_dict["phat2"] - ll2) / ll2)

        nosync_rtxs = count_rtxs(None, None, strategy=NoSyncStrategy)
        nosync_acc1s.append(abs(nosync_rtxs["10.0.5.1 > 10.0.7.2"]/return_dict["pkts1"] - ll1) / ll1)
        nosync_acc2s.append(abs(nosync_rtxs["10.0.6.1 > 10.0.7.2"]/return_dict["pkts2"] - ll2) / ll2)

        sync_diffs.append(return_dict["phat1"] - return_dict["phat2"])
        nosync_diffs.append(nosync_rtxs["10.0.5.1 > 10.0.7.2"]/return_dict["pkts1"] - nosync_rtxs["10.0.6.1 > 10.0.7.2"]/return_dict["pkts2"])

    acc1_plot.append(np.average(acc1s))
    acc2_plot.append(np.average(acc2s))

    nosync_acc1_plot.append(np.average(nosync_acc1s))
    nosync_acc2_plot.append(np.average(nosync_acc2s))

    sync_diff_plot.append(np.average(sync_diffs))
    nosync_diff_plot.append(np.average(nosync_diffs))

# %%

plt.figure(figsize=(12, 8))
plt.title("Difference against throughputs in flow 1")
plt.plot(throughputs, acc1_plot, label="flow 1 sync")
# plt.plot(ll2s, acc2_plot, label="flow 2 sync")

plt.plot(throughputs, nosync_acc1_plot, label="flow 1 no sync")
# plt.plot(ll2s, nosync_acc2_plot, label="flow 2 no sync")

plt.legend()

plt.figure(figsize=(12, 8))
plt.title("Difference against throughputs in flow 1")
plt.plot(throughputs, acc2_plot, label="flow 1 sync")
# plt.plot(ll2s, acc2_plot, label="flow 2 sync")

plt.plot(throughputs, nosync_acc2_plot, label="flow 1 no sync")
# plt.plot(ll2s, nosync_acc2_plot, label="flow 2 no sync")

plt.legend()

plt.figure(figsize=(12, 8))
plt.title("Link loss diferrence (flow 1 - flow 2)")

plt.plot(throughputs, sync_diff_plot, label="sync")


plt.plot(throughputs, nosync_diff_plot, label="no sync")
plt.plot(throughputs, [ll1 - ll2 for _ in throughputs], label="correct")

plt.legend()


# %%
# [SETUP
con1s = []
con2s = []

phat1s = []
phat2s = []

prec=0.001

# %%
# [EXPERIMENT] - vary traffic pattern and which packets lost, check that results do not vary too much


for _ in range(50):
    results = run_classic_mode(20, 20, 40, 0.01, 0.01, 0.02, E=prec, random_loss=0, tcpFlows=20)
    n = results["n"]
    rtx1 = results["rtx1"]
    rtx2 = results["rtx2"]
    phat1 = rtx1 / n
    phat2 = rtx2 / n

    confidence1 = get_confidence(n, prec, phat1)
    confidence2 = get_confidence(n, prec, phat2)

    con1s.append(confidence1)
    con2s.append(confidence2)
    phat1s.append(phat1)
    phat2s.append(phat2)



# %%
new_con1s = [change_confidence(0.0001, 0.0005, c) for c in con1s]
new_con2s = [change_confidence(0.0001, 0.0005, c) for c in con2s]


print("p_hat 1 is {}, mean {}, std is {}".format(0.01, np.mean(phat1s),np.std(phat1s)))
print("p_hat 2 is {}, mean {}, std is {}".format(0.02, np.mean(phat2s),np.std(phat2s)))
print("con 1 mean {}, std is {}".format(np.mean(new_con1s),np.std(new_con1s)))
print("con 2 mean {}, std is {}".format(np.mean(new_con2s), np.std(new_con2s)))


# %%

plt.figure(figsize=(12, 8))
plt.xlim(left=0.99, right=1)
plt.title("Histogram of confidences for flow 1")
plt.hist(new_con1s)

plt.figure(figsize=(12, 8))
plt.xlim(left=0.98, right=0.99)
plt.title("Histogram of confidences for flow 2")
plt.hist(new_con2s, bins=30)

plt.figure(figsize=(12, 8))
plt.xlim(left=0, right=0.02)
plt.title("Histogram of $\hat p$ for flow 1")
plt.hist(phat1s)

plt.figure(figsize=(12, 8))
plt.xlim(left=0.01, right=0.03)
plt.title("Histogram of $\hat p$ for flow 2")
plt.hist(phat2s)


# %%
# [EXPERIMENT] - check if precision correct
prec = 0.0005

ll1 = 0.01
ll2 = 0.02

runs = 50
correct_runs1 = 0
correct_runs2 = 0
new_con1s = [change_confidence(0.0001, prec, c) for c in con1s]
new_con2s = [change_confidence(0.0001, prec, c) for c in con2s]
for i in range(runs):
    # results = run_classic_mode(20, 20, 40, 0.01, ll1, ll2, E=prec, random_loss=0, tcpFlows=20)
    # n = results["n"]
    # rtx1 = results["rtx1"]
    # rtx2 = results["rtx2"]
    # phat1 = rtx1 / n
    # phat2 = rtx2 / n

    # confidence1 = get_confidence(n, prec, phat1)
    # confidence2 = get_confidence(n, prec, phat2)
    phat1 = phat1s[i]
    phat2 = phat2s[i]

    if (phat1 <= (ll1 + prec)) and (phat1 >= (ll1 - prec)):
        correct_runs1 += 1

    if (phat2 <= (ll2 + prec)) and (phat2 >= (ll2 - prec)):
        correct_runs2 += 1
    
print("flow 1 predicted confidence: {} real confidence: {}".format(np.mean(new_con1s), correct_runs1 / runs))
print("flow 2 predicted confidence: {} real confidence: {}".format(np.mean(new_con2s), correct_runs2 / runs))


# %%
# [EXPERIMENT DATA] - time unbound check time

alpha = 0.95
E = 0.001
ll1 = 0.01
ll2 = 0.02
delta = 0.01
required_times = []

for _ in range(50):
    results = time_unbound_mode(20, 20, delta, ll1, ll2, E, alpha, tcpFlows=20, seconds=40)
    required_times.append(results["time"])

# %%

expected_lambda1 = 3000
expected_lambda2 = 3000
advice_time = advice_run_test_time(alpha, E, expected_lambda1, expected_lambda2, delta, max_phat=0.2)


# %%
print(alpha_to_zscore(0.95) / 0.001)


