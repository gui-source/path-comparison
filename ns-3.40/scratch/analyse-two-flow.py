# %%
import math
import random
import re
import shutil
import os
from enum import Enum

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm

import scipy.stats as st

# %%
# INITALISE DICT
REAL_COUNTS = {}

# %%
# configs
Lambda = 5000
seconds = 100
delta = 0.001

# cwd
os.chdir(r'/Users/guichengtong/Desktop/ns-allinone-3.40/ns-3.40')

print(os.getcwd())

# %%

class EventType(Enum):
    ENQUEUE = 1
    SEND_OR_DEQUEUE = 2
    RECEIVE = 3

def run_simple(Lambda=1000):
    os.system("rm -rf traces")
    os.system("mkdir traces")
    os.system("./ns3 run scratch/two-flow.cc -- -lambda_1={} > traces/two-flow-logs.txt 2>&1".format(Lambda))
    print("run for \lambda={}".format(Lambda))

def run_random(l, save_data = True):
    os.system("rm -rf traces")
    os.system("mkdir traces")
    os.system("./ns3 run scratch/two-flow.cc -- -lambda_1={} -seed={} > traces/two-flow-logs.txt 2>&1".format(l, random.randrange(1,2**32)))
    print("run for \lambda={}, seconds={}".format(l, seconds))
    # save data
    if save_data:
        cwd = os.getcwd()
        os.system("rm -rf {}/data/l_{}-s_{}/traces".format(cwd, l, seconds))
        shutil.copytree("{}/traces".format(cwd), "{}/data/l_{}-s_{}/traces".format(cwd, l, seconds))

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
    id_matcher = re.search("id \d*", line)
    id_no = int(line[id_matcher.span()[0]+3:id_matcher.span()[1]])
    
    return {"link": link,
            "event": event, 
            "time": time, 
            "len":length,
            "src_dest":src_dest_str,
            "id":id_no}

def parse_link(host_node, int1, int2, events, Lambda, use_saved_data=True):
    # returns the link string
    filename = "traces/-{}-Int{}->{}.tr".format(host_node, int1, int2)

    if use_saved_data:
        filename = "data/l_{}-s_{}/".format(Lambda, seconds) + filename
        # print("using saved data from", filename)

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


# %%
# remake all data in range of Lambda 
for i in range(1, 11):
    Lambda = i * 1000
    run_random(Lambda)

# %%
# parse link on Se1, plot inter-arrival times between packets
def plot_inter_arrival_times(Lambda):
    # only use saved data
    events = []
    parse_link("Router1", "R1", "Se1", events, Lambda)

    s_e =  sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, events), key=lambda x: x["time"])

    first_p_time =  s_e[0]["time"]
    prev_time = first_p_time
    last_p_time = s_e[-1]["time"]

    print(first_p_time, last_p_time)

    times = []

    for i in range(1, len(s_e)):
        times.append(s_e[i]["time"] - prev_time)
        prev_time = s_e[i]["time"]
    
    
    plt.figure()
    plt.title("Lamba = {}".format(Lambda))
    plt.xlabel("Inter-arrival times")
    plt.ylabel("Frequency")

    plt.hist(times, bins=50)

    plt.show()

plot_inter_arrival_times(10000)

# %%
# proves that our process is poisson - graph frequency of packets per second
def plot_distribution():
    events = []
    parse_link("Router1", "R1", "Se1", events, Lambda)

    s_e =  sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, events), key=lambda x: x["time"])

    first_p_time =  s_e[0]["time"]
    last_p_time = s_e[-1]["time"]

    print(first_p_time, last_p_time)

    # number of bins = number of seconds
    num_bins = math.ceil(last_p_time - first_p_time) 
    print("Number of bins: {}".format(num_bins))
    bins = range(num_bins)
    counts = [0 for _ in bins]

    # need to adjust times to be time - first_p_time

    for e in s_e:
        # count number of packets in each second
        # print(math.floor(e["time"] - first_p_time))
        counts[math.floor(e["time"] - first_p_time)] += 1

    print("total of {} packets".format(sum(counts)))
    
    plt.figure()
    plt.title("Frequency Distribution for $\lambda$={}".format(Lambda))
    plt.xlabel("Number of packets per second")
    plt.ylabel("Frequency")
    # ignore the last second
    plt.hist(counts[:-1], bins=20)

    plt.show()

Lambda=1000
plot_distribution()

# %%
# helper functions to count synchronised packets

class DefaultStrategy:
    # memory of 1, only sync if packets are within delta
    def __init__(self, sorted_events, link1, link2):
        self.events = sorted_events
        self.link1 = link1
        self.link2 = link2

    def count_sync(self, delta):
        sync_count = 0
        prev_e = None
        for e in self.events:
            if prev_e is None:
                prev_e = e
            
            else:
                # sync condition - e is of different stream within prev_e[time] + delta
                if e["link"] == prev_e["link"] and e["time"] < prev_e["time"] + delta:
                    sync_count += 1
                    prev_e = None
                # otherwise update prev event
                else:
                    prev_e = e
        
        return sync_count
    
def count_sync(Lambda, delta, seconds, strategy=DefaultStrategy, re_init_data =False):
    # return number of packets that can be synced per second
    # re_init_data - option to remake data in dictionary
    global REAL_COUNTS
    save_key = "{}-{}".format(Lambda, delta)

    # if we do not want to remake data and we do have data, we just return 
    if (not re_init_data) and save_key in REAL_COUNTS.keys():
        return REAL_COUNTS[save_key] 
    
    print("MAKING DATA FOR {}".format(save_key))

    s1_events = []
    l1 = parse_link("Router1", "R1", "Se1", s1_events, Lambda)
    s2_events = []
    l2 = parse_link("Router1", "R1", "Se2", s2_events, Lambda)

    events =  sorted(filter(lambda x: x["event"] == EventType.SEND_OR_DEQUEUE, s1_events + s2_events), key=lambda x: x["time"])

    counter = strategy(events, l1, l2).count_sync(delta)

    REAL_COUNTS[save_key] = counter / seconds

    return counter / seconds

count_sync(3000, delta, seconds, re_init_data=True)

# %%
# helper functions for plotting 

def plot_delta_dist(Lambda, seconds, remake_data=False):
    # plot the number of packets synced for different values of delta
    deltas = [i/10000 for i in range(1, 10)]
    counts = [count_sync(Lambda, d, seconds, re_init_data=remake_data) for d in deltas]

    plt.plot(deltas, counts, label="$\lambda$ = {}".format(Lambda))

def get_expected(Lambda, Delta):
    if (1/Delta) < Lambda:
        return math.floor(Lambda / 2)
    else:
        return math.floor((1/Delta) * ((1- math.e ** (-Delta*Lambda)) ** 2)) 

def plot_expected_sync(Lambda):
    # plot the number of packets synced for different values of delta
    deltas = [i/100000 for i in range(1, 100)]
    expecteds = [get_expected(Lambda, d) for d in deltas]

    plt.plot(deltas, expecteds, label="$\lambda$ = {}".format(Lambda))

    # to show the point of delta = 1/Lambda
    plt.scatter([1/Lambda], [get_expected(Lambda, 1/Lambda)], c='black', marker='o')

def get_confidence(Lambda, delta, seconds, E=0.001, l1 = True):
    # returns total packets on both paths
    n = count_sync(Lambda, delta, seconds) * seconds 
    p_hat = 0.1

    z = math.sqrt(((E**2) * n) / (p_hat * (1 - p_hat)))

    # since we are two-tailed, we want to remove the -z score
    confidence = st.norm.cdf(z) - st.norm.cdf(-z)

    return confidence

def plot_confidence_against_delta(Lambda, seconds):
    deltas = [i/10000 for i in range(1, 10)]
    confidence = [get_confidence(Lambda, d, seconds) for d in deltas]

    plt.plot(deltas, confidence, label="$\lambda$ = {}".format(Lambda))

# %%
# plot number of packets synced against delta
    # first run takes ~30 min to make new data

fig = plt.figure()
ax = plt.subplot(111)
plt.title("Number of packets synced against delta")

for i in range(1, 11):
    plot_delta_dist(i * 1000, seconds)

plt.xlabel("Delta, $s$")
plt.ylabel("Number of packets synced per second")
ax.legend(bbox_to_anchor=(1.05, 1))
plt.show()

# %%
# plot expected number of packets against delta

fig = plt.figure()
ax = plt.subplot(111)
plt.title("Expected number of packets synced against delta")
plt.ylim(top=7400, bottom=-200)

for i in range(1, 11):
    plot_expected_sync(i * 1000)

plt.xlabel("Delta, $s$")
plt.ylabel("Number of packets synced per second")
plt.legend()
ax.legend(bbox_to_anchor=(1.4, 1))
plt.show()

# %%
# plot confidence against delta

fig = plt.figure()
ax = plt.subplot(111)
plt.title("Lower bound on confidence for simulation")
for i in range(1, 11):
    plot_confidence_against_delta(i * 1000, seconds)
plt.xlabel("Delta, $s$")
plt.ylabel("Confidence")
plt.legend()
ax.legend(bbox_to_anchor=(1.4, 1))
plt.show()

# %%

def get_optimal_delta(Lambda, seconds, max_diff=100):
    deltas = [i/10000 for i in range(1, 10)]
    counts = [count_sync(Lambda, d, seconds) for d in deltas]
    count_diffs = [counts[i+1] - counts[i]  for i in range(len(counts) - 1)]
    # if count_diffs < 100, take as optimal
    opt_delta = [i for i in range(len(count_diffs)) if count_diffs[i] < max_diff][0]
    
    print("optimal delta value for lambda {} is {}".format(Lambda,deltas[opt_delta]))
    return deltas[opt_delta]

def plot_optimal_deltas(seconds):
    # plot the optimal delta value for each value of Lambda
    lambdas = np.array(range(1000, 10001, 1000))

    plt.figure()
    plt.title("Optimal delta against lambda for different thresholds")

    for t in range(50, 100, 10):
            
        opt_deltas = []
        for l in lambdas:
            opt_deltas.append(get_optimal_delta(l, seconds, t))

        # plot line of best fit
        a, b = np.polyfit(lambdas, np.array(opt_deltas), 1)
        plt.plot(lambdas, a * lambdas + b, label="t={}".format(t))
        # plt.scatter(lambdas, opt_deltas)


    plt.ylabel("Delta")
    plt.xlabel("Lambda")
    plt.legend()
    plt.show()

plot_optimal_deltas(seconds)


