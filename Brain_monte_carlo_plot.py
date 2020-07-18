import matplotlib.pyplot as plt
import copy as cp
import collections as clct
import ast

with open('avlch_stats.txt', 'r') as stats_file:
    raw_stats = ast.literal_eval(stats_file.read())


washed_stats = {}
for key, value in raw_stats.items() :
    if key <= 200 :
        washed_stats[key] = cp.deepcopy(raw_stats[key])

ordered_washed_stats = clct.OrderedDict(sorted(raw_stats.items()))


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlabel('Avalanche Scale')
ax1.set_ylabel('Occurence observed')
ax1.set_xscale('log')
ax1.set_yscale('log')

for key, value in ordered_washed_stats.items() :
    point = ax1.scatter(key, value)

plt.show()
