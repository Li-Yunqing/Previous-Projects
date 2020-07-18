import random as rdm
import matplotlib.pyplot as plt
import scipy.io as sio
import numpy as np
import copy as cp
import ast
import datetime

start_time = datetime.datetime.now()

nbr_o_steps = 10000
N = 90

r1 = 0.01
r2 = 0.1
Threshold = 0.025

class Node:
    state = 'q'
    r_q_delay = 0


def if_all_quiescent(a_node_list):
    for i in range(len(a_node_list)):
        if a_node_list[i].state == 'e':
            return False
    return True


DTI_data = sio.loadmat('DTI_gunter')['DTI_gunter']
weighted_connection_map = DTI_data * 0.0005

fresh_node_list = []
for i in range(N):
    fresh_node_list.append(Node())


result_storage = []
current_node_list = cp.deepcopy(fresh_node_list)
for step in range(nbr_o_steps):
    new_node_list = cp.deepcopy(fresh_node_list)
    #print(id(current_node_list), id(new_node_list), id(fresh_node_list))
    for i in range(N):
        #print( current_node_list[i].state, new_node_list[i].state, 'before')
        if current_node_list[i].state == 'q' :
            weighted_input = 0
            for j in range(N):
                if current_node_list[j].state == 'e':
                    #print(weighted_connection_map[i,j])
                    weighted_input = weighted_input + weighted_connection_map[i,j]
            indicator = rdm.random()
            #print(r1, indicator)
            if (indicator < r1 or weighted_input > Threshold):
                new_node_list[i].state = 'e'
            else:
                new_node_list[i].state = 'q'
        elif current_node_list[i].state == 'e' :
            new_node_list[i].state = 'r'
        elif current_node_list[i].state == 'r' :
            if (current_node_list[i].r_q_delay == 0 ):
                indicator = rdm.random()
                if (indicator < r2 ):
                    new_node_list[i].state = 'r'
                    new_node_list[i].r_q_delay = rdm.randint(5,40)
                else:
                    new_node_list[i].state = 'q'
                    new_node_list[i].r_q_delay = 0
            elif (current_node_list[i].r_q_delay == 1):
                new_node_list[i].state = 'q'
                new_node_list[i].r_q_delay = 0
            else:
                new_node_list[i].state = 'r'
                new_node_list[i].r_q_delay = current_node_list[i].r_q_delay - 1
        else:
            print('Warning! Error!')
        #print( current_node_list[i].state, new_node_list[i].state, 'after' )
    current_node_list = cp.deepcopy(new_node_list)
    result_storage.append(cp.deepcopy(new_node_list))







with open('avlch_stats.txt', 'r') as stats_file:
    existent_stats = ast.literal_eval(stats_file.read())
#print(type(existent_stats))



avalanche_statistics = cp.deepcopy(existent_stats)
avalanche_size = 0

for step in range(nbr_o_steps):
    if if_all_quiescent(result_storage[step]):
        if (avalanche_size in avalanche_statistics) :
            avalanche_statistics[avalanche_size] = avalanche_statistics[avalanche_size] + 1
        else:
            avalanche_statistics[avalanche_size] = 1
        avalanche_size = 0
    for i in range(N):
        #print(step, i, result_storage[step][i].state, result_storage[step][i].r_q_delay)
        if result_storage[step][i].state == 'e':
            avalanche_size = avalanche_size + 1


#print(avalanche_statistics)

with open('avlch_stats.txt', 'w') as stats_file:
    stats_file.write(str(avalanche_statistics))



'''
state_storage = np.zeros((nbr_o_steps, N))
avarage_state = np.zeros(N)
correlation_accumulated = np.zeros((N, N))

for step in range(nbr_o_steps) :
    for i in range(N) :
        if result_storage[step][i].state == 'e':
            state_storage[step, i] = 1
        else:
            state_storage[step, i] = 0

for i in range(N) :
    avarage_state[i] = np.average(state_storage[:, i])

for step in range(nbr_o_steps) :
    for i in range(N) :
        for j in range(N) :
            correlation_accumulated[i, j] = correlation_accumulated[i, j] + \
            (state_storage[step, i] - avarage_state[i]) * ((state_storage[step, j] - avarage_state[j]))

#print(avarage_state, correlation_accumulated)
    
end_time = datetime.datetime.now()
print((end_time - start_time).seconds)
'''








