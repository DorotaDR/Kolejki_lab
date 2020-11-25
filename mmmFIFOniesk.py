import numpy as np
import queue
import copy
import matplotlib.pyplot as plt
from itertools import combinations

from scipy.special import binom
from math import factorial

m = 3  # kanały równoległę - lądowanie przeglad itd w jedno
#
# pasy przeznaczone do startu
# pasy przeznaczone do lądowania
# stanowiska przeglądu
# stanowiska tankowania
# stanowiska serwisu


# ile samolotow na godz
lambdas = []
lambd = 2
# intensywność
mi_ = [0.2, 0.4, 0.6]
ro_ = [lambd / mi for mi in mi_]


def calculate_SK(k, ro_list):
    sum_sk = 0

    # k elementowe kobinacje bez powtorzen indeksow
    comb = combinations(range(len(ro_list)), k)

    # iloczyn
    for pair in list(comb):
        prod = 1
        for elem_index in pair:
            prod = prod * ro_list[elem_index]
        # dodaj iloczyn do sumy
        sum_sk += prod

    return sum_sk

    # prob_0 = ???


def probability_k_in_sys(k, ro_list):
    m = len(ro_list)
    if k == 0:
        val_to_sum = [calculate_SK(k, ro_list) / (factorial(k) * binom(m, k)) for k in range(m - 1)]
        sum_ = sum(val_to_sum)
        fraction = calculate_SK(m, ro_list) * calculate_SK(m - 1, ro_list) / (
                factorial(m) * (calculate_SK(m - 1, ro_list) - calculate_SK(m, ro_list)))
        prob_k = 1 / (1 + sum_ + fraction)
    elif k > 0 and k < m:
        prob_k = probability_k_in_sys(0, ro_list) * calculate_SK(k, ro_list) / (factorial(k) * binom(m, k) ** (k - m))
    else:
        prob_k = probability_k_in_sys(0, ro_list) * calculate_SK(k, ro_list) / (
                factorial(k) * calculate_SK(m - 1, ro_list) ** (k - m))
    return prob_k


def calculate_K(ro_list):
    # Średnia liczba zgłoszeń w systemie:
    m = len(ro_list)
    tmp_numerator = probability_k_in_sys(0, ro_list) * calculate_SK(m - 1, ro_list)
    tmp_denominator = factorial(m) * (calculate_SK(m - 1, ro_list) / calculate_SK(m, ro_list) - 1) ** 2

    tmp_numerator_2 = m * probability_k_in_sys(m, ro_list)
    tmp_denominator_2 = calculate_SK(m - 1, ro_list)

    K = tmp_numerator / tmp_denominator + tmp_numerator_2 / tmp_denominator_2

    return K


prob = probability_k_in_sys(5, ro_)

K = calculate_K(ro_)
print(prob)
print(K)



# kod MM1 FIFO niesk
# # Input Parameters
# total_time = int(input("Enter time for simulation (Hours): "))
# IAT_rate = int(input("Enter Job Arrival Rate (/Hour): "))
# ST_rate = int(input("Enter Job Service Rate (/Hour): "))
# rho = IAT_rate / ST_rate
#
# # Initialize Parameters
# qu = queue.Queue()
# curr_process = None
# IAT = []
# ST = []
# AT = []
# wait_time = []
# server_busy = False
# list_wait = []
# list_delay = []
#
# num_processes = int(np.random.poisson(IAT_rate) * total_time)
# num_processes_served = 0
#
# # Populate Inter-Arrival-Times (IAT)
# for i in range(num_processes):
#     temp = np.random.exponential(1 / IAT_rate) * 60 * 60
#     if i == 0:
#         IAT.append(0)
#     else:
#         IAT.append(int(temp - temp % 1))
#
# # Populate Service-Times (ST) (where ST[i]!=0)
# while not len(ST) == num_processes:
#     temp = np.random.exponential(1 / ST_rate) * 60 * 60
#     if not int(temp - temp % 1) < 1:
#         ST.append(int(temp - temp % 1))
#
# # Save a copy of ST
# ST_copy = copy.deepcopy(ST)
#
# # Get Arrival-Times (AT) from IAT starting at t=0
# # and initialize Waiting-Times to 0
# for i in range(num_processes):
#     if i == 0:
#         AT.append(0)
#     else:
#         AT.append(AT[i - 1] + IAT[i])
#     wait_time.append(0)
#
# # Simulation of M/M/1 Queue (i represents current time)
#
# for i in range(total_time * 60 * 60):
#
#     if server_busy:
#         for item in list(qu.queue):
#             wait_time[item] = wait_time[item] + 1
#         ST[curr_process] = ST[curr_process] - 1
#         if ST[curr_process] == 0:
#             server_busy = False
#             num_processes_served = num_processes_served + 1
#
#     for j in range(num_processes):
#         if i == AT[j]:
#             qu.put(j)
#
#     if not server_busy and not qu.empty():
#         curr_process = qu.get()
#         server_busy = True
#
#     sum_wait = 0
#     sum_delay = 0
#
#     for i in range(m):
#
#         for i in range(num_processes_served):
#             sum_wait = sum_wait + wait_time[i]
#             sum_delay = sum_delay + wait_time[i] + ST_copy[i]
#
#         if num_processes_served == 0:
#             list_wait.append(0)
#             list_delay.append(0)
#         else:
#             list_wait.append(sum_wait / (num_processes_served * 60 * 60))
#             list_delay.append(sum_delay / (num_processes_served * 60 * 60))
#
# plt.plot([i + 1 for i in range(total_time * 60 * 60)], list_wait)
# plt.ylabel("Avg Wait Times")
# plt.show()
#
# plt.plot([i + 1 for i in range(total_time * 60 * 60)], list_delay)
# plt.ylabel("Avg Delay Times")
# plt.show()
