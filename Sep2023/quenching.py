from __plt__ import *
from __reader__ import *

import numpy as np
import math

def func(x, A, B):
    return np.log10(A) + (-1*(10**x)/B)*math.log10(math.exp(1))

def func_growth(x, A, B):
    return np.log10(A) + ((10**x)/B)*math.log10(math.exp(1))

source_path = r'E:/LICENSE/ProgsData/main/GAMAforOleg.txt'
input_path = r'E:/backup/backup_BPT/GAMA_ETG_OLA.csv'
output_path = r'E:/databases/Merged.csv'
# merge_phys_databases(source_path, input_path, output_path)
# calculating_mdms(output_path)

DataFrame = pd.read_csv(output_path)

def search_algo(ages_list, age):
    ages = np.full(ages_list.shape, age)
    diff = np.abs(ages - ages_list)
    minima = np.min(diff)
    return np.where(diff == minima)

width = 0.05
ages_list = np.arange(8.85, 9.95, width)
active_list = np.zeros(np.shape(ages_list))
full_list = np.zeros(np.shape(ages_list))

for i in range(len(DataFrame['ager_percentile50'])):
    j = search_algo(ages_list, DataFrame['ager_percentile50'][i])
    full_list[j] += 1
    if DataFrame['WHAN'][i] in ['ELR', 'LLR', 'NER']:
        active_list[j] += 1

k = active_list
n = full_list
error=np.sqrt((k+1.)*(k+2.)/((n+2.)*(n+3.)) - ((k+1.)**2)/((n+2)**2.))
print(error)
Y = active_list/full_list
X = ages_list

popt, pcov = curve_fit(func_growth, X, Y, sigma=error)

A = popt[0]
B = popt[1]
A_err = np.sqrt(np.diag(pcov))[0]
B_err = np.sqrt(np.diag(pcov))[1]
print(f'{A}+-{A_err}; {B/10**9}+-{B_err/10**9}')

plt.bar(X, Y, color='blue', width = width, edgecolor ='grey', align='center')

x = np.arange(8.8, 9.95, 0.01)
plt.plot(x, func_growth(x, A, B), color='black')
plt.ylim(-0.1, 1.1)
plt.show()

# add weights to the data-points
    