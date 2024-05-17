from __plt__ import *
from __reader__ import *

import numpy as np
import math

def func(x, A, B):
    return np.log10(A) + (-1*(10**x)/B)*math.log10(math.exp(1))

def func_growth(x, A, B):
    return np.log10(A) + ((10**x)/B)*math.log10(math.exp(1))

def search_algo(ages_list, age):
    ages = np.full(ages_list.shape, age)
    diff = np.abs(ages - ages_list)
    minima = np.min(diff)
    return np.where(diff == minima)

def plotting_fitting(function, keys, filename, title, ylabel, color, x_text, y_text):
    output_path = r'E:/databases/Merged.csv'
    DataFrame = pd.read_csv(output_path)
    width = 0.10
    ages_list = np.arange(8.85, 9.95, width)
    active_list = np.zeros(np.shape(ages_list))
    full_list = np.zeros(np.shape(ages_list))

    for i in range(len(DataFrame['ager_percentile50'])):
        j = search_algo(ages_list, DataFrame['ager_percentile50'][i])
        full_list[j] += 1
        if DataFrame['WHAN'][i] in keys:
            active_list[j] += 1

    k = active_list
    n = full_list
    error=np.sqrt((k+1.)*(k+2.)/((n+2.)*(n+3.)) - ((k+1.)**2)/((n+2)**2.))
    print(error)
    Y = active_list/full_list
    X = ages_list

    popt, pcov = curve_fit(function, X, Y, sigma=error, absolute_sigma=True)

    A = popt[0]
    B = popt[1]
    A_err = np.sqrt(np.diag(pcov))[0]
    B_err = np.sqrt(np.diag(pcov))[1]
    print(f'{A}+-{A_err}; {B/10**9}+-{B_err/10**9}')

    plt.bar(X, Y, color=color, width = width, edgecolor ='black', align='center')

    x = np.arange(8.8, 9.95, 0.01)
    plt.plot(x, function(x, A, B), color='black')
    plt.errorbar(X, Y, yerr=error, fmt='o', color='black')
    plt.ylim(-0.1, 1.1)
    plt.text(x_text, y_text, r'$\tau$ = {:.3}+-{:.3} Gyr'.format(B/10**9, B_err/10**9))
    plt.title(title)
    plt.ylabel(ylabel)
    plt.savefig(filename)
    plt.show()

if __name__ == '__main__':
    plotting_fitting(func, ['sAGN', 'wAGN', 'SF'], 'active_quenching.pdf', 'AGNs and SFGs quenching', r'$f_{active}$', 'mediumorchid', 9.6, 1.0)
    plotting_fitting(func_growth, ['NER', 'ELR', 'LLR'], 'passive_quenching.pdf', 'RGs evolution', r'$f_{RGs}$', 'peru', 9.0, 1.0)
# add weights to the data-points
    