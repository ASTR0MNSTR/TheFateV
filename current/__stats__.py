from scipy.optimize import curve_fit
import random
from scipy.stats import bootstrap
from scipy.stats import spearmanr
import numpy as np

def linear_function(x, a, b):
    return a*x + b

def median_position(xs, ys, par, x_bids):
    y_values = empty(x_bids)
    x_values = empty(x_bids)
    for j, item in enumerate(par):
        for i, pair in enumerate(x_bids):
            if item >= pair[0] and item < pair[1]:
                y_values[i].append(ys[j])
                x_values[i].append(xs[j])
    
    X = [np.median(np.array(item)) if len(item) > 1 else -99 for item in x_values]    
    Y = [np.median(np.array(item)) if len(item) > 1 else -99 for item in y_values]    
    ages = [np.mean(np.array(item)) for item in x_bids]
    
    return X, Y, ages

def empty(bins):
    listed = []
    for item in bins:
        listed.append([])
    return listed

def width_estimation(x, y_mid):
    x_data = np.array(x)
    y_data = np.array(y_mid)
    ssfr = y_data - x_data
    popt, pcov = curve_fit(linear_function, x_data, y_data)
    
    delta_SFR = np.abs(y_data - linear_function(x_data, *popt))
    return popt, np.std(delta_SFR), np.median(ssfr), np.median(ssfr) - np.percentile(ssfr, 16), np.percentile(ssfr, 84) - np.median(ssfr)

def monte_carlo(x, y_mid, y_up, y_down, x_bids):

        y_values = empty(x_bids)
        stmeaner = []
        stmean = []
        bidding = [(pair[0] + pair[1])/2 for pair in x_bids]
        # res = spearmanr(x, y_mid)
        
        for j, item in enumerate(x):
            for i, pair in enumerate(x_bids):
                if item >= pair[0] and item < pair[1]:
                    y_values[i].append([y_mid[j] - y_down[j], y_up[j] - y_mid[j], y_mid[j]])
                    break
        
        length = [len(item) for item in y_values]

        for mean in y_values:
            if len(mean) <= 10:
                stmean.append(-99)
                stmeaner.append(0)
            else:
                medians = []
                for i in range(1001):
                    data = []
                    for pair in mean:
                        if random.random() > 0.5:
                            data.append(pair[2] + abs(random.gauss(0, pair[1])))
                        else:
                            data.append(pair[2] - abs(random.gauss(0, pair[0])))
                    medians.append(np.median(data))

                medians = np.asarray(medians)
                stmean.append(np.median(medians))
                stmeaner.append([[np.median(medians) - np.percentile(medians, 16)], [np.percentile(medians, 84) - np.median(medians)]])
        # return bidding, stmean, stmeaner, length, res
        return bidding, stmean, stmeaner, length
    
def bootstrapper(x, y_mid, y_up, y_down, x_bids):
        y_values = empty(x_bids)
        stmeaner = []
        stmean = []
        x_values = [(pair[0] + pair[1])/2 for pair in x_bids]

        for j, item in enumerate(x):
            for i, pair in enumerate(x_bids):
                if item >= pair[0] and item < pair[1]:
                    y_values[i].append([y_mid[j] - y_down[j], y_up[j] - y_mid[j], y_mid[j]])
                    break
        
        length = [len(item) for item in y_values]

        for mean in y_values:
            if len(mean) <= 3:
                stmean.append(-99)
                stmeaner.append(0)
            else:
                data = [pair[2] for pair in mean]
                values = (data,)
                res = bootstrap(values, np.median, n_resamples=1000)
                medians = res.bootstrap_distribution
                # print(medians)
                stmean.append(np.median(medians))
                medians = np.sort(medians)
                stmeaner.append([[np.median(medians) - medians[160]], [medians[840] - np.median(medians)]])
        return x_values, stmean, stmeaner, length

def empty_out(bins):
    listed = []
    for item in bins:
        listed.append([])
    return listed

def up_lim_analysis(x, ks, x_bids):
    y_values = empty_out(x_bids)
    result = []

    for j, item in enumerate(x):
        for i, pair in enumerate(x_bids):
            if item >= pair[0] and item < pair[1]:
                y_values[i].append(ks[j])
                break
                
    for bin in y_values:
        k = 0
        for element in bin:
            if element == 1:
                k += 1
        result.append(k < 0.5*len(bin))
    return result