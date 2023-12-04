from scipy.optimize import curve_fit
import random
from scipy.stats import bootstrap
import numpy as np

def monte_carlo(x, y_mid, y_up, y_down, x_bids):
        y_values = [[], [], [], [], [], []]
        stmeaner = []
        stmean = []
        medians = [(pair[0] + pair[1])/2 for pair in x_bids]

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
                data = []
                for i in range(100):
                    av = []
                    for pair in mean:
                        if random.random() > 0.5:
                            av.append(pair[2] + abs(random.gauss(0, pair[1])))
                        else:
                            av.append(pair[2] - abs(random.gauss(0, pair[0])))
                    data.append(np.mean(av)) #from st to np
                
                data.sort()
                stmean.append(data[49])
                stmeaner.append([[data[49] - data[15]], [data[83] - data[49]]])
        return medians, stmean, stmeaner, length
    
def bootstrapper(x, y_mid, y_up, y_down, x_bids):
        y_values = [[], [], [], [], [], []]
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