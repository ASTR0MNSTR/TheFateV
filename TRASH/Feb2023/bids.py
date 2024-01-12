import csv
import matplotlib.pyplot as plt 
import numpy as np

class Main:
    def __init__(self, database):
        self.database = database
        self.data_dict = []
        self.ages = []
        self.UNC = [0, 0, 0, 0, 0, 0, 0]
        self.NOAGN = [0, 0, 0, 0, 0, 0, 0]
        self.YESAGN = [0, 0, 0, 0, 0, 0, 0]
        self.NOEL = [0, 0, 0, 0, 0, 0, 0]
        self.NDA = [0, 0, 0, 0, 0, 0, 0]
        self.HBTRICK = [0, 0, 0, 0, 0, 0, 0]
        self.SAMPLE = [0, 0, 0, 0, 0, 0, 0]

        self.UNC_perc = [0, 0, 0, 0, 0, 0, 0]
        self.NOAGN_perc = [0, 0, 0, 0, 0, 0, 0]
        self.YESAGN_perc = [0, 0, 0, 0, 0, 0, 0]
        self.NOEL_perc = [0, 0, 0, 0, 0, 0, 0]
        self.NDA_perc = [0, 0, 0, 0, 0, 0, 0]
        self.HBTRICK_perc = [0, 0, 0, 0, 0, 0, 0]


    def reading(self):
        with open(self.database, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                self.data_dict.append({'age' : float(row['age']), 'AGN' : row['AGN']})
    
    def counter(self):
        self.ages = [8.9, 9.1, 9.3, 9.5, 9.7, 9.9, 10.1]
        for i in range(len(self.ages)):
            for item in self.data_dict:
                if abs(item['age'] - self.ages[i]) <= 0.1:
                    self.SAMPLE[i] += 1
                    
                    if item['AGN'] == 'NO':
                        self.NOAGN[i] += 1
                    elif item['AGN'] == 'HBTRICK':
                        self.HBTRICK[i] += 1
                    elif item['AGN'] == 'NOEL':
                        self.NOEL[i] += 1
                    elif item['AGN'] == 'YES':
                        self.YESAGN[i] += 1
                    elif item['AGN'] == 'UNC':
                        self.UNC[i] += 1
                    elif item['AGN'] == 'NDA':
                        self.NDA[i] += 1
        
        for j in range(len(self.SAMPLE)):
            try:
                self.NDA_perc[j] = self.NDA[j]/self.SAMPLE[j]
                self.NOEL_perc[j] = (self.NDA[j]+self.NOEL[j])/self.SAMPLE[j]
                self.NOAGN_perc[j] = (self.NDA[j]+self.NOEL[j]+self.NOAGN[j])/self.SAMPLE[j]
                print(self.NOAGN[j]/self.SAMPLE[j])
                self.UNC_perc[j] = (self.NDA[j]+self.NOAGN[j]+self.NOEL[j]+self.UNC[j])/self.SAMPLE[j]
                self.YESAGN_perc[j] = (self.NDA[j]+self.UNC[j]+self.NOAGN[j]+self.NOEL[j]+self.YESAGN[j])/self.SAMPLE[j]
            except:
                pass
                        
    def plotter(self):
        barWidth = 0.25
        br1 = np.arange(len(self.SAMPLE))
 
        # Make the plot
        bar1 = plt.bar(br1, self.YESAGN_perc, color ='midnightblue', width = barWidth,
        edgecolor ='grey', label ='YES')
        plt.bar(br1, self.UNC_perc, color ='springgreen', width = barWidth, edgecolor ='grey', label ='UNC')
        plt.bar(br1, self.NOAGN_perc, color ='mediumvioletred', width = barWidth, edgecolor ='grey', label ='NOAGN')
        plt.bar(br1, self.NOEL_perc, color ='orchid', width = barWidth, edgecolor ='grey', label ='NOEL')
        
 
        # Adding Xticks
        plt.xlabel('Age, Gyr', fontweight ='bold', fontsize = 15)
        plt.ylabel('Number', fontweight ='bold', fontsize = 15)
        plt.xticks([r for r in range(len(self.SAMPLE))],
        ['[8.8-9.0]', '[9.0-9.2]', '[9.2-9.4]', '[9.4-9.6]', '[9.6-9.8]', '[9.8-10.0]', '[10.0-10.2]'])

        k = 0
        for rect in bar1:
            height = rect.get_height()
            plt.text(rect.get_x() + rect.get_width() / 2.0, height, str(self.SAMPLE[k]), ha='center', va='bottom')
            k+=1

        plt.legend()
        plt.savefig('AGN_AGE.png')
        plt.show()
    
if __name__ == '__main__':
    obj = Main('out.csv')
    obj.reading()
    obj.counter()
    obj.plotter()

