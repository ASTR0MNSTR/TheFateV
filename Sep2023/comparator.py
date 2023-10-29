import csv

class Main:
    def __init__(self, data_csv, data_txt):
        self.data_csv = data_csv
        self.data_txt = data_txt

        self.data_dict = {}

    def reading_csv(self):
        with open(self.data_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                self.data_dict.update({int(row['GAMAID']) : [row['AGN'], row['SC_WHAN']]})
    
    def reading_txt(self):
        count = 0
        count_err = 0
        count_mis_index = 0
        with open(self.data_txt, 'r') as input:
            while True:
                line = input.readline()
                if not line:
                    break
                line_out = line.strip()
                galaxy_pars = line_out.split()

                GAMAID = int(galaxy_pars[0])
                SC_BPT = galaxy_pars[-2]
                SC_WHAN = galaxy_pars[-1]

                try: 
                    if self.data_dict[GAMAID][0] == SC_BPT and self.data_dict[GAMAID][1] == SC_WHAN:
                        count += 1
                    else:
                        count_err += 1
                except:
                    with open('txt.txt', 'a') as f:
                        line = str(GAMAID) + ' ' + str(galaxy_pars[740]) + ' ' +  str(galaxy_pars[742]) + ' ' +  str(galaxy_pars[743]) + ' ' + SC_BPT + ' ' + SC_WHAN + '\n'
                        f.write(line)
        
        print('True galaxies: ', count)
        print('Err index galaxies: ', count_err)


if __name__ == '__main__':
    obj = Main('E:/backup\ALMA9\GAMA_ETG_OLA.csv', 'E:\LICENSE\ProgsData\main\GAMAforOleg_1.txt')
    obj.reading_csv()
    obj.reading_txt()