

class Main:
    
    def __init__(self, extractor, gama):
        self.extractor = extractor
        self.gama = gama

        self.GAMAIDs = []
        self.lines_out = []

    def extraction(self):
        with open(r'E:\backup\backup_BPT\PRZEMO\list_of_probs.txt', 'r') as f:
            lines = f.readlines()
            lines_strip = [line.strip() for line in lines]
            for line in lines_strip:
                list_str = line.split(',')
                for id in list_str:
                    if id != '':
                        self.GAMAIDs.append(int(id))
        print(self.GAMAIDs)

    def read_gama(self):
        count = 0
        with open(self.gama, 'r') as input:
            header = input.readline()
            self.lines_out.append(header + '\n')
            while True:
                line = input.readline()
                if not line:
                    break
                count += 1
                Main.dividing(self, line)
        print(count)

    def dividing(self, line):
        line = line.strip()
        try:
            GAMAID = int(line.split()[0])
        except:
            GAMAID = -1
        if GAMAID in self.GAMAIDs:
            self.lines_out.append(line + '\n')
        
    def file_out(self):
        with open('PROB_PRZEMO.txt', 'w') as f:
            f.writelines(self.lines_out)

if __name__ == '__main__':
    obj = Main(r'E:\backup\backup_BPT\PRZEMO\list_of_probs.txt', r'E:\LICENSE\ProgsData\main\GAMAv3.txt')
    obj.extraction()
    obj.read_gama()
    obj.file_out()