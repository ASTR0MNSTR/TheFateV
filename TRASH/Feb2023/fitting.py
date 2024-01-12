import pandas as pd
import csv

class Main:

    def __init__(self, mytable, out):
        self.mytable = mytable
        self.out = out

        self.data_mytable = []
        self.data_out = []

    def data_reading(self):
        mytable_data = pd.read_csv(self.mytable)
        out_data= pd.read_csv(self.out)
        print(mytable_data.shape)
        new_data = pd.merge(mytable_data, out_data, how="left", on=["GAMAID"])
        print(new_data.shape)

        cols = list(new_data.columns.values)
        cols.pop(cols.index('GAMAID'))
        cols.pop(cols.index('0I'))
        new_data = new_data[['GAMAID', '0I']+cols]
        new_data.to_csv('my_table_2.csv')


if __name__ == '__main__':
    obj = Main('sample_out_spectra.csv', 'out.csv')
    obj.data_reading()
