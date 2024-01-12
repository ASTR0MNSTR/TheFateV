import csv
import json


def rewrite_csv(path_csv, path_json):
    with open(path_csv) as csvfile:
        reader = csv.DictReader(csvfile)
        dict_general = []
        for row in reader:
            dict_general.append(
                {'index': int(row['index']), 
                    'name': row['name'], 
                    'flux_1': float(row['flux_1']), 
                    'flux_1_er': float(row['flux_1_er']), 
                    'flux_1_c': float(row['flux_1_c']), 
                    'flux_1_c_er': float(row['flux_1_c_er']), 
                    'freq_f': float(row['freq_f']), 
                    'freq_c': float(row['freq_c']), 
                    'd': float(row['d']), 
                    'z': float(row['z']), 
                    'SFR': float(row['SFR']), 
                    'M_H2': float(row['M_H2']), 
                    'M_H2_er': float(row['M_H2_er']),
                    'Td' : float(row['Td'])})

        with open(path_json, "w") as jsonfile:
            json.dump(dict_general, jsonfile, skipkeys=True,
                      indent=2, separators=(',', ': '))
            jsonfile.close()


if __name__ == '__main__':
    rewrite_csv('ring_galaxies.csv',
                'RG1.json')