cd_WHAN = {
            'UNC' : ['springgreen', 9, '.'],
            'LLR' : ['maroon', 9, '.'],
            'ELR' : ['sandybrown', 9, '.'],
            'SFG' : ['mediumvioletred', 9, '.'],
            'wAGN' : ['blue', 15, '.'],
            'sAGN' : ['midnightblue', 15, '.'],
            'NER' : ['chocolate', 9, '.']
        }

colors_WHAN = {
    'sAGN' : 'midnightblue',
    'wAGN' : 'blue',
    'UNC' : 'springgreen',
    'SFG' : 'mediumvioletred',
    'LLR' : 'maroon',
    'ELR' : 'sandybrown',
    'NER' : 'chocolate'
}

color_dict_BPT ={
            'AGNXY' : ['midnightblue', 15, '.'],
            'AGNX' : ['dodgerblue', 15, '.'],
            'AGNY' : ['black', 15, '.'],
            'UNCXY' : ['springgreen', 9, '.'],
            'UNCX' : ['darkgreen', 9, '.'],
            'UNCY' : ['limegreen', 9, '.'],
            'SFGXY' : ['mediumvioletred', 9, '.'],
            'SFGX' : ['crimson', 9, '.'],
            'SFGY' : ['fuchsia', 9, '.'],
            'NOEL' : ['silver', 9, '.'],
        }

colors_BPT = {
    'AGNXY' : 'midnightblue', 
    'AGNX' : 'dodgerblue',
    # 'AGNY' : 'black',
    'UNCXY' : 'springgreen',
    'UNCX' : 'darkgreen',
    'UNCY' : 'limegreen', 
    'SFGXY' : 'mediumvioletred', 
    'SFGX' : 'crimson',
    'SFGY' : 'fuchsia',
    'NOEL' : 'silver',
    }

color_dict_leg ={
            'AGN' : ['midnightblue', 20, 'o'],
            'UNC' : ['springgreen', 20, 'o'],
            'SFG' : ['mediumvioletred', 20, 'o'],
            'NOEL' : ['silver', 20, 'o'],
        }

cd_WHAN_leg = {
        #    'NDA' : ['gold', 30, '.'],
            'sAGN' : ['midnightblue',20, 'o'],
            'wAGN' : ['blue', 20, 'o'],
            'UNC' : ['springgreen', 20, 'o'],
            'SFG' : ['mediumvioletred', 20, 'o'],
            'ELR' : ['sandybrown', 20, 'o'],
            'NER' : ['chocolate', 20, 'o'],
            'LLR' : ['maroon', 20, 'o']
        }

list_names_BPT = ['AGN', 'UNC', 'SFG', 'NOEL']
list_names_BPT_1 = ['AGNXY', 'AGNX', 'UNCXY', 'UNCX', 'SFGXY', 'SFGX', 'NOEL']
list_names_WHAN = ['sAGN', 'wAGN', 'SFG','ELR', 'NER', 'LLR', 'ALL']


WHAN_color_plt = {
    'sAGN' : ['gray', 150, 'P'],
    'wAGN' : ['gray', 150, 'o'],
    'SFG' : ['gray', 150, '*'],
    'ELR' : ['gray', 150, 'D'],
    'NER' : ['gray', 150, '^'],
    'LLR' : ['gray', 150, 'v']
}

BPT_color_plt ={
    'AGNXY' : ['gray', 150, 'P'],
    'UNCXY' : ['gray', 150, 'H'],
    'SFGXY' : ['gray', 150, '*'],
}