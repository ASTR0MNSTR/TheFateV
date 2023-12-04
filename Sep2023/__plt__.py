

def plotter_BPT(self, x, y, up, down, AGN_key, bids):
        yes_temp = []
        UNC_temp = []
        no_temp = []
        noel_temp = []

        yes_temp_age = []
        UNC_temp_age = []
        no_temp_age = []
        noel_temp_age = []

        yes_temp_up = []
        UNC_temp_up = []
        no_temp_up = []
        noel_temp_up = []

        yes_temp_down = []
        UNC_temp_down = []
        no_temp_down = []
        noel_temp_down = []

        tot_age = []
        tot = []
        tot_up = []
        tot_down = []

        for item in self.data_dict:
            X = item[x]
            Y = item[y]
            Y_up = item[up]
            Y_down = item[down]
            AGN = item[AGN_key]
            self.ax4.scatter(X, Y, alpha= 0.4, color = self.color_dict_BPT[AGN][0], marker='.', s = self.color_dict_BPT[AGN][1])
            tot_age.append(X)
            tot.append(Y)
            tot_up.append(Y_up)
            tot_down.append(Y_down)
            #marker = self.color_dict[AGN][2]
            if AGN in ['AGNXY', 'AGNX', 'AGNY']:
                    yes_temp.append(Y)
                    yes_temp_age.append(X)
                    yes_temp_up.append(Y_up)
                    yes_temp_down.append(Y_down)
            elif AGN in ['UNCXY', 'UNCX', 'UNCY']:
                    UNC_temp.append(Y)
                    UNC_temp_age.append(X)
                    UNC_temp_up.append(Y_up)
                    UNC_temp_down.append(Y_down)
            elif AGN in ['SFXY', 'SFX', 'SFY']:
                    no_temp.append(Y)
                    no_temp_age.append(X)
                    no_temp_up.append(Y_up)
                    no_temp_down.append(Y_down)
            elif AGN == 'NOEL':
                    noel_temp.append(Y)
                    noel_temp_age.append(X)
                    noel_temp_up.append(Y_up)
                    noel_temp_down.append(Y_down)
        
        class_list = [[yes_temp_age, yes_temp, 'midnightblue', yes_temp_down, yes_temp_up, 'P'], [UNC_temp_age, UNC_temp, 'springgreen', UNC_temp_down, UNC_temp_up, 'H'], [no_temp_age, no_temp, 'mediumvioletred', no_temp_down, no_temp_up, '*'], [noel_temp_age, noel_temp, 'orchid', noel_temp_down, noel_temp_up, 'o']]

        self.means = []
        self.errs = []
        self.ages = []
        if bids == True:
            age_bids = [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]]
            ages_const = np.arange(8.8, 10.2, 0.2)
        else:
            age_bids = [[10.0, 10.25], [10.25, 10.5], [10.5, 10.75], [10.75, 11], [11, 11.25], [11.25, 11.5]]
            ages_const = np.arange(7, 12, 1)

        for item in class_list:
            X_plot = []
            Y_plot = []
            X, Y, err, length = bootstrapper(item[0], item[1], item[3], item[4], age_bids)
            self.ages = X
            self.means.append(Y)
            self.errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    self.ax4.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[2], fmt=item[5], ms = 12)
                    #self.ax4.text(X[j], Y[j], length[j], c = 'red')
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
            self.ax4.plot(X_plot, Y_plot, alpha = 1, color=item[2])

        #for key in self.color_dict_leg.keys():
        #        self.ax4.scatter(-99, -99, alpha= 1, color = self.color_dict_leg[key][0], marker = self.color_dict_leg[key][2], s = self.color_dict_leg[key][1], label=key)
        self.ax4.set_xticks(ages_const)
        for j, item in enumerate(class_list):
            self.ax4.scatter(-99, -99, alpha = 1, color=item[2], marker=item[5], s = 150, label=self.list_names_BPT[j])
        self.ax4.legend(fontsize="13")

        #self.ax4.set_xlim(8, 10.1)
        #self.ax4.set_ylim(14, 26)
        #self.fig.savefig(name_file + '.pdf')
        #plt.show()
    
def plotter_WHAN(self, x, y, up, down, AGN_key, bids): 
        tot_age = []
        tot = []
        tot_up = []
        tot_down = []

        yes_temp = []
        UNC_temp = []
        no_temp = []
        noel_temp = []

        yes_temp_age = []
        UNC_temp_age = []
        no_temp_age = []
        noel_temp_age = []

        yes_temp_up = []
        UNC_temp_up = []
        no_temp_up = []
        noel_temp_up = []

        yes_temp_down = []
        UNC_temp_down = []
        no_temp_down = []
        noel_temp_down = []

        wAGN_temp = []
        wAGN_temp_age = []
        wAGN_temp_down = []
        wAGN_temp_up = []

        llr_temp = []
        llr_temp_age = []
        llr_temp_down = []
        llr_temp_up = []


        for item in self.data_dict:
            X = item[x]
            Y = item[y]
            Y_up = item[up]
            Y_down = item[down]
            AGN = item[AGN_key]
            self.ax5.scatter(X, Y, alpha= 0.4, color = self.color_dict_WHAN[AGN][0], marker='.', s = self.color_dict_WHAN[AGN][1])
            tot_age.append(X)
            tot.append(Y)
            tot_up.append(Y_up)
            tot_down.append(Y_down)
            #marker = self.color_dict[AGN][2]
            if AGN in ['sAGN']:
                    yes_temp.append(Y)
                    yes_temp_age.append(X)
                    yes_temp_up.append(Y_up)
                    yes_temp_down.append(Y_down)
            elif AGN in ['ELR']:
                    UNC_temp.append(Y)
                    UNC_temp_age.append(X)
                    UNC_temp_up.append(Y_up)
                    UNC_temp_down.append(Y_down)
            elif AGN in ['SF']:
                    no_temp.append(Y)
                    no_temp_age.append(X)
                    no_temp_up.append(Y_up)
                    no_temp_down.append(Y_down)
            elif AGN in ['RG']:
                    noel_temp.append(Y)
                    noel_temp_age.append(X)
                    noel_temp_up.append(Y_up)
                    noel_temp_down.append(Y_down)
            elif AGN in ['LLR']:
                    llr_temp.append(Y)
                    llr_temp_age.append(X)
                    llr_temp_up.append(Y_up)
                    llr_temp_down.append(Y_down)
            elif AGN in ['wAGN']:
                    wAGN_temp.append(Y)
                    wAGN_temp_age.append(X)
                    wAGN_temp_up.append(Y_up)
                    wAGN_temp_down.append(Y_down)
        
        #class_list = [[yes_temp_age, yes_temp, 'midnightblue'], [UNC_temp_age, UNC_temp, 'red'], [no_temp_age, no_temp, 'mediumvioletred'], [noel_temp_age, noel_temp, 'orchid'], [llr_age, llr, 'maroon']]
        class_list = [[yes_temp_age, yes_temp, 'midnightblue', yes_temp_down, yes_temp_up, 'P'], [wAGN_temp_age, wAGN_temp, 'blue', wAGN_temp_down, wAGN_temp_up, 'P'], [no_temp_age, no_temp, 'mediumvioletred', no_temp_down, no_temp_up, '*'], [UNC_temp_age, UNC_temp, 'sandybrown', UNC_temp_down, UNC_temp_up, 'D'], [noel_temp_age, noel_temp, 'chocolate', noel_temp_down, noel_temp_up, 'o'], [llr_temp_age, llr_temp, 'maroon', llr_temp_down, llr_temp_up, 'o']]

        self.means = []
        self.errs = []
        self.ages = []
        if bids == True:
            age_bids = [[8.8, 9.0], [9.0, 9.2], [9.2, 9.4], [9.4, 9.6], [9.6, 9.8], [9.8, 10.0]]
            ages_const = np.arange(8.8, 10.2, 0.2)
        else:
            age_bids = [[10.0, 10.25], [10.25, 10.5], [10.5, 10.75], [10.75, 11], [11, 11.25], [11.25, 11.5]]
            ages_const = np.arange(7, 12, 1)

        for item in class_list:
            X_plot = []
            Y_plot = []
            X, Y, err, length = bootstrapper(item[0], item[1], item[3], item[4], age_bids)
            self.ages = X
            self.means.append(Y)
            self.errs.append(err)
            for j in range(len(X)):
                if Y[j] != -99:
                    self.ax5.errorbar(X[j], Y[j], alpha = 1, xerr=0, yerr= err[j], color=item[2], fmt=item[5], ms = 12)
                    #self.ax5.text(X[j], Y[j], length[j], c = 'red')
                    X_plot.append(X[j])
                    Y_plot.append(Y[j])
            self.ax5.plot(X_plot, Y_plot, alpha = 1, color=item[2])

        #for key in self.color_dict_WHAN.keys():
        #    self.ax5.scatter(-99, -99, alpha= 0.7, color = self.color_dict_WHAN[key][0], marker = self.color_dict_WHAN[key][2], s = self.color_dict_WHAN[key][1], label=key)
        self.ax5.set_xticks(ages_const)
        for j, item in enumerate(class_list):
            self.ax5.scatter(-99, -99, alpha = 1, color=item[2], marker=item[5], s = 150, label=self.list_names_WHAN[j])
        self.ax5.legend(fontsize="13")