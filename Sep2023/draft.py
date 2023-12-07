    def plotting_WHAN(self):

        self.ax6 = self.fig.add_subplot(self.gs_top[1,0])
        self.ax7 = self.fig.add_subplot(self.gs_top[1,1], sharey=self.ax6)

        gs_top = plt.GridSpec(2, 1, hspace=0)
        self.fig = plt.figure(figsize=(8,8), tight_layout=True)

        self.ax6 = self.fig.add_subplot(gs_top[0,:])
        self.ax7 = self.fig.add_subplot(gs_top[1,:], sharex=self.ax6)

        self.topaxes = [self.ax7, self.ax6]

        self.ax7.tick_params(top=True, labeltop=False, bottom=True, labelbottom=True, right=True, direction='in')
        self.ax6.tick_params(top=True, labeltop=False, bottom=True, labelbottom=False, right=True, direction='in')
        #for ax in self.topaxes[1:]:
        #plt.setp(ax.get_xticklabels(), visible=False)

        for ax in self.topaxes:    
            ax.set_ylabel(r"$log(EW_{H\alpha})$")
            ax.set_xlim([-2, 2.5])
            ax.set_ylim([-3, 3])
            ax.axhline(y = 0.47712, color = 'black', linestyle='dashed')
            ax.axhline(y = -0.301, color = 'black', linestyle='dotted')
            ax.text(1.5, 0, 'ELR')
            ax.text(1.5, -2, 'LLR')

            X_wAGN = np.arange(-0.4, 2.5, 0.01)
            ax.plot(X_wAGN, 0.77815125+X_wAGN*0, 'black')
            ax.text(-1.5, 2, 'SF')
            ax.text(1.5, 0.5, 'wAGN')
            ax.text(1.5, 2, 'sAGN')

            Y_sAGN = np.arange(0.47712, 3, 0.01)
            ax.plot(-0.4+Y_sAGN*0, Y_sAGN, 'black')
        
        self.ax7.set_yticks(np.arange(-3, 2.1, 1))

        self.ax7.set_xlabel(r"$log(N[II]/H\alpha)$")

        k = 0
        # norm = mpl.colors.Normalize(vmin=8.8,vmax=10.25)
        # choose a colormap
        # c_m = mpl.cm.jet
        # create a ScalarMappable and initialize a data structure
        # self.s_m = mpl.cm.ScalarMappable(cmap=c_m, norm=norm)
        # self.s_m.set_array([])

        for pars in self.flux_er_mod9:
            coord = pars[-1][11]
            x = coord[0]
            pair_x_flags = coord[1]
            y = coord[2]
            pair_y_flags = coord[3]
            age = pars[-2]
            AGN = pars[-1][6]
            SC_WHAN = pars[-1][12]

            if AGN[-1] == '!':
                AGN = AGN[:-1]
                self.ax6.scatter(x, y, color='none', edgecolors='crimson', s=20)
                self.ax7.scatter(x, y, color='none', edgecolors='crimson', s=20)
                

            if SC_WHAN[-1] == '!':
                SC_WHAN = SC_WHAN[:-1]
                self.ax6.scatter(x, y, color='none', edgecolors='black', s=50)
                self.ax7.scatter(x, y, color='none', edgecolors='black', s=50)

            if y >= -3 and y <= 3 and x >= -3.5 and x <= 2:
                k += 1
            if len(pair_x_flags) == 0 and len(pair_y_flags) == 0:
                self.ax7.scatter(x, y, s=self.color_dict[AGN][1], color=self.color_dict[AGN][0], alpha=1, marker=self.color_dict[AGN][2])
                self.ax6.scatter(x, y, s=self.cd_WHAN[SC_WHAN][1], color=self.cd_WHAN[SC_WHAN][0], marker = '.', alpha=1)
                #if len(pair_x_flags) != 0:
                #    Main.plotting_arrows(
                #        self, self.ax6, x, y, pair_x_flags, [], AGN)
                # self.ax6.scatter(plots[i][0], plots[i][1], color=self.s_m.to_rgba(age), marker="x")
            else:
                Main.plotting_arrows_WHAN(
                    self, self.ax6, self.cd_WHAN, x, y, pair_x_flags, pair_y_flags, SC_WHAN, m_y = 1, m_x = 0.3)
                Main.plotting_arrows_WHAN(
                    self, self.ax7, self.color_dict, x, y, pair_x_flags, pair_y_flags, AGN, m_y = 1, m_x = 0.3)
        
        print(k)

        for key in self.cd_WHAN_leg.keys():
            self.ax6.scatter(-99, -99, alpha= 1, color = self.cd_WHAN_leg[key][0], marker = self.cd_WHAN_leg[key][2], s = self.cd_WHAN_leg[key][1], label=key)

        self.ax7.scatter(-99, -99, alpha= 1, color = 'midnightblue', label='AGN', s = 30, marker='o')
        self.ax7.scatter(-99, -99, alpha= 1, color = 'springgreen', label='UNC', s = 30, marker='o')
        self.ax7.scatter(-99, -99, alpha= 1, color = 'mediumvioletred', label='SF', s = 30, marker='o')
        #self.ax6.scatter(-99, -99, color='none', edgecolors='crimson', s=20, label='Abs. lines BPT')
        #self.ax6.scatter(-99, -99, color='none', edgecolors='black', s=20, label='Abs. lines WHAN')
        #self.ax7.scatter(-99, -99, color='none', edgecolors='crimson', s=20, label='Abs. lines BPT')
        #self.ax7.scatter(-99, -99, color='none', edgecolors='black', s=20, label='Abs. lines WHAN')

        self.ax6.legend(loc=3, fontsize="13")
        self.ax7.legend(loc=3, fontsize="13")
        self.fig.savefig('./FIGURES/WHAN.pdf')
        #self.fig.savefig('WHAN.pdf')

        # plt.show()
