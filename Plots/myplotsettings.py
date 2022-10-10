from matplotlib import rc, rcParams
rc_params = {'axes.labelsize': 18,
             'axes.titlesize': 15,
             'font.size': 15,
             'lines.linewidth' : 2,
             'legend.fontsize': 12,
             'xtick.labelsize': 14,
             'ytick.labelsize': 14
            }
rcParams.update(rc_params)

rc('text.latex', preamble='\\usepackage{txfonts}')
rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='times')
rc('mathtext', default='sf')
rc('lines', markeredgewidth=1)
rc('lines', linewidth=2)
rc('image', cmap='turbo')

myc = ['blue', 'orange', 'green', 'red', 'lawngreen', 'magenta', 'cyan', 'yellow', 'k']
myls = ['-', '--', ':', '-.', (0, (1, 10)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 5, 1, 5, 1, 5)), '-', '--']
mym = ['.', '+', '*', 's', 'x', 'v', 'o', '2', '.', '+']

cornerplot_settings = {'levels':[0.1,0.5,0.9], 'smooth':0.9, 'plot_datapoints':False, 'plot_density':False, 'fill_contours':True, 'plot_contours':True}
