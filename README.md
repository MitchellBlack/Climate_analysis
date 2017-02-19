# Code for analysing and visualising climate data

For a full list of plotting options, use the help command (`-h`)

e.g., `python plot_hist.py -h` provides the following user manual:

```
usage: plot_hist.py [-h] [--infile [INFILE [INFILE ...]]] [--varname VARNAME]
                    [--units UNITS] [--color COLOR] [--lgnd_label LGND_LABEL]
                    [--plot_kde {True,False}] [--binwidth BINWIDTH]
                    [--alpha ALPHA] [--title TITLE] [--ylabel YLABEL]
                    [--xlabel XLABEL] [--ylim MIN MAX] [--xlim MIN MAX]
                    [--y_mnrticks Y_MNRTICKS] [--x_mnrticks X_MNRTICKS]
                    [--lgnd_loc {1,2,3,4,5,6,7,8,9,10}]
                    [--lgnd_size {xx-small,x-small,small,medium,large,x-large,xx-large}]
                    [--label_mean {True,False}]
                    [--plot_threshold PLOT_THRESHOLD] [--ofile OFILE]

required arguments:
  --infile [INFILE [INFILE ...]]
                        Input file name(s) - can use wildcard '*.nc' or infile
                        '.list'
  --varname VARNAME     Input file variable
  --units UNITS         Units for plotting
  --color COLOR         Color of timeseries (inc. hexcode, e.g. ~hash~31a354
  --lgnd_label LGND_LABEL
                        Timeseries label for figure legend
  --plot_kde {True,False}
                        plot kernel density estimate [default: False]

optional arguments:
  -h, --help            show this help message and exit
  --binwidth BINWIDTH   Bin width [default: 1.0]
  --alpha ALPHA         Fill shape transparency [default: 0.3]
  --title TITLE         Figure title. Use '_' instead of ' '. [default: None]
  --ylabel YLABEL       Y-axis label. Use '_' instead of ' '. [default: None]
  --xlabel XLABEL       X-axis label. Use '_' instead of ' '. [default: None]
  --ylim MIN MAX        Y-axis bounds [default: auto]
  --xlim MIN MAX        X-axis bounds [default: auto]
  --y_mnrticks Y_MNRTICKS
                        No. minor ticks on primary y-axis [default: None]
  --x_mnrticks X_MNRTICKS
                        No. minor ticks on x-axis [default: None]
  --lgnd_loc {1,2,3,4,5,6,7,8,9,10}
                        Location of the primary figure legend [default: no
                        legend]
  --lgnd_size {xx-small,x-small,small,medium,large,x-large,xx-large}
                        Size of the legend text [default: medium]
  --label_mean {True,False}
                        Indicate histogram mean values on plot
  --plot_threshold PLOT_THRESHOLD
                        Plot threshold as vertical line [default: None]
  --ofile OFILE         Name of ouput file [default: None]



Usage:
    python plot_hist.py -h

    legend options:
       location:   1 upper right, 2 upper left, 3 lower left, 4 lower right, 5 right, 6 center left,
                   7 center right, 8 lower center, 9 upper center, 10 center, None no legend
    label options:
       keys:       ~openp~ ( , ~closep~ ), ~deg~ '$^\circ$', ~hash~ # 
Author:
  Mitchell Black, mtblack@student.unimelb.edu.au

```
