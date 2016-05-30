"""
Filename:     plot_hist.py
Author:       Mitchell Black, mtblack@student.unimelb.edu.au
Description:  Plot histograms from n 1-D NetCDF timeseries
"""

import numpy as np
import cdms2
import MV2
import os
import sys
import copy
import argparse
from dateutil.parser import *
import datetime 
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats.kde import gaussian_kde
from matplotlib.dates import DateFormatter
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

EU = os.path.expanduser
sys.path.append(EU("/Users/mtblack/Documents/University.dir/PhD.dir/Code.dir/modules.dir/"))
sys.path.append(EU("/home/mblack/PhD_Code/modules.dir/"))

import modules_netcdf as nio
import modules_general as gio
import modules_plot as pio

def main(inargs):
    
    # Check input arguments
    plot_dict = vars(inargs)    
    del inargs
    
    mand_argnames = ['varname','units','color','lgnd_label','plot_kde']
    pio.check_inargs_mand_rpt(plot_dict, args2check = mand_argnames)
 
    # Initialise figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Get attributes for plotting
    alpha = plot_dict['alpha'][0]
    binwidth = plot_dict['binwidth'][0]
    label_mean = plot_dict['label_mean'][0]
    
    bins=np.arange(-10000, 10000 , binwidth)
    
    x_min=0
    x_max=0
    delta=0
    # Plot each tseries
    for i in range(len(plot_dict['infile'])):
        infiles    = pio.read_infile_list(plot_dict['infile'][i])
        varname    = plot_dict['varname'][i][0]
        units      = pio.fix_label(plot_dict['units'][i][0])
        color      = pio.fix_label(plot_dict['color'][i][0])
        label      = plot_dict['lgnd_label'][i][0]
        plot_kde   = plot_dict['plot_kde'][i][0]
         
        pio.print_plotting_status(i,varname, units)      
        
        var_tseries = nio.concat_nctseries2array(infiles,varname,units,check_axis='False',ensure_taxis_same='False')
        var_tseries = var_tseries.flatten()
        var_tseries = var_tseries[~np.isnan(var_tseries)]

        if plot_kde == 'True':
            kde = gaussian_kde(var_tseries)
            xvals = np.linspace(np.min(var_tseries), np.max(var_tseries), 100)
            ax.plot(xvals,kde(xvals), color=color, alpha=alpha, label=label)
        else:
            ax.hist(var_tseries, bins=bins, normed=1, histtype='stepfilled', facecolor=color, alpha=alpha, label=label)

        if np.nanmin(var_tseries) < x_min :
            x_min = np.nanmin(var_tseries)

        if np.nanmax(var_tseries) > x_max :
            x_max = np.nanmax(var_tseries)

        if label_mean == 'True':
            label = pio.fix_label(label)
            var_mean = np.nanmean(var_tseries)
            var_std  = np.nanstd(var_tseries)
            lbl_mean = 'Mean '+label+': '+str(round(var_mean,2))
            lbl_std  = 'StDev '+label+': '+str(round(var_std,2))
            ax.text(0.05, (0.95 - delta), lbl_mean, transform=ax.transAxes, fontsize='small')
            ax.text(0.05, (0.90 - delta), lbl_std, transform=ax.transAxes, fontsize='small')
            delta = delta + 0.10

    # Add figure display options
    plot_threshold = plot_dict['plot_threshold'][0]
    lgnd_size = plot_dict['lgnd_size'][0]
    lgnd_loc  = plot_dict['lgnd_loc'][0]
    xlim = plot_dict['xlim']
    ylim = plot_dict['ylim']
    xlabel = plot_dict['xlabel'][0]
    ylabel = plot_dict['ylabel'][0]
    title = plot_dict['title'][0]
    y_mnrticks = plot_dict['y_mnrticks']
    x_mnrticks = plot_dict['x_mnrticks']
    ofile = plot_dict['ofile'][0]
  
    pio.plot_legend_no_repeats(ax, font_size = lgnd_size, loc = lgnd_loc)
    
    # set xlim using default min/max
    ax.set_xlim([x_min, x_max])

    if xlim:
        assert xlim[0] < xlim[1]
        ax.set_xlim([xlim[0], xlim[1]])
   
    if ylim:
        assert ylim[0] < ylim[1]
        ax.set_ylim([ylim[0], ylim[1]])

    if xlabel:
        ax.set_xlabel(pio.fix_label(xlabel))

    if ylabel:
        ax.set_ylabel(pio.fix_label(ylabel))

    if title:
        ax.set_title(pio.fix_label(title))
    
    if y_mnrticks:
        pio.set_yaxis_minorticks(ax,y_mnrticks)
    
    if x_mnrticks:
        pio.set_xaxis_minorticks(ax,x_mnrticks)

    if plot_threshold:
        ax.plot([plot_threshold,plot_threshold], [ax.get_ylim()[0],ax.get_ylim()[1]],'--', color='black')

    # Save figure to file
    fig.savefig(ofile)



if __name__ == '__main__':

    extra_info = """
Usage:
    python plot_hist.py -h

    legend options:
       location:   1 upper right, 2 upper left, 3 lower left, 4 lower right, 5 right, 6 center left,
                   7 center right, 8 lower center, 9 upper center, 10 center, None no legend
    label options:
       keys:       ~openp~ ( , ~closep~ ), ~deg~ '$^\circ$', ~hash~ # 
Author:
  Mitchell Black, mtblack@student.unimelb.edu.au
"""
    description='Plot histogram from NetCDF file(s)'
    parser = argparse.ArgumentParser(fromfile_prefix_chars='@',
                                     description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.convert_arg_line_to_args = pio.convert_arg_line_to_args

    # Create group for required arguments
    parser_required = parser.add_argument_group('required arguments')
        
    # Mandantory input arguments [repeated as desired]
    parser_required.add_argument("--infile", type=str, nargs='*', action='append', help="Input file name(s) - can use wildcard '*.nc' or infile '.list'")
    parser_required.add_argument("--varname", type=str, nargs=1, action='append', help="Input file variable")
    parser_required.add_argument("--units", type=str, nargs=1, action='append', help="Units for plotting")
    parser_required.add_argument("--color", type=str, nargs=1, action='append', help="Color of timeseries (inc. hexcode, e.g. ~hash~31a354")
    parser_required.add_argument("--lgnd_label", type=str, nargs=1, action='append', help="Timeseries label for figure legend")
    parser_required.add_argument("--plot_kde", type=str, nargs=1, action='append', choices=['True','False'], help="plot kernel density estimate [default: False]")   
    # Optional input arguments [not repeated]
    parser.add_argument("--binwidth", type=float, nargs=1, default=[1.0], help="Bin width [default: 1.0]")
    parser.add_argument("--alpha", type=float, nargs=1, default=[0.3], help="Fill shape transparency [default: 0.3]")
    parser.add_argument("--title", type=str, nargs=1, default=[None], help="Figure title. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--ylabel", type=str, nargs=1, default=[None], help="Y-axis label. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--xlabel", type=str, nargs=1, default=[None], help="X-axis label. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--ylim", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="Y-axis bounds [default: auto]")
    parser.add_argument("--xlim", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="X-axis bounds [default: auto]")
    parser.add_argument("--y_mnrticks", type=int, nargs=1, default=None, help="No. minor ticks on primary y-axis [default: None]")
    parser.add_argument("--x_mnrticks", type=int, nargs=1, default=None, help="No. minor ticks on x-axis [default: None]")
    parser.add_argument("--lgnd_loc", type=int, nargs=1, default=[1], choices=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], help="Location of the primary figure legend [default: no legend]")
    parser.add_argument("--lgnd_size", type=str, nargs=1, default=['medium'], choices=['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], help="Size of the legend text [default: medium]")
    parser.add_argument("--label_mean", type=str, nargs=1, default=['False'], choices=['True', 'False'], help="Indicate histogram mean values on plot")
    parser.add_argument("--plot_threshold", type=float, nargs=1, default=[None], help="Plot threshold as vertical line [default: None]")

    parser.add_argument("--ofile", type=str, nargs=1, default=['plot.png'], help="Name of ouput file [default: None]")

    args = parser.parse_args()
    
    main(args)
