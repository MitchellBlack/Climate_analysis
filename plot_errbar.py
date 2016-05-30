"""
Filename:     plot_box.py
Author:       Mitchell Black, mtblack@student.unimelb.edu.au
Description:  Plot boxplots from n 1-D NetCDF timeseries
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
    
    mand_argnames = ['varname','units','color','xlabel']
    pio.check_inargs_mand_rpt(plot_dict, args2check = mand_argnames)
 
    # Initialise figure
    alpha = plot_dict['alpha'][0]
    grid  = plot_dict['grid'][0]
    ylim = plot_dict['ylim']
    ylabel = plot_dict['ylabel'][0]
    title = plot_dict['title'][0]
    y_mnrticks = plot_dict['y_mnrticks']
    x_mnrticks = plot_dict['x_mnrticks']
    ofile = plot_dict['ofile'][0]
 
    fig = plt.figure()
    ax = fig.add_subplot(111)
    labels = []
    xvals  = []

    for i in range(len(plot_dict['infile'])):
        infiles    = pio.read_infile_list(plot_dict['infile'][i])
        varname    = plot_dict['varname'][i][0]
        units      = pio.fix_label(plot_dict['units'][i][0])
        color      = pio.fix_label(plot_dict['color'][i][0])
        label      = pio.fix_label(plot_dict['xlabel'][i][0])
        plot_type   = plot_dict['plot_type'][0]

        pio.print_plotting_status(i,varname, units)      
        var_tseries = nio.concat_nctseries2array(infiles,varname,units,ensure_taxis_same='False')
        var_tseries = var_tseries.flatten()
        
        mean = gio.calc_nanmean(var_tseries, axis=0)
        pctl_lwr = gio.calc_nanquantile(var_tseries, q=plot_dict['conf_range'][0], axis=0)
        pctl_upr = gio.calc_nanquantile(var_tseries, q=plot_dict['conf_range'][1], axis=0)

        labels.append(label)
        xvals.append(i)
        
        if plot_type == 'errbar':
            ax.plot(i,mean,marker='o',markersize=plot_dict['marker_size'][0],markerfacecolor=color,markeredgecolor=color)
            ax.plot([i,i],[pctl_lwr,pctl_upr],c=color,linewidth=plot_dict['line_width'][0])
        elif plot_type == 'points':
            var_tseries.flatten()
            for n in xrange(0,len(var_tseries)):
                ax.plot(i,var_tseries[n],marker='o',markersize=plot_dict['marker_size'][0],markerfacecolor=color,markeredgecolor=color)

    plt.xticks(xvals,labels,rotation=plot_dict['xlabel_rotation'][0])
    ax.tick_params(axis='x', labelsize=plot_dict['xlabel_size'][0])
    ax.set_xlim(-1,len(plot_dict['infile']))
    plt.subplots_adjust(bottom=0.15)
    
    if grid == 'on':
        ax.yaxis.grid(True)
     
    if ylim:
        assert ylim[0] < ylim[1]
        ax.set_ylim([ylim[0], ylim[1]])

    if ylabel:
        ax.set_ylabel(pio.fix_label(ylabel))

    if title:
        ax.set_title(pio.fix_label(title))
    
    if y_mnrticks:
        pio.set_yaxis_minorticks(ax,y_mnrticks)
    
    if x_mnrticks:
        pio.set_xaxis_minorticks(ax,x_mnrticks)

    # Save figure to file
    fig.savefig(ofile)

if __name__ == '__main__':

    extra_info = """
Usage:
    python plot_box.py -h

    legend options:
       location:   1 upper right, 2 upper left, 3 lower left, 4 lower right, 5 right, 6 center left,
                   7 center right, 8 lower center, 9 upper center, 10 center, None no legend
    label options:
       keys:       ~openp~ ( , ~closep~ ), ~deg~ '$^\circ$', ~hash~ # 
Author:
  Mitchell Black, mtblack@student.unimelb.edu.au
"""
    description='Plot boxplots from NetCDF file(s)'
    parser = argparse.ArgumentParser(description=description,
                                     epilog=extra_info,
                                     argument_default=argparse.SUPPRESS,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    # Create group for required arguments
    parser_required = parser.add_argument_group('required arguments')
        
    # Mandantory input arguments [repeated as desired]
    parser_required.add_argument("--infile", type=str, nargs='*', action='append', help="Input file name(s) - can use wildcard '*.nc' or infile '.list'")
    parser_required.add_argument("--varname", type=str, nargs=1, action='append', help="Input file variable")
    parser_required.add_argument("--units", type=str, nargs=1, action='append', help="Units for plotting")
    parser_required.add_argument("--color", type=str, nargs=1, action='append', help="Color of timeseries (inc. hexcode, e.g. ~hash~31a354")
    parser_required.add_argument("--xlabel", type=str, nargs=1, action='append', help="Timeseries label for figure legend")
   
    # Optional input arguments [not repeated]
    parser.add_argument("--marker_size", type=float, nargs=1, default=[7.0], help="Marker size [default: 7.0]")
    parser.add_argument("--line_width", type=float, nargs=1, default=[1.5], help="Line width [default: 1.5]")
    parser.add_argument("--conf_range", type=float, nargs=2, default=(0.05, 0.95), metavar=('MIN', 'MAX'), help="Range for confidence bounds [default: 0.05 0.95]")
    parser.add_argument("--alpha", type=float, nargs=1, default=[0.3], help="Fill shape transparency [default: 0.3]")
    parser.add_argument("--title", type=str, nargs=1, default=[None], help="Figure title. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--xlabel_rotation", type=str, nargs=1, default=['horizontal'], help="X-axis label rotation ('horizontal', 'vertical', or degrees). [default: horizontal]")
    parser.add_argument("--xlabel_size", type=float, nargs=1, default=[10], help="X-axis label font sizerotation [default: 10]")
    parser.add_argument("--ylabel", type=str, nargs=1, default=[None], help="Y-axis label. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--ylim", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="Y-axis bounds [default: auto]")
    parser.add_argument("--y_mnrticks", type=int, nargs=1, default=None, help="No. minor ticks on primary y-axis [default: None]")
    parser.add_argument("--x_mnrticks", type=int, nargs=1, default=None, help="No. minor ticks on x-axis [default: None]")
    parser.add_argument("--grid", type=str, nargs=1, default=['on'], choices=['on','off'], help="adding horizontal grid lines [default: on]") 
    parser.add_argument("--ofile", type=str, nargs=1, default=['plot_errbar.png'], help="Name of ouput file [default: None]")
    parser.add_argument("--plot_type", type=str, nargs=1, default=['errbar'], choices=['errbar','points'], help="type of plot [default: errbar]")
 
    args = parser.parse_args()
    
    main(args)
