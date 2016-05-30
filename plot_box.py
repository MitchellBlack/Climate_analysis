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
    
    mand_argnames = ['varname','units','color','lgnd_label']
    pio.check_inargs_mand_rpt(plot_dict, args2check = mand_argnames)

    alpha = plot_dict['alpha'][0]
    grid  = plot_dict['grid'][0]
    lgnd_size = plot_dict['lgnd_size'][0]
    lgnd_loc  = plot_dict['lgnd_loc'][0]
    xlim = plot_dict['xlim']
    ylim = plot_dict['ylim']
    xlabel = plot_dict['xlabel'][0]
    ylabel = plot_dict['ylabel'][0]
    title = plot_dict['title'][0]
    y_mnrticks = plot_dict['y_mnrticks']
    x_mnrticks = plot_dict['x_mnrticks']
    whisker = plot_dict['whisker'][0]
    outliers = plot_dict['outliers'][0]
    ofile = plot_dict['ofile'][0]
    widths = plot_dict['box_width'][0]
    show_mean = (plot_dict['show_mean'][0] == 'True')
    xtick_rotation = plot_dict['xtick_rotation'][0]

    whis_override = False   
    if whisker == '1.5':
        whis=1.5
        print "whiskers are 1.5* interquartile range"
    elif whisker == 'range':
        whis='range'
        print "whiskers span min and max values"
    elif whisker == '5-95':
        whis=0
        outliers='off'
        pctle=[5,95]
        whis_override = True
        print "whiskers span 5th to 95th percentiles"
    elif whisker == '10-90':
        whis=0
        outliers='off'
        pctle=[10,90]
        whis_override = True
        print "whiskers span 10th to 90th percentiles"

    # read each set of infiles

    box_data=[]
    box_color=[]
    box_label=[]
    pctle_low=[]
    pctle_high=[]
    pctle_25=[]
    pctle_75=[]

    for i in range(len(plot_dict['infile'])):
        infiles    = pio.read_infile_list(plot_dict['infile'][i])
        varname    = plot_dict['varname'][i][0]
        units      = pio.fix_label(plot_dict['units'][i][0])
        color      = pio.fix_label(plot_dict['color'][i][0])
        label      = pio.fix_label(plot_dict['lgnd_label'][i][0])

        pio.print_plotting_status(i,varname, units)      
        var_tseries = nio.concat_nctseries2array(infiles,varname,units,ensure_taxis_same='False')
        
        if whis_override == True :
            pctle_low.append( np.percentile(var_tseries,pctle[0]) )
            pctle_high.append( np.percentile(var_tseries,pctle[1]) )
            pctle_25.append( np.percentile(var_tseries,25) )
            pctle_75.append( np.percentile(var_tseries,75) )

        var_tseries = var_tseries.flatten()

        box_data.append(var_tseries)
        box_color.append(color)
        box_label.append(label)

    # Initialise figure

    fig = plt.figure()
    ax = fig.add_subplot(111)

    if grid == 'on':
        ax.yaxis.grid(True)
        ax.set_axisbelow(True)

    if outliers == 'on':
        showfliers=True
        print "Box plot showing outliers"
    else:
        showfliers=False
        print "Box plot with outliers hidden"

    box = plt.boxplot(box_data, patch_artist=True, showmeans=show_mean, whis=whis, showfliers=showfliers, widths=widths) # add additional args here, such as whis options

    for patch, color in zip(box['boxes'], box_color):
        patch.set_facecolor(color)
        patch.set_alpha(alpha)

    if whis_override == True:
        wdth = 1./3. * widths
        plt.setp(box['caps'],linewidth=0.0)
        plt.setp(box['medians'], linewidth=2.5, color='black')
        plt.scatter(range(1,len(pctle_low)+1),pctle_low,marker='_')
        plt.scatter(range(1,len(pctle_high)+1),pctle_high,marker='_')
        plt.vlines(range(1,len(pctle_low)+1),pctle_low,pctle_25)
        plt.vlines(range(1,len(pctle_high)+1),pctle_high,pctle_75)
        plt.hlines(pctle_low,np.array(range(1,len(pctle_low)+1))-wdth,np.array(range(1,len(pctle_low)+1))+wdth)
        plt.hlines(pctle_high,np.array(range(1,len(pctle_high)+1))-wdth,np.array(range(1,len(pctle_high)+1))+wdth)

    plt.setp(ax, xticks=[y+1 for y in range(len(box_data))], xticklabels=box_label)

    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=xtick_rotation)
    plt.tight_layout()
    
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
   
    # Optional input arguments [not repeated]
    parser.add_argument("--alpha", type=float, nargs=1, default=[1.0], help="Fill shape transparency [default: 1.0]")
    parser.add_argument("--title", type=str, nargs=1, default=[None], help="Figure title. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--ylabel", type=str, nargs=1, default=[None], help="Y-axis label. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--xlabel", type=str, nargs=1, default=[None], help="X-axis label. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--ylim", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="Y-axis bounds [default: auto]")
    parser.add_argument("--xlim", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="X-axis bounds [default: auto]")
    parser.add_argument("--y_mnrticks", type=int, nargs=1, default=None, help="No. minor ticks on primary y-axis [default: None]")
    parser.add_argument("--x_mnrticks", type=int, nargs=1, default=None, help="No. minor ticks on x-axis [default: None]")
    parser.add_argument("--xtick_rotation", type=float, nargs=1, default=[0.], help="Rotation angle for x tick labels [default: 0.]")
    parser.add_argument("--lgnd_loc", type=int, nargs=1, default=[1], choices=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], help="Location of the primary figure legend [default: no legend]")
    parser.add_argument("--lgnd_size", type=str, nargs=1, default=['medium'], choices=['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], help="Size of the legend text [default: medium]")
    parser.add_argument("--grid", type=str, nargs=1, default=['on'], choices=['on','off'], help="adding horizontal grid lines [default: on]") 
    parser.add_argument("--show_mean", type=str, nargs=1, default=['False'], choices=['True','False'], help="show mean value on plot [default: False]")
    parser.add_argument("--box_width", type=float, nargs=1, default=[0.5], help="width of box [default: 0.5")
    parser.add_argument("--whisker", type=str, nargs=1, default=['5-95'], choices=['1.5','range','5-95','10-90'], help="set reach of whiskers past the first and third quartiles")
    parser.add_argument("--outliers", type=str, nargs=1, default=['on'], choices=['on','off'], help="turn on or off plotting of outliers")
    parser.add_argument("--ofile", type=str, nargs=1, default=['plot_box.png'], help="Name of ouput file [default: None]")

    args = parser.parse_args()
    
    main(args)
