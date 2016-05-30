"""
Filename:     plot_returntimes.py
Author:       Mitchell Black, mtblack@student.unimelb.edu.au
Description:  Plot returntimes from n 1-D NetCDF timeseries
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
    
    mand_argnames = ['infile','varname','units','plot_conf','color','lgnd_label','period']
    pio.check_inargs_mand_rpt(plot_dict, args2check = mand_argnames, check_against='infile')
 
    # Initialise figure
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Get attributes for plotting
    conf_range = plot_dict['conf_range']
    marker_size = plot_dict['marker_size'][0]
    direction = plot_dict['direction'][0]
    bsn = plot_dict['bsn'][0]
    highlight_val = plot_dict['highlight_val']
    alpha = plot_dict['alpha'][0]

    zoom = False
    if plot_dict['zoom_fctr']:
        for i in ['zoom_fctr','zoom_xlim','zoom_ylim']:
            if not plot_dict[i]:
                raise ValueError, 'need to specify argument --'+i
        
        zoom = True
        zoom_fctr = plot_dict['zoom_fctr'][0]
        zoom_loc  = plot_dict['zoom_loc'][0]
        zoom_xlim = plot_dict['zoom_xlim']
        zoom_ylim = plot_dict['zoom_ylim']
        zoom_cnrs = plot_dict['zoom_line_cnrs']
        zoom_xyticklbs = plot_dict['zoom_xyticklbs'][0]

        assert zoom_xlim[0] < zoom_xlim[1]
        assert zoom_ylim[0] < zoom_ylim[1]

        axins = zoomed_inset_axes(ax, zoom=zoom_fctr, loc=zoom_loc, borderpad=1.5) 

    # Plot each tseries
    for i in range(len(plot_dict['infile'])):
        infiles    = pio.read_infile_list(plot_dict['infile'][i])
        varname    = plot_dict['varname'][i][0]
        units      = pio.fix_label(plot_dict['units'][i][0])
        plot_conf  = plot_dict['plot_conf'][i][0]
        color      = pio.fix_label(plot_dict['color'][i][0])
        label      = plot_dict['lgnd_label'][i][0]
        period     = plot_dict['period'][i][0]

        pio.print_plotting_status(i,varname, units)      
        
        var_tseries = nio.concat_nctseries2array(infiles,varname,units,ensure_taxis_same='False')

        xvals, yvals = pio.calc_return_times(var_tseries, direction=direction, period=period)
     
        ax.semilogx(xvals,yvals,'ko',color=color,marker='o',markersize=marker_size,markeredgecolor=color,label=label)
    
        if zoom:
            axins.semilogx(xvals,yvals,'ko',color=color,marker='o',markersize=marker_size,markeredgecolor=color,label=label)

        if highlight_val:
            highlight(ax, xvals, yvals, highlight_val, direction, marker_size)
            if zoom:
                highlight(axins, xvals, yvals, highlight_val, direction, marker_size)

        if plot_conf == 'True':
            conf_vals = pio.calc_return_time_confidences(var_tseries, direction=direction, c = [ conf_range[0], conf_range[1] ], bsn = bsn)
            ax.fill_between(xvals, conf_vals[0], conf_vals[1], facecolor=color, alpha = alpha, edgecolor=color)

            if zoom:
                axins.fill_between(xvals, conf_vals[0], conf_vals[1], facecolor=color, alpha = alpha, edgecolor=color)
           

    # Add figure display options
   
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

    markerscale = 6. / marker_size
        
    pio.plot_legend_no_repeats(ax, font_size = lgnd_size, loc = lgnd_loc, markerscale=markerscale)
    
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

    if zoom:
        axins.set_xlim(zoom_xlim[0],zoom_xlim[1])
        axins.set_ylim(zoom_ylim[0],zoom_ylim[1])
        mark_inset(ax, axins, loc1= zoom_cnrs[0], loc2= zoom_cnrs[1], fc="none", ec="0.5") 
       
        if zoom_xyticklbs == 'off':
            axins.set_xticklabels([])
            axins.set_yticklabels([])


    # Save figure to file
    fig.savefig(ofile)



def highlight(ax, xvals, yvals, highlight_val, direction, marker_size):
    if direction == "ge":
        fi = np.interp([highlight_val], yvals[::-1], xvals[::-1])
    elif direction =="le":
        fi = np.interp([highlight_val], yvals, xvals)
    else:
        raise ValueError, 'direction not recognised'

    ax.plot(fi, highlight_val, color='black',marker='*',markersize=8.0,markeredgecolor='black')

if __name__ == '__main__':

    extra_info = """
Usage:
    python plot_returntimes.py -h

    legend options:
       location:   1 upper right, 2 upper left, 3 lower left, 4 lower right, 5 right, 6 center left,
                   7 center right, 8 lower center, 9 upper center, 10 center, None no legend
    label options:
       keys:       ~openp~ ( , ~closep~ ), ~deg~ '$^\circ$', ~hash~ # 
Author:
  Mitchell Black, mtblack@student.unimelb.edu.au
"""
    description='Plot timeseries from NetCDF file'
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
    parser_required.add_argument("--plot_conf", type=str, nargs=1, action='append', choices=['True','False'], help="Plot confidence bounds?")
    parser_required.add_argument("--color", type=str, nargs=1, action='append', help="Color of timeseries (inc. hexcode, e.g. ~hash~31a354")
    parser_required.add_argument("--lgnd_label", type=str, nargs=1, action='append', help="Timeseries label for figure legend")
    parser_required.add_argument("--period", type=int, nargs=1, action='append', help="Normalise return time by period. (no. tvals per year) [default: 1]")
    
    # Optional input arguments [not repeated]
    parser.add_argument("--marker_size", type=float, nargs=1, default=[2.0], help="Marker size [default: 2.0]")
    parser.add_argument("--direction", type=str, nargs=1, default=['ge'], choices=['le','ge'], help="Calculating return times for values exceeding (ge) or below (le) thresholds [default: 'ge']")
    parser.add_argument("--bsn", type=int, nargs=1, default=[1e4], help="bootstrap number [default: 1e4]")
    parser.add_argument("--highlight_val", type=float, nargs=1, default=[None], help="Marker size [default: None]")
    parser.add_argument("--conf_range", type=float, nargs=2, default=(0.05, 0.95), metavar=('MIN', 'MAX'), help="Range for confidence bounds [default: 0.05 0.95]")

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

    parser.add_argument("--zoom_fctr", type=float, nargs=1, default=None, help="Zoom magnification [default: None (no zoom box)]")
    parser.add_argument("--zoom_xlim", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="Zoom x-axis bounds [default: None]")
    parser.add_argument("--zoom_ylim", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="Zoom y-axis bounds [default: None]")
    parser.add_argument("--zoom_line_cnrs", type=int, nargs=2, default=[1,3], metavar=('cnr1', 'cnr2'), help="Specify corners for lines joining zoom location and zoom inset (tl: 1, tr: 2, bl: 3, br: 4) [default: 1, 3]")
    parser.add_argument("--zoom_loc", type=int, nargs=1, default=[4], choices=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], help="Location of zoom inset [default: 1]")
    parser.add_argument("--zoom_xyticklbs", type=str, nargs=1, default=['on'], choices=['off', 'on'], help="Turn on xy tick labels for zoom inset [default: on]")

    parser.add_argument("--ofile", type=str, nargs=1, default=['plot.png'], help="Name of ouput file [default: None]")

    args = parser.parse_args()
    
    main(args)
