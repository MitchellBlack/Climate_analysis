"""
Filename:     plot_tseries.py
Author:       Mitchell Black, mtblack@student.unimelb.edu.au
Description:  Plot 1-D timeseries from NetCDF file(s)
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
import scipy.stats as stats
from matplotlib.dates import DateFormatter


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
    
    mand_argnames = ['infile','varname','units','plot_type','color','lgnd_label']
    pio.check_inargs_mand_rpt(plot_dict, args2check = mand_argnames, check_against='infile')
    
    pio.check_inargs_opt_rpt(plot_dict,'line_style', default = 'solid')    
    pio.check_inargs_opt_rpt(plot_dict,'line_width', default = 1.5)  
    pio.check_inargs_opt_rpt(plot_dict,'marker_type', default = 'o')  
    pio.check_inargs_opt_rpt(plot_dict,'marker_size', default = 1.5)  
    pio.check_inargs_opt_rpt(plot_dict,'yaxis', default = 'primary')  
  
    # Initialise figure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax2 = None

    if ['secondary'] in plot_dict['yaxis']:
        ax2 = ax.twinx()
    
    datasets=[]
    
    # Plot each tseries
    for i in range(len(plot_dict['infile'])):
        infiles    = pio.read_infile_list(plot_dict['infile'][i])
        varname    = plot_dict['varname'][i][0]
        units      = pio.fix_label(plot_dict['units'][i][0])
        line_style = pio.line_alias2style(plot_dict['line_style'][i][0])
        line_color = pio.fix_label(plot_dict['color'][i][0])
        line_width = plot_dict['line_width'][i][0]
        line_label = plot_dict['lgnd_label'][i][0]
        plot_type  = plot_dict['plot_type'][i][0]
        yaxis      = plot_dict['yaxis'][i][0]
        marker_size= plot_dict['marker_size'][i][0]
        ensure_taxis_same = plot_dict['ensure_taxis_same'][0]
        plot_anomaly = plot_dict['plt_anomaly'][0]
        show_bias = plot_dict['show_bias'][0]
        show_corr = plot_dict['show_corr'][0]

        if yaxis == 'primary':
            axn = ax
        else:
            axn = ax2
      
        pio.print_plotting_status(i,varname, units)      
        
        var_tseries,xdates = nio.concat_nctseries2array(infiles,varname,units,return_dates=True,ensure_taxis_same=ensure_taxis_same)

        datasets.append(var_tseries)

        if plot_anomaly == 'True':
            var_tseries = var_tseries - np.nanmean(var_tseries)

        if (len(var_tseries[0]) == 1) and (plot_type == 'tseries_conf'):
                raise ValueError, 'cannot plot confidence bounds for single timeseries'

        if (plot_type in ['tseries_line','tseries_marker']):
            if (plot_type == 'tseries_line'):
                fmt = line_style
            else:
                fmt = plot_dict['marker_type'][i][0]

            for mbr in range(0,var_tseries.shape[0]):
                axn.plot_date(xdates,var_tseries[mbr,:],fmt=fmt,color=line_color,linewidth=line_width,label=line_label, markersize = marker_size)
                 
        elif (plot_type in ['tseries_conf_mean','tseries_conf_median']):
            
            if plot_type == 'tseries_conf_mean':
                thick = np.ma.mean(var_tseries, axis=0)
            elif plot_type == 'tseries_conf_median':
                thick = np.ma.median(var_tseries, axis=0)

            pctl_lwr = gio.calc_nanquantile(var_tseries, q=plot_dict['conf_range'][0], axis=0)
            pctl_upr = gio.calc_nanquantile(var_tseries, q=plot_dict['conf_range'][1], axis=0)
            axn.plot_date(xdates,thick,fmt=line_style,color=line_color,linewidth=line_width,label=line_label)
            axn.fill_between(xdates,pctl_lwr,pctl_upr,facecolor=line_color, alpha=plot_dict['alpha'][0], edgecolor=line_color)

    # Add figure display options
    
    if show_bias == 'True':
        bias = np.ma.mean(datasets[0]) - np.ma.mean(datasets[1])
        lbl_bias = 'Bias: '+str(round(bias,2))
        ax.text(0.05, 0.05, lbl_bias, transform=ax.transAxes, fontsize='large')


    if show_corr == 'True':
        r, p = stats.pearsonr(calc_ens_tmean(datasets[0]),calc_ens_tmean(datasets[1]))   
        lbl_corr = 'r: '+str(round(r,2))+' (p '+str(round(p,2))+')'
        ax.text(0.70, 0.05, lbl_corr, transform=ax.transAxes, fontsize='large')


    if plot_dict['lgnd_loc']:
        markerscale = 6. / marker_size
        pio.plot_legend_no_repeats(ax, ax2, font_size = plot_dict['lgnd_size'][0], loc = plot_dict['lgnd_loc'][0], markerscale=markerscale)

    if plot_dict['xlim']:
        xmin = nio.get_datetime([plot_dict['xlim'][0],])[0]
        xmax = nio.get_datetime([plot_dict['xlim'][1],])[0]
        assert xmax > xmin
        
        ax.set_xlim([xmin, xmax])
    
    if plot_dict['yplim']:
        ax.set_ylim([plot_dict['yplim'][0], plot_dict['yplim'][1]])

    if plot_dict['yslim'] and ax2:
        ax2.set_ylim([plot_dict['yslim'][0], plot_dict['yslim'][1]])

    if plot_dict['xlabel']:
        ax.set_xlabel(pio.fix_label(plot_dict['xlabel'][0]))

    if plot_dict['yplabel']:
        ax.set_ylabel(pio.fix_label(plot_dict['yplabel'][0]))

    if plot_dict['yslabel'] and ax2:
        ax2.set_ylabel(pio.fix_label(plot_dict['yslabel'][0]))
 
    if plot_dict['title']:
        ax.set_title(pio.fix_label(plot_dict['title'][0]))
    
    if plot_dict['yp_mnrticks']:
        pio.set_yaxis_minorticks(ax,plot_dict['yp_mnrticks'][0])
    
    if plot_dict['ys_mnrticks'] and ax2:
        pio.set_yaxis_minorticks(ax2,plot_dict['ys_mnrticks'][0])
    
    if plot_dict['x_mnrticks']:
        pio.set_xaxis_minorticks(ax,plot_dict['x_mnrticks'][0])

    if plot_dict['date_format'][0]:
        ax.xaxis.set_major_formatter(DateFormatter(pio.fix_label(plot_dict['date_format'][0])))

    # Save figure to file
    fig.savefig(plot_dict['ofile'][0])

def calc_ens_tmean(data):
    if len(data.shape) != 1 : 
        data = np.ma.mean(data, axis=0)
    return data

if __name__ == '__main__':

    extra_info = """
Usage:
    python plot_timeseries.py -h

    legend options:
       location:   1 upper right, 2 upper left, 3 lower left, 4 lower right, 5 right, 6 center left,
                   7 center right, 8 lower center, 9 upper center, 10 center, None no legend
    label options:
       keys:       ~openp~ ( , ~closep~ ), ~deg~ '$^\circ$', ~hash~ # 

    date sting format:
                   see https://docs.python.org/2/library/datetime.html#strftime-strptime-behavior

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
    parser_required.add_argument("--plot_type", type=str, nargs=1, action='append', choices=['tseries_line','tseries_conf_mean', 'tseries_conf_median','tseries_marker'], help="type of plot - tseries or tseries_conf")
    parser_required.add_argument("--color", type=str, nargs=1, action='append', help="Color of timeseries (inc. hexcode, e.g. ~hash~31a354")
    parser_required.add_argument("--lgnd_label", type=str, nargs=1, action='append', help="Timeseries label for figure legend")

    # Optional input arguments [repeated as desired] 
    parser.add_argument("--line_style", type=str, nargs=1, action='append', choices=['solid', 'dashdash', 'dotdash', 'dot'], help="Line style [default: 'solid'")
    parser.add_argument("--line_width", type=float, nargs=1, action='append', help="Line width [default: 1.5]")
    parser.add_argument("--marker_type", type=str, nargs=1, action='append', help="Marker type [default: 'o']")
    parser.add_argument("--marker_size", type=float, nargs=1, action='append', help="Marker size [default: 1.5]")
    parser.add_argument("--yaxis", type=str, nargs=1, action='append', choices=['primary','secondary'], help="specify y-axis for plotting [default: 'primary'")
   
    # Optional input arguments [not repeated]
    parser.add_argument("--conf_range", type=float, nargs=2, default=(0.05, 0.95), metavar=('MIN', 'MAX'), help="Range for confidence bounds [default: 0.05 0.95]")
    parser.add_argument("--alpha", type=float, nargs=1, default=[0.3], help="Fill shape transparency [default: 0.3]")
    parser.add_argument("--title", type=str, nargs=1, default=None, help="Figure title. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--yplabel", type=str, nargs=1, default=None, help="Primary y-axis label. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--yslabel", type=str, nargs=1, default=None, help="Secondary y-axis label. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--xlabel", type=str, nargs=1, default=None, help="X-axis label. Use '_' instead of ' '. [default: None]")
    parser.add_argument("--yplim", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="Primary y-axis bounds [default: auto]")
    parser.add_argument("--yslim", type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help="Secondary y-axis bounds [default: auto]")
    parser.add_argument("--xlim", type=str, nargs=2, default=None, metavar=('STRTDATE', 'ENDDATE'), help="X-axis date range [YYYY-MM-DD[ HH:MM:S.S]] [default: auto]")
    parser.add_argument("--yp_mnrticks", type=int, nargs=1, default=None, help="No. minor ticks on primary y-axis [default: None]")
    parser.add_argument("--ys_mnrticks", type=int, nargs=1, default=None, help="No. minor ticks on secondary y-axis [default: None]")
    parser.add_argument("--x_mnrticks", type=int, nargs=1, default=None, help="No. minor ticks on x-axis [default: None]")
    parser.add_argument("--lgnd_loc", type=int, nargs=1, default=[1], choices=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10], help="Location of the primary figure legend [default: no legend]")
    parser.add_argument("--lgnd_size", type=str, nargs=1, default=['medium'], choices=['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], help="Size of the legend text [default: medium]")
    parser.add_argument("--date_format", type=str, nargs=1, default=[None], help="Date format for x-axis, e.g. '~pcnt~Y~endash~~pcnt~M~endash~~pcnt~D' [default: None]")
    parser.add_argument("--ensure_taxis_same", type=str, nargs=1, default=['True'], choices=['True','False'],help="Ensure axis are same for groups of files (inc. time axis) [default: True]")
    parser.add_argument("--plt_anomaly", type=str, nargs=1, default=['False'], choices=['True','False'],help="Plot as anomalies (relative to mean of each set of input files) [default: False]")
    parser.add_argument("--show_bias", type=str, nargs=1, default=['False'], choices=['True','False'],help="Indicate the bias between timeseries input sets #1 and #2 (mean 1 - mean 2) [default: False]")
    parser.add_argument("--show_corr", type=str, nargs=1, default=['False'], choices=['True','False'],help="Indicate correlation between input sets #1 and #2.  [default: False]")
    parser.add_argument("--ofile", type=str, nargs=1, default=['plot.png'], help="Name of ouput file [default: None]")

    args = parser.parse_args()
    
    main(args)
