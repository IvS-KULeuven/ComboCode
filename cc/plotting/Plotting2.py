# -*- coding: utf-8 -*-

"""
Preliminary plotting package in Pylab.

Author: R. Lombaert

"""

import pylab as pl
import math
import types
from scipy import array, zeros
from scipy import argmax

from cc.tools.io import DataIO



def plotTiles(data,dimensions,cfg='',**kwargs):
    
    '''
    Plot data in tiles in a single figure. 
    
    The number of tiles in the x and y direction can be specified by the 
    dimensions keyword. 
    
    The data are then given in the form of a dictionary, in which each entry 
    describes the contents of a single tile. These dictionaries are passed in 
    a list of which the length is x-dim * y-dim. Can be less, but then no data 
    will be plotted in the final # tiles. 
    
    If length is less, and the keytags
    keyword is included, the keytags will be put in the final tile. This only 
    works if all data dicts have the same amount of data lists.
    
    The same line types are used for all tiles. They can be specified as well 
    as left to the default.
    
    General properties can also be passed as extra keywords or as a cfg file, 
    as with the plotCols() method.
    
    @param data: The data to be plotted. This list contains dictionaries that 
                 include all information of a single tile in the plot. 
                 
                 Must be included in the dictionary:
                 
                 x: list of lists giving the x-coords
                 
                 y: list of lists giving the y-coords
                 
                 Possible keywords in the dicts are:
                 
                 xerr: List of arrays of x-errors. if [] no errors are included 
                 on the plot. must have same len() as x, include None in the 
                 list for data that have to be plotted without errors.
                 If upper and lower limits differ, instead of single array for 
                 xerri, a list with two elements can be given as well with two 
                 lists: one for the upper limit and one for the lower limit
                 
                 yerr: List of arrays of y-errors. if [] no errors are included 
                 on the plot. must have same len() as y, include None in the 
                 list for data that have to be plotted without errors.
                 If upper and lower limits differ, instead of single array for 
                 yerri, a list with two elements can be given as well with two 
                 lists: one for the upper limit and one for the lower limit
                 
                 labels: (labels, string, x-, y-position)
                 for label in a list, default []
                 
                 histoplot: list of indices of data lists to plot as histogram,
                 default []
                 
                 xmin, xmax, ymin and ymax: floats, default None
                 
                 line_labels: (string,x-pos,same type-integer (eg molecule fi),
                 vibrational?)
                 for label indicating emission lines, default []
                 
                 xaxis: name of x axis (TeX enabled), if keyword not included 
                 here, the method's xaxis is used. If not wanted, use empty 
                 string
                 
                 yaxis: name of y axis (TeX enabled), if keyword not included 
                 here, the method's yaxis is used. If not wanted, use empty 
                 string
                 
    @type data: list[dict]
    @param dimensions: The number of tiles in the x and y direction is given:
                       (x-dim,y-dim)
    @type dimensions: tuple(int)
    @keyword cfg: config filename read as a dictionary, can replace any keyword 
                  given to plotCols. Can also be a dictionary itself, in which
                  case no file is read and the kwargs are updated with the 
                  content of cfg
                  
                  (default: '')
    @type cfg: string/dict
    @keyword filename: filename of the figure, without extension and with path, 
                       if None, the plot is just shown and not saved
                       
                       (default: None)
    @type filename: string
    @keyword extension: extension of the plot filename, adds dot if not present
                        If None, three outputfiles are created: png, eps, pdf
                        Multiple extensions can be requested at the same time 
                        through a list.
                        
                        (default: None)
    @type extension: string/list
    @keyword figsize: the size of the figure, default is A4 page ratio
                      
                      (default: (20.*math.sqrt(2.), 20.) )
    @type figsize: tuple(float,float)
    @keyword show_plot: show fig before saving (hit enter to continue to save)
                   
                   (default: 0)
    @type show_plot: bool
    @keyword xaxis: name of x axis (TeX enabled)
                    
                    (default: r'$\lambda\ (\mu m)$')
    @type xaxis: string
    @keyword yaxis: name of y axis (TeX enabled)
                    
                    (default: r'$F_\\nu (Jy)$')
    @type yaxis: string
    @keyword fontsize_key: fontsize of the keys
                           
                           (default: 16)
    @type fontsize_key: int
    @keyword fontsize_axis: fontsize axis labels
                            
                            (default: 24)
    @type fontsize_axis: int
    @keyword fontsize_title: fontsize title
                            
                            (default: 26)
    @type fontsize_title: int
    @keyword fontsize_label: fontsize of the labels
                             
                             (default: 12)
    @type fontsize_label: int
    @keyword fontsize_ticklabels: fontsize of axis tick labels
                                  
                                  (default: 20)
    @type fontsize_ticklabels: int
    @keyword bold_ticklabels: boldface axis tick labels
                              
                              (default: 0)
    @type bold_ticklabels: bool
    @keyword size_ticklines: The relative size of the tick lines
                              
                             (default: 10)
    @type size_ticklines: float
    @keyword legend_numpoints: Number of points in the legend lines. 
    
                               (default: 1)
    @type legend_numpoints: int
    @keyword keytags: if default no keys, else a key for every dataset in y in 
                      all tiles. Is included in the final tile, so len of dicts
                      in data has to be x-dim * y-dim - 1
                      
                      (default: [])
    @type keytags: list[string] 
    @keyword line_label_types: line types for line labels if requested. 
                               Overridden if number doesn't match number of 
                               different labels requested. 
                               
                               (default: [])
    @type line_label_types: list[string]
    @keyword line_label_color: give different color to a linelabel based on 
                               an integer. Black if off. If line_label_types
                               are given, this keyword is ignored.
                               
                               (default: 0)
    @type line_label_color: bool
    @keyword line_label_lines: put vertical lines where linelabels occur 
                               will be put at the x-position
                               
                               (default: 0)
    @type line_label_lines: bool
    @keyword baseline_line_labels: linelabels at the baseline instead of at 
                                   the bottom
                                   
                                   (default: 0)
    @type baseline_line_labels: bool
    @keyword no_line_label_text: Don't show text for line labels, only lines are
                                 shown if requested.
                               
                                 (default: 0)
    @type no_line_label_text: bool
    @keyword short_label_lines: The label lines are short and at the top of plot
                                
                                (default:0)
    @type short_label_lines: bool   
    @keyword line_label_spectrum: linelabels are set at the top and bottom of 
                                  the spectrum, as for PACS spectra. If 2, all
                                  labels are set at the top.
                                  
                                  (default: 0)
    @type line_label_spectrum: bool
    @keyword line_label_linewidth: The line width of line label lines.
    
                                   (default: 2)
    @type line_label_linewidth: int
    @keyword line_label_dashedvib: Use dashed lines for vibrational transitions
                                   in line label lines. Only applied if the 
                                   ground state label line is a full line.
                             
                                   (default: 0)
    @type line_label_dashedvib: bool
    @keyword linewidth: width of all the lines in the plot
                        
                        (default: 1)
    @type linewidth: int
    @keyword thick_lw_data: Use twice the linewidth for data points (through
                            histoplot keyword).
                            
                            (default: 0)
    @type thick_lw_data: bool
    @keyword xlogscale: set logarithmic scale of x-axis
                        
                        (default: 0)
    @type xlogscale: bool
    @keyword ylogscale: set logarithmic scale of y-axis
                   
                        (default: 0)
    @type ylogscale: bool
    @keyword xmin: if default then autoscaling is done, otherwise min x value.
                   Only used when ddict entries are not defined.
                   
                   (default: None)
    @type xmin: float
    @keyword xmax: if default then autoscaling is done, otherwise max x value
                   Only used when ddict entries are not defined.
                   
                   (default: None)
    @type xmax: float
    @keyword ymin: if default then autoscaling is done, otherwise min y value
                   Only used when ddict entries are not defined.
                   
                   (default: None)
    @type ymin: float
    @keyword ymax: if default then autoscaling is done, otherwise max y value
                   Only used when ddict entries are not defined.
                   
                   (default: None)
    @type ymax: float
    @keyword transparent: for a transparent background
                          
                          (default: 0)
    @type transparent: bool
    @keyword removeYvalues: remove all Y tickmarks on the Y axis
                            
                            (default: 0)
    @type removeYvalues: bool
    @keyword removeXvalues: remove all X tickmarks on the X axis
                            
                            (default: 0)
    @type removeXvalues: bool
    @keyword line_types: if empty, standard line types are used 
                         if an entry is 0, standard line types are used
                         line types are pythonesque ('-k', line style + color) 
                         if list is not empty, take len equal to len(x)
                         
                         (default: [])
                         
                         Colors:
                            - r: red
                            - y: yellow
                            - b: blue
                            - c: cyan
                            - g: green
                            - k: black
                            - m: magenta
                            - [0 and 1]: grayscale 
                              (always as a 3-digit float, e.g. 0.75!!)

                            
                         Styles:
                            - -: line
                            - s: squares
                            - o: large filled circles
                            - .: small filled circles
                            - --: stripes
                            - -.: stripe-point line
                            - x: crosses
                            - +: pluses
                            - p: pentagons
                            - d: filled circle + vertical line
                            - |: vertical line
                            - h,H: different hexagons
                            - *: stars
                            - 2,3,4: "triple crosses" in different orientations
                            - v,>,<,^: Triangles in different orientations
    @type line_types: list[string]
    @keyword wspace: Adjust the width between subplots horizontally
    
                     (default: 0.2)
    @type wspace: float
    @keyword hspace: Adjust the width between subplots vertically
    
                     (default: 0.2)
    @type hspace: float
    @keyword ws_bot: Adjust the amount of white space at the bottom
    
                     (default: 0.1)
    @type ws_bot: float
    @keyword ws_top: Adjust the amount of white space at the top
    
                     (default: 0.9)
    @type ws_top: float
    @keyword ws_left: Adjust the amount of white space at the left
    
                      (default: 0.125)
    @type ws_left: float
    @keyword ws_right: Adjust the amount of white space at the bottom
    
                       (default: 0.9)
    @type ws_right: float
    @keyword landscape: Save the plot in landscape orientation
    
                        (default: 0)
    @type landscape: bool
    @keyword markeredgewidth: Increase the linewidth of marker edges (eg black 
                              circles around colored points) or the thickness 
                              of crosses, dots, etc... Also used for the size 
                              of error bar caps.
                              
                              (default: 1)
    @type markeredgewidth: int
    
    @return: the plotfilename with extension is returned
    @rtype: string
    
    '''

    if cfg:
        if type(cfg) is types.DictType:
            kwargs.update(cfg)
        else:
            kwargs.update(DataIO.readDict(cfg,convert_lists=1,convert_floats=1))
    filename=kwargs.get('filename',None)
    extension=kwargs.get('extension',None)
    figsize=kwargs.get('figsize',(20.*math.sqrt(2.), 20.))
    show_plot=kwargs.get('show_plot',0)
    xaxis=kwargs.get('xaxis',r'$\lambda$ ($\mu$m)')
    yaxis=kwargs.get('yaxis',r'$F_\nu$ (Jy)')
    fontsize_key=kwargs.get('fontsize_key',20)
    fontsize_axis=kwargs.get('fontsize_axis',24)
    fontsize_title=kwargs.get('fontsize_title',26)
    fontsize_label=kwargs.get('fontsize_label',12)
    fontsize_ticklabels=kwargs.get('fontsize_ticklabels',20)
    bold_ticklabels=kwargs.get('bold_ticklabels',0)
    linewidth=kwargs.get('linewidth',1)
    xlogscale=kwargs.get('xlogscale',0)
    ylogscale=kwargs.get('ylogscale',0)
    transparent=kwargs.get('transparent',0)
    removeYvalues=kwargs.get('removeYvalues',0)
    removeXvalues=kwargs.get('removeXvalues',0)
    keytags=kwargs.get('keytags',[])
    legend_numpoints=kwargs.get('legend_numpoints',1)
    line_types=kwargs.get('line_types',[])
    xmin=kwargs.get('xmin',None)
    xmax=kwargs.get('xmax',None)
    ymin=kwargs.get('ymin',None)
    ymax=kwargs.get('ymax',None)
    wspace = kwargs.get('wspace',0.2)
    hspace = kwargs.get('hspace',0.2)
    ws_bot = kwargs.get('ws_bot',0.1)
    ws_top = kwargs.get('ws_top',0.9)
    ws_left = kwargs.get('ws_left',0.125)
    ws_right = kwargs.get('ws_right',0.9)
    baseline_line_labels=kwargs.get('baseline_line_labels',0)
    line_label_types=kwargs.get('line_label_types',[])
    line_label_color=kwargs.get('line_label_color',0)
    line_label_lines=kwargs.get('line_label_lines',0)
    line_label_spectrum=kwargs.get('line_label_spectrum',0)  
    line_label_dashedvib=kwargs.get('line_label_dashedvib',0)
    no_line_label_text=kwargs.get('no_line_label_text',0)
    line_label_linewidth=kwargs.get('line_label_linewidth',2)
    size_ticklines=kwargs.get('size_ticklines',10)
    landscape = kwargs.get('landscape',0)
    short_label_lines = kwargs.get('short_label_lines',0)
    thick_lw_data = kwargs.get('thick_lw_data',0)
    markeredgewidth = kwargs.get('markeredgewidth',1)
    
    xdim = dimensions[0]
    ydim = dimensions[1]
    if keytags:
        if not len(data) <= (xdim*ydim)-1:
            print 'Too many dicts in data to be accomodated by requested ' + \
                  'dimensions, including keytags on the last tile. Aborting.'
            return
        else:
            xlens = array([len(ddict['x']) for ddict in data])
            imaxlen = argmax(xlens)
            maxlen = max(xlens)
            if not len(keytags) == maxlen:
                print 'Number of keytags does not equal number of datasets. '+\
                      'Leaving them out.'
                keytags = []
            else:
                #- copy the data list such that it does not mess with the original
                data = list(data) 
                data.append(dict([('x',[[0]]*maxlen),('y',[[0]]*maxlen)]))
    else:
        if not len(data) <= xdim*ydim:
            print 'Too many dicts in data to be accomodated by requested ' + \
                  'dimensions. Aborting.'
            return
    for i,ddict in enumerate(data):
        if not ddict.has_key('x') or not ddict.has_key('y'):
            print 'Dictionary in data list with index %i does not '%i + \
                  'have either the x or the y keyword. Aborting.'
            return
        if len(ddict['x']) != len(ddict['y']):
            print 'Dictionary in data list with index %i does not '%i + \
                  'the same dimensions for the x and y lists. Aborting.'
        if not ddict.has_key('histoplot'):
            ddict['histoplot'] = []
        if not ddict.has_key('labels'):
            ddict['labels'] = []
        if not ddict.has_key('xmin'):
            ddict['xmin'] = xmin
        if not ddict.has_key('xmax'):
            ddict['xmax'] = xmax
        if not ddict.has_key('ymin'):
            ddict['ymin'] = ymin
        if not ddict.has_key('ymax'):
            ddict['ymax'] = ymax
        if not ddict.has_key('line_labels'):
            ddict['line_labels'] = []
        if not ddict.has_key('xaxis'):
            ddict['xaxis'] = xaxis
        if not ddict.has_key('yaxis'):
            ddict['yaxis'] = yaxis
        if not ddict.has_key('yerr'):
            ddict['yerr'] = [None]*len(ddict['x'])
        if not ddict.has_key('xerr'):
            ddict['xerr'] = [None]*len(ddict['x'])
            
    pl.rc('text', usetex=True)
    pl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    # Figure properties
    figprops = dict(figsize=figsize)
    # New figure
    fig = pl.figure(**figprops) 
    fig.subplots_adjust(wspace=wspace)
    fig.subplots_adjust(hspace=hspace)
    fig.subplots_adjust(bottom=ws_bot)
    fig.subplots_adjust(top=ws_top)
    fig.subplots_adjust(right=ws_right)
    fig.subplots_adjust(left=ws_left)
    
    if transparent: 
        fig.patch.set_alpha(0.5)
    if not line_types:
        linestyles = ['-','-','-','-','-','-','-',\
                    '--','--','--','--','--','--','--',\
                    '.-','.-','.-','.-','.-','.-','.-',\
                    '.','.','.','.','.','.','.',\
                    'x','x','x','x','x','x','x',\
                    'o','o','o','o','o','o','o']
        colors = ['r','b','k','g','m','y','c']
        line_types = [ls + col for ls,col in zip(linestyles,6*colors)]
    for ddict,itile in zip(data,xrange(xdim*ydim)):
        sub = pl.subplot(ydim,xdim,itile+1)
        these_data = []
        [these_data.append([xi,yi,lp,xerri,yerri]) 
             for xi,yi,lp,xerri,yerri in zip(ddict['x'],ddict['y'],line_types,\
                                             ddict['xerr'],ddict['yerr'])
             if list(yi) and yi <> None]
        for index,(xi,yi,lp,xerri,yerri) in enumerate(these_data):
            ls,col = splitLineStyle(lp)
            if index in ddict['histoplot']:
                leg = sub.step(xi,yi,ls,where='mid',color=col,\
                        linewidth=(thick_lw_data and linewidth*2 or linewidth))
            else:
                leg = sub.plot(xi,yi,ls,linewidth= linewidth,color=col)
            if '--' in lp:
                leg[0].set_dashes([15,5])
            if '.-' in lp or '-.' in lp:
                leg[0].set_dashes([15,5,2,5])
            if xerri <> None: 
                try:
                    test = len(xerri[0])
                    lls = xerri[0]
                    if len(xerri) == 1:
                        uls = xerri[0]
                    else:
                        uls = xerri[1]
                except TypeError: 
                    lls = xerri
                    uls = xerri
                for (ll,ul,xii,yii) in zip(lls,uls,xi,yi):
                    if ll != 0 or ul != 0:
                        sub.errorbar(x=[xii],y=[yii],xerr=[[ll],[ul]],\
                                     ecolor=col,\
                                     lolims=ll==0,\
                                     uplims=ul==0,\
                                     fmt=None,\
                                     capsize=5,\
                                     markeredgewidth=markeredgewidth,\
                                     elinewidth=linewidth/2.,\
                                     barsabove=True)#,zorder=zo,alpha=alph)
            if yerri <> None:
                try:
                    test = len(yerri[0])
                    lls = yerri[0]
                    if len(yerri) == 1:
                        uls = yerri[0]
                    else:
                        uls = yerri[1]
                except TypeError: 
                    lls = yerri
                    uls = yerri
                for (ll,ul,xii,yii) in zip(lls,uls,xi,yi):
                    if ll != 0 or ul != 0:
                        sub.errorbar(x=[xii],y=[yii],yerr=[[ll],[ul]],\
                                     ecolor=col,\
                                     lolims=ll==0,\
                                     uplims=ul==0,\
                                     fmt=None,\
                                     capsize=5,\
                                     markeredgewidth=markeredgewidth,\
                                     elinewidth=linewidth/2.,\
                                     barsabove=False)#,zorder=zo,alpha=alph)



        if keytags and itile == len(data)-1:
            prop = pl.matplotlib.font_manager.FontProperties(size=fontsize_key)
            lg = pl.legend(tuple(keytags),loc=(0,0),prop=prop,\
                           numpoints=legend_numpoints)
            lg.legendPatch.set_alpha(0.0)
            ax = pl.gca()
            ax.get_xaxis().set_visible(False)
            ax.get_yaxis().set_visible(False)
            pl.xlim(xmax=-1)
        else:
            sub.autoscale_view(tight=True,scaley=False)
            if xlogscale:
                sub.set_xscale('log')
            if ylogscale:
                sub.set_yscale('log')
            pl.ylabel(ddict['yaxis'],fontsize=fontsize_axis)
            pl.xlabel(ddict['xaxis'],fontsize=fontsize_axis)
            ax = pl.gca()
            if ddict['labels']:
                for s,x,y in ddict['labels']:
                    pl.text(x,y,s,transform=ax.transAxes,fontsize=fontsize_label)
            if ddict['line_labels']:
                y_pos = max([max(yi) for yi in ddict['y'] if yi.size])*1.3
                xmin = min([min(xi) for xi in ddict['x'] if xi.size])
                xmax = max([max(xi) for xi in ddict['x'] if xi.size])
                setLineLabels(line_labels=ddict['line_labels'],y_pos=y_pos,\
                              line_label_spectrum=line_label_spectrum,\
                              line_label_color=line_label_color,\
                              line_label_lines=line_label_lines,\
                              baseline_line_labels=baseline_line_labels,\
                              fontsize_label=fontsize_label,\
                              no_line_label_text=no_line_label_text,\
                              short_label_lines=short_label_lines,\
                              line_label_types=line_label_types,\
                              linewidth=line_label_linewidth,xmin=xmin,\
                              xmax=xmax,\
                              line_label_dashedvib=line_label_dashedvib)
            if removeYvalues:
                ax.set_yticks([])
            if removeXvalues:
                ax.set_xticks([])
            for label in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
                label.set_fontsize(fontsize_ticklabels)
                if bold_ticklabels:
                    label.set_fontweight('bold')
            for tl in sub.get_xticklines() + sub.get_yticklines():
                tl.set_markersize(size_ticklines)
                tl.set_markeredgewidth(1.2)
            for tl in sub.yaxis.get_minorticklines() + sub.xaxis.get_minorticklines():
                tl.set_markersize(size_ticklines/2.)
                tl.set_markeredgewidth(1.2)
            pl.axes(ax) 
            if transparent:
                ax.patch.set_alpha(0.6)
            if ddict['xmin'] <> None:
                pl.xlim(xmin=ddict['xmin']) # min([min(xi) for xi in x])
            if ddict['xmax'] <> None:
                pl.xlim(xmax=ddict['xmax'])
            if ddict['ymin'] <> None:
                pl.ylim(ymin=ddict['ymin']) # min([min(xi) for xi in x])
            if ddict['ymax'] <> None:
                pl.ylim(ymax=ddict['ymax'])
    if filename <> None:
        if extension is None:
            extension = ['.eps','.png','.pdf']
        elif type(extension) is types.StringType:
            extension = [extension]
        for iext,ext in enumerate(extension):
            if ext[0] != '.': 
                extension[iext] = '.' + ext
        [pl.savefig(filename+ext,\
                    orientation=(landscape and 'landscape' or 'portrait'))
         for ext in extension]
    if show_plot or filename is None:
        pl.subplots_adjust(bottom=0.45)
        pl.subplots_adjust(top=0.95)
        pl.subplots_adjust(right=0.7)
        pl.subplots_adjust(left=0.05)
        pl.show()
    pl.close('all')
    
    return filename <> None and filename+extension[0] or None
    
    

def plotCols(x=[],y=[],xerr=[],yerr=[],cfg='',**kwargs):
    
    '''
    Plot a spectrum, with or without a number of subplots.
    
    Pylab can be made to show the plot before saving, and will do so anyway if 
    the filename is not defined, in which case no file will be saved.
    
    The Pylab figure types can be used here, by setting the extension keyword.
    
    @keyword x: List of lists/arrays of x-values, requires inputfilenames if []
                
                (default: [])
    @type x: list[list/arrays/...]
    @keyword y: List of lists/arrays of y-values, requires inputfilenames if []
                
                (default: [])
    @type y: list[list/arrays/...]    
    @keyword xerr: lists with syntax: xerr[xerri] or xerr[[xerri_low,xerri_up]]
                   if [] no errors are included on the plot
                   Must have same len() as x, include None in the list for data
                   that have to be plotted without errors, while other data do
                   have errors associated with them.
                   If upper and lower limits differ, instead of single 
                   array/list for xerri, a list with two elements can be given
                   as well with two lists: one for the upper limit and one for
                   the lower limit
                   For now errorbar color is same as data color
                   
                   (default: [])
    @type xerr: list[list/arrays/...]
    @keyword yerr: List of lists/arrays of y-errs 
                   if [] no errors are included on the plot
                   Must have same len() as y, include None in the list for data
                   that have to be plotted without errors, while other data do
                   have errors associated with them.
                   For now errorbar color is same as data color
                   If upper and lower limits differ, instead of single 
                   array/list for yerri, a list with two elements can be given
                   as well with two lists: one for the upper limit and one for
                   the lower limit
                   
                   (default: [])
    @type yerr: list[list/arrays/...]
    @keyword twiny_x: If a second list of datasets has to be plotted with a 
                      different yaxis, include the second x axis data here. 
                      Only works for two different y axes, with a single x axis.
    
                      (default: [])
    @type twiny_x: list[array]
    @keyword twiny_y: If a second list of datasets has to be plotted with a 
                      different yaxis, include the second y axis data here. 
                      Only works for two different y axes, with a single x axis.
     
                      (default: [])
    @type twiny_y: list[array]
    @keyword inputfiles: list of filenames with two up to six input columns 
                         (x|y) or (x|y|xerr) or (x|y|xerr|yerr) or 
                         (x|y|xerr_low|xerr_up|yerr_low|yerr_up)
                         If no yerr, three columns will do. If yerr, xerr must
                         be included as well. Use 0's if there is no error on x.
                         If upper and lower limits differ, then include 6 
                         columns in total. Cols 3 and 4 are xerr_low and 
                         xerr_up respectively. Cols 5 and 6 are yerr_low and 
                         yerr_up respectively. All 6 must be included in this 
                         case! No combinations possible.
                         
                         (default: [])
    @type inputfiles: list[string]
    @keyword cfg: config filename read as a dictionary, can replace any keyword 
                  given to plotCols. Can also be a dictionary itself, in which
                  case no file is read and the kwargs are updated with the 
                  content of cfg
                  
                  (default: '')
    @type cfg: string/dict
    @keyword filename: filename of the figure, without extension and with path, 
                       if None, the plot is just shown and not saved
                       
                       (default: None)
    @type filename: string
    @keyword extension: extension of the plot filename, adds dot if not present
                        If None, three output files are created: png, eps, pdf
                        Multiple extensions can be requested at the same time 
                        through a list.
                        
                        (default: None)
    @type extension: string/list
    @keyword number_subplots: #subplots in which to plot in vertical direction
                              
                              (default: 1)
    @type number_subplots: int
    @keyword figsize: the size of the figure, default is A4 page ratio
                      
                      (default: (20.*math.sqrt(2.), 20.) )
    @type figsize: tuple(float,float)
    @keyword show_plot: show fig before saving (hit enter to continue to save)
                   
                   (default: 0)
    @type show_plot: bool
    @keyword xaxis: name of x axis (TeX enabled)
                    
                    (default: r'$\lambda\ (\mu m)$')
    @type xaxis: string
    @keyword all_xaxislabels: Include xaxis labels for every subplot, in case
                              number_subplots != 1. If off, only the bottom 
                              subplot gets an xaxis label.
                              
                              (default: 0)
    @type all_xaxislabels: bool
    @keyword yaxis: name of y axis (TeX enabled)
                    
                    (default: r'$F_\\nu (Jy)$')
    @type yaxis: string
    @keyword twinyaxis: name of twin y axis (TeX enabled) if twinx and twiny
                        are not empty. If this is None then the twiny_x and 
                        twin_y are ignored. Also ignored if number_subplots!=1
                    
                        (default: None)
    @type twinyaxis: string
    @keyword plot_title: plot title (TeX enabled), if empty string no title is 
                         used
                        
                         (default: '')
    @type plot_title: string
    @keyword fontsize_key: fontsize of the keys
                           
                           (default: 16)
    @type fontsize_key: int
    @keyword fontsize_axis: fontsize axis labels
                            
                            (default: 24)
    @type fontsize_axis: int
    @keyword fontsize_title: fontsize title
                            
                            (default: 26)
    @type fontsize_title: int
    @keyword fontsize_label: fontsize of the labels
                             
                             (default: 12)
    @type fontsize_label: int
    @keyword fontsize_localized_label: fontsize of the localized labels
                             
                             (default: 12)
    @type fontsize_localized_label: int
    @keyword fontsize_ticklabels: fontsize of axis tick labels
                                  
                                  (default: 20)
    @type fontsize_ticklabels: int
    @keyword bold_ticklabels: boldface axis tick labels
                              
                              (default: 0)
    @type bold_ticklabels: bool
    @keyword legend_numpoints: Number of points in the legend lines. 
    
                               (default: 1)
    @type legend_numpoints: int
    @keyword keytags: if default no keys, else a key for every dataset in y
                      
                      (default: [])
    @type keytags: list[string] 
    @keyword twiny_keytags: as keytags, but only if twinyaxis <> None.
                          
                            (default: [])
    @type twiny_keytags: list[string] 
    @keyword size_ticklines: The relative size of the tick lines
                              
                             (default: 10)
    @type size_ticklines: float
    @keyword labels: string, x-, y-position for label, positions in fig coords:
                    0.0 bottom left, 1,1 upper right
                    
                    (default: [])
    @type labels: list[(string,float,float)]
    @keyword localized_labels: Same as labels, but with the xpos given in  
                               axis coordinates instead of figure coordinates. 
                               Useful for labels you only want to appear in 
                               multi-subplot figures in certain subplots. In 
                               addition the color of the text is given as an 
                               extra entry.
                    
                               (default: [])
    @type localized_labels: list[(string,float,float,string)]   
    @keyword line_labels: (string,x-pos and same type-integer (eg molecule fi),
                          vibrational?)
                          for the label,specifically to indicate emission lines
                          
                          (default: [])
    @type line_labels: list[(string,float,int)]
    @keyword line_label_types: line types for line labels if requested. 
                               Overridden if number doesn't match number of 
                               different labels requested. 
                               
                               (default: [])
    @type line_label_types: list[string]
    @keyword line_label_color: give different color to a linelabel based on 
                               an integer. Black if off. If line_label_types
                               are given, this keyword is ignored.
                               
                               (default: 0)
    @type line_label_color: bool
    @keyword line_label_lines: put vertical lines where linelabels occur 
                               will be put at the x-position
                               
                               (default: 0)
    @type line_label_lines: bool
    @keyword baseline_line_labels: linelabels at the baseline instead of at 
                                   the bottom
                                   
                                   (default: 0)
    @type baseline_line_labels: bool
    @keyword line_label_spectrum: linelabels are set at the top and bottom of 
                                  the spectrum, as for PACS spectra. If 2, all
                                  labels are set at the top.
                                  
                                  (default: 0)
    @type line_label_spectrum: bool
    @keyword no_line_label_text: Don't show text for line labels, only lines are
                                 shown if requested.
                               
                                 (default: 0)
    @type no_line_label_text: bool
    @keyword short_label_lines: The label lines are short and at the top of plot
                                
                                (default:0)
    @type short_label_lines: bool  
    @keyword line_label_linewidth: The line width of line label lines.
    
                                   (default: 2)
    @type line_label_linewidth: int
    @keyword line_label_dashedvib: Use dashed lines for vibrational transitions
                                   in line label lines. Only applied if the 
                                   ground state label line is a full line.
                             
                                   (default: 0)
    @type line_label_dashedvib: bool
    @keyword linewidth: width of all the lines in the plot
                        
                        (default: 1)
    @type linewidth: int
    @keyword markeredgewidth: Increase the linewidth of marker edges (eg black 
                              circles around colored points) or the thickness 
                              of crosses, dots, etc... Also used for the size 
                              of error bar caps.
                              
                              (default: 1)
    @type markeredgewidth: int
    @keyword thick_lw_data: Use twice the linewidth for data points (through
                            histoplot keyword).
                            
                            (default: 0)
    @type thick_lw_data: bool
    @keyword err_linewidth: width of all the errorbars in the plot
                            
                            (default: 1)
    @type err_linewidth: int    
    @keyword key_location: location of the key in the plot in figure coords.
                           Can be 'best' as well to allow python to figure out
                           the best location
                           
                           (default: 'best' )
    @type key_location: tuple/str
    @keyword xlogscale: set logarithmic scale of x-axis
                        
                        (default: 0)
    @type xlogscale: bool
    @keyword ylogscale: set logarithmic scale of y-axis
                   
                        (default: 0)
    @type ylogscale: bool
    @keyword xmin: if default then autoscaling is done, otherwise min x value
                   
                   (default: None)
    @type xmin: float
    @keyword xmax: if default then autoscaling is done, otherwise max x value
                   
                   (default: None)
    @type xmax: float
    @keyword ymin: if default then autoscaling is done, otherwise min y value
                   
                   (default: None)
    @type ymin: float
    @keyword ymax: if default then autoscaling is done, otherwise max y value
                   
                   (default: None)
    @type ymax: float
    @keyword twiny_ymin: if default then autoscaling is done, otherwise min y 
                         value for the twiny axis
                   
                         (default: None)
    @type twiny_ymin: float
    @keyword twiny_ymax: if default then autoscaling is done, otherwise min y 
                         value for the twiny axis
                   
                         (default: None)
    @type twiny_ymax: float
    @keyword transparent: for a transparent background
                          
                          (default: 0)
    @type transparent: bool
    @keyword removeYvalues: remove all Y tickmarks on the Y axis
                            
                            (default: 0)
    @type removeYvalues: bool
    @keyword removeXvalues: remove all X tickmarks on the X axis
                            
                            (default: 0)
    @type removeXvalues: bool
    @keyword histoplot: plot as histogram for data indices given in this list 
                        respective to their positions in the x and y inputlists
                        
                        (default: [])
    @type histoplot: list
    @keyword line_types: if empty, standard line types are used 
                         if an entry is 0, standard line types are used
                         line types are pythonesque ('-k', line style + color) 
                         if list is not empty, take len equal to len(x)
                         
                         (default: [])
                         
                         Colors:
                            - r: red
                            - y: yellow
                            - b: blue
                            - c: cyan
                            - g: green
                            - k: black
                            - m: magenta
                            - [0 and 1]: grayscale 
                              (always as a 3-digit float, e.g. 0.75!!)
                            
                         Styles:
                            - -: line
                            - s: squares
                            - o: large filled circles
                            - .: small filled circles
                            - --: stripes
                            - -.: stripe-point line
                            - x: crosses
                            - +: pluses
                            - p: pentagons
                            - d: filled circle + vertical line
                            - |: vertical line
                            - h,H: different hexagons
                            - *: stars
                            - 2,3,4: "triple crosses" in different orientations
                            - v,>,<,^: Triangles in different orientations
    @type line_types: list[string]
    @keyword markersize: Give different markersize for each x-plot here. If a 
                         single number, the markersize is used for all datasets
                         Default values are set to 5.
                         
                         (default: [])
    @type markersize: list
    @keyword zorder: List that matches the x and y grids. Number gives the 
                     index of the layer that is requested for the respective 
                     datasets. Higher number means the datasets will be on top
                     of lower numbered datasets. Default goes for no preferred
                     ordering.
                     
                     (default: [])
    @type zorder: list    
    @keyword alpha: List that matches the x and y grids. Values give the 
                    opacity of the plotted curve. 
    
                    (default: [])
    @type alpha: list    
    @keyword twiny_line_types: As line_types, but for the twiny data.
    
                              (default: [])
    @type twiny_line_types: list[string]
    @keyword horiz_lines: a list of y-coords to put full black horizontal lines
                        
                          (default: [])
    @type horiz_lines: list[float]
    @keyword horiz_rect: a list of tuples to set a horizontal colored 
                          transparent rectangle on the plot. tuple(y1,y2,color)
                        
                          (default: [])
    @type horiz_rect: list[tuple]
    @keyword vert_lines: a list of x-coords to put full black vertical lines
    
                         (default: [])
    @type vert_lines: list[float]
    @keyword vert_rect: a list of tuples to set a vertical colored transparent 
                        rectangle on the plot. tuple(x1,x2,color)
                        
                        (default: [])
    @type vert_rect: list[tuple]
    @keyword ws_bot: Adjust the amount of white space at the bottom
    
                     (default: 0.1)
    @type ws_bot: float
    @keyword ws_top: Adjust the amount of white space at the top
    
                     (default: 0.9)
    @type ws_top: float
    @keyword ws_left: Adjust the amount of white space at the left
    
                      (default: 0.125)
    @type ws_left: float
    @keyword ws_right: Adjust the amount of white space at the bottom
    
                       (default: 0.9)
    @type ws_right: float
    @keyword hspace: Adjust the width between subplots vertically
    
                     (default: 0.2)
    @type hspace: float
    @keyword landscape: Save the plot in landscape orientation
    
                        (default: 0)
    @type landscape: bool
    @keyword arrows: Draw arrows in a plot. (x0,y0,delta(x),delta(y),width,col,zorder)
                     works like localized_labels.
                     
                     (default: [])
    @type arrows: list[list]
    
    @return: the plotfilename with extension is returned
    @rtype: string

    '''
    
    if cfg:
        if type(cfg) is types.DictType:
            kwargs.update(cfg)
        else:
            kwargs.update(DataIO.readDict(cfg,convert_lists=1,convert_floats=1))
    inputfiles=kwargs.get('inputfiles',[])
    filename=kwargs.get('filename',None)
    extension=kwargs.get('extension',None)
    number_subplots=int(kwargs.get('number_subplots',1))
    figsize=kwargs.get('figsize',(20.*math.sqrt(2.), 20.))
    show_plot=kwargs.get('show_plot',0)
    twiny_x = kwargs.get('twiny_x',[])
    twiny_y = kwargs.get('twiny_y',[])
    xaxis=kwargs.get('xaxis',r'$\lambda$\ ($\mu m$)')
    yaxis=kwargs.get('yaxis',r'$F_\nu$ ($Jy$)')
    twinyaxis=kwargs.get('twinyaxis',None)
    plot_title=kwargs.get('plot_title','')
    fontsize_key=kwargs.get('fontsize_key',16)
    fontsize_axis=kwargs.get('fontsize_axis',24)
    fontsize_title=kwargs.get('fontsize_title',26)
    fontsize_label=kwargs.get('fontsize_label',12)
    fontsize_localized_label=kwargs.get('fontsize_localized_label',12)
    fontsize_ticklabels=kwargs.get('fontsize_ticklabels',20)
    bold_ticklabels=kwargs.get('bold_ticklabels',0)
    size_ticklines=kwargs.get('size_ticklines',10)
    keytags=kwargs.get('keytags',[])
    legend_numpoints=kwargs.get('legend_numpoints',1)
    twiny_keytags=kwargs.get('twiny_keytags',[])
    labels=kwargs.get('labels',[])
    localized_labels=kwargs.get('localized_labels',[])
    linewidth=kwargs.get('linewidth',1)
    markeredgewidth=kwargs.get('markeredgewidth',1)
    err_linewidth=kwargs.get('err_linewdth',1)
    key_location=kwargs.get('key_location','best')
    xlogscale=kwargs.get('xlogscale',0)
    ylogscale=kwargs.get('ylogscale',0)
    xmin=kwargs.get('xmin',None)
    xmax=kwargs.get('xmax',None)
    ymin=kwargs.get('ymin',None)
    ymax=kwargs.get('ymax',None)
    twiny_ymin=kwargs.get('twiny_ymin',None)
    twiny_ymax=kwargs.get('twiny_ymax',None)
    transparent=kwargs.get('transparent',0)
    removeYvalues=kwargs.get('removeYvalues',0)
    removeXvalues=kwargs.get('removeXvalues',0)
    histoplot=kwargs.get('histoplot',[])
    line_types=kwargs.get('line_types',[])
    twiny_line_types=kwargs.get('twiny_line_types',[])
    baseline_line_labels=kwargs.get('baseline_line_labels',0)
    line_labels=kwargs.get('line_labels',[])
    line_label_types=kwargs.get('line_label_types',[])
    line_label_color=kwargs.get('line_label_color',0)
    line_label_lines=kwargs.get('line_label_lines',0)
    line_label_spectrum=kwargs.get('line_label_spectrum',0)
    line_label_linewidth=kwargs.get('line_label_linewidth',2)
    line_label_dashedvib=kwargs.get('line_label_dashedvib',0)
    no_line_label_text=kwargs.get('no_line_label_text',0)
    horiz_lines=kwargs.get('horiz_lines',[])
    vert_lines=kwargs.get('vert_lines',[])
    horiz_rect=kwargs.get('horiz_rect',[])
    vert_rect=kwargs.get('vert_rect',[])
    ws_bot = kwargs.get('ws_bot',0.1)
    ws_top = kwargs.get('ws_top',0.9)
    ws_left = kwargs.get('ws_left',0.125)
    ws_right = kwargs.get('ws_right',0.9)
    hspace = kwargs.get('hspace',0.2)
    landscape = kwargs.get('landscape',0)
    thick_lw_data = kwargs.get('thick_lw_data',0)
    short_label_lines = kwargs.get('short_label_lines',0)
    all_xaxislabels = kwargs.get('all_xaxislabels',0)
    markersize = kwargs.get('markersize',[])
    zorder = kwargs.get('zorder',[])
    alpha = kwargs.get('alpha',[])
    arrows = kwargs.get('arrows',[])
    if inputfiles:
        x,y, xerr, yerr = [],[],[],[]
        read_input = [DataIO.readCols(f) for f in inputfiles]
        for inputfile in read_input:
            x.append(inputfile[0])
            y.append(inputfile[1])
            if len(inputfile) == 2:
                xerr.append(None)
                yerr.append(None)
            elif len(inputfile) == 3:
                xerr.append(inputfile[2])
                yerr.append(None)
            elif len(inputfile) == 4:
                if set(list(inputfile[3])) != set([0]): 
                    yerr.append(inputfile[3])
                else:
                    yerr.append(None)
                if set(list(inputfile[2])) != set([0]): 
                    xerr.append(inputfile[2])
                else:
                    xerr.append(None)
            elif len(inputfile) == 6:
                if set(list(inputfile[2])) != set([0]) \
                        or set(list(inputfile[3])) != set([0]):
                    xerr.append([inputfile[2],inputfile[3]])
                else:
                    xerr.append(None)
                if set(list(inputfile[4])) != set([0]) \
                        or set(list(inputfile[5])) != set([0]):
                    yerr.append([inputfile[4],inputfile[5]])
                else:
                    yerr.append(None)
            else:
                print 'WARNING! Inputfiles have wrong number of columns. ' + \
                      'Aborting.'
                return
    if not list(x) or not list(y):
        print 'WARNING: No x and/or y input defined! Skipping plotting.'
        return
    pl.rc('text', usetex=True)
    pl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    try:
        dummy = len(y[0])
        y = [array(yi) for yi in y]
    except TypeError:
        y = [array(y)]
    try:
        dummy = len(x[0])      
        if len(y) == len(x): 
            x = [array(xi) for xi in x]
        else:
            raise TypeError
    except TypeError:
        x = [array(x) for yi in y]
    keytags = list(keytags)
    #if len(keytags) != len(x) and keytags: keytags = []
    keytags = [k.replace(';',',') for k in keytags]
    labels = [(l1.replace(';',','),l2,l3) for l1,l2,l3 in labels]
    labels = [(l1.replace('{','('),l2,l3) for l1,l2,l3 in labels]
    labels = [(l1.replace('}',')'),l2,l3) for l1,l2,l3 in labels]
    if len(xerr) != len(x) and xerr: xerr = []
    if len(yerr) != len(y) and yerr: yerr = []
    if len(xerr) > len(yerr): yerr = [None]*len(xerr)
    if len(yerr) > len(xerr): xerr = [None]*len(yerr)
    if not xerr and not yerr: xerr,yerr = [None]*len(x), [None]*len(y)
    if number_subplots != 1:
        delta = (x[0][-1] - x[0][0])/float(number_subplots)
    # Figure properties
    figprops = dict(figsize=figsize)
    # Subplot properties
    #adjustprops=dict(left=0.1,bottom=0.1,right=0.97,top=0.93,wspace=0.2,hspace=0.2)         
    # New figure
    fig = pl.figure(**figprops)     
    fig.subplots_adjust(hspace=hspace)
    fig.subplots_adjust(bottom=ws_bot)
    fig.subplots_adjust(top=ws_top)
    fig.subplots_adjust(right=ws_right)
    fig.subplots_adjust(left=ws_left)
    if transparent: 
        fig.patch.set_alpha(0.7)
    # Tunes the subplot layout 
    #fig.subplots_adjust(**adjustprops)                                                                            
    linestyles = ['-','-','-','-','-','-','-',\
                  '--','--','--','--','--','--','--',\
                  ':',':',':',':',':',':',':',\
                  'x','x','x','x','x','x','x',\
                  'o','o','o','o','o','o','o',\
                  ':',':',':',':',':',':',':',\
                  's','s','s','s','s','s','s',\
                  '.','.','.','.','.','.','.',\
                  '-.','-.','-.','-.','-.','-.','-.',\
                  '+','+','+','+','+','+','+',\
                  'h','h','h','h','h','h','h',\
                  'd','d','d','d','d','d','d',\
                  '|','|','|','|','|','|','|',\
                  'p','p','p','p','p','p','p',\
                  '2','2','2','2','2','2','2']
    colors = ['r','b','k','g','m','y','c']
    extra_line_types = [ls + col for ls,col in zip(linestyles,6*colors)]
    line_types,extra_line_types = setLineTypes(x,line_types,extra_line_types)
    if number_subplots != 1: 
        twinyaxis = None
        xmin = None
        xmax = None
    if twinyaxis <> None:
        twiny_line_types,extra_line_types = setLineTypes(twiny_x,\
                                                         twiny_line_types,\
                                                         extra_line_types)
    if type(markersize) is types.ListType and len(markersize) != len(x):
        markersize = [5]*len(x)
    elif not type(markersize) is types.ListType:
        markersize = [markersize]*len(x)
    
    if not zorder or len(zorder) != len(x):
        zorder = range(len(x))
    
    if not alpha or len(alpha) != len(x):
        alpha = [1]*len(x)
        
    for i in xrange(number_subplots):
        sub = pl.subplot(number_subplots, 1, i+1)
        these_data = []
        #- Only if number_subplots is not 1, this splitting
        #- of data is done. Otherwise, normal plots are made.
        if number_subplots != 1:
            [these_data.append(\
                [xi[abs(xi-(xi[0]+(2*(i+1)-1)*delta/2.)) <= delta/2.],\
                 yi[abs(xi-(xi[0]+(2*(i+1)-1)*delta/2.)) <= delta/2.],\
                 lp,\
                 ms,\
                 zo,\
                 alph,\
                 xerri <> None \
                    and (type(xerri[0]) is not types.ListType \
                        and list(xerri[abs(xi-(xi[0]+(2*(i+1)-1)*delta/2.))\
                                       <= delta/2.]) \
                        or [list(xerrii[abs(xi-(xi[0]+(2*(i+1)-1)*delta/2.))\
                                        <= delta/2.])
                            for xerrii in xerri]) \
                    or None,\
                 yerri <> None \
                    and (type(xerri[0]) is not types.ListType \
                        and list(yerri[abs(xi-(xi[0]+(2*(i+1)-1)*delta/2.))\
                                       <= delta/2.]) \
                        or [list(yerrii[abs(xi-(xi[0]+(2*(i+1)-1)*delta/2.))\
                                        <= delta/2.])
                            for yerrii in yerri]) \
                    or None])
             for xi,yi,lp,ms,zo,alph,xerri,yerri in zip(x,y,line_types,markersize,zorder,alpha,xerr,yerr) 
             if list(yi)]
        else:
            [these_data.append([xi,yi,lp,ms,zo,alph,xerri,yerri]) 
             for xi,yi,lp,ms,zo,alph,xerri,yerri in zip(x,y,line_types,markersize,zorder,alpha,xerr,yerr)
             if list(yi)]
        no_err = []
        legends = []
        for index,(xi,yi,lp,ms,zo,alph,xerri,yerri) in enumerate(these_data):
            ls,col = splitLineStyle(lp)
            if index in histoplot:
                leg, = sub.step(xi,yi,ls,where='mid',ms=ms,\
                               linewidth=(thick_lw_data and linewidth*2. \
                                                        or linewidth),\
                               markeredgewidth=markeredgewidth,zorder=zo,\
                               alpha=alph,color=col)
            else:
                leg, = sub.plot(xi,yi,ls,linewidth=linewidth,ms=ms,\
                               markeredgewidth=markeredgewidth,zorder=zo,\
                               alpha=alph,color=col)
            if '--' in lp:
                leg.set_dashes([8,3])
            if '.-' in lp or '-.' in lp:
                leg.set_dashes([8,3,2,3])
            legends.append(leg)
            if xerri <> None: 
                try:
                    test = len(xerri[0])
                    lls = xerri[0]
                    if len(xerri) == 1:
                        uls = xerri[0]
                    else:
                        uls = xerri[1]
                except TypeError: 
                    lls = xerri
                    uls = xerri
                for (ll,ul,xii,yii) in zip(lls,uls,xi,yi):
                    if ll != 0 or ul != 0:
                        sub.errorbar(x=[xii],y=[yii],xerr=[[ll],[ul]],\
                                     ecolor=col,\
                                     lolims=ll==0,\
                                     uplims=ul==0,\
                                     fmt=None,\
                                     capsize=5,\
                                     markeredgewidth=markeredgewidth,\
                                     elinewidth=linewidth/2.,\
                                     barsabove=True,zorder=zo,alpha=alph)
            if yerri <> None:
                try:
                    test = len(yerri[0])
                    lls = yerri[0]
                    if len(yerri) == 1:
                        uls = yerri[0]
                    else:
                        uls = yerri[1]
                except TypeError: 
                    lls = yerri
                    uls = yerri
                for (ll,ul,xii,yii) in zip(lls,uls,xi,yi):
                    if ll != 0 or ul != 0:
                        sub.errorbar(x=[xii],y=[yii],yerr=[[ll],[ul]],\
                                     ecolor=col,\
                                     lolims=ll==0,\
                                     uplims=ul==0,\
                                     fmt=None,\
                                     capsize=5,\
                                     markeredgewidth=markeredgewidth,\
                                     elinewidth=linewidth/2.,\
                                     barsabove=False,zorder=zo,alpha=alph)
        sub.autoscale_view(tight=True,scaley=False)
        if xlogscale:
            sub.set_xscale('log')
        if ylogscale:
            sub.set_yscale('log')
        if line_labels:
            xtemp = [xi for xi in x if xi.size][0]
            if number_subplots != 1:
                xmin = xtemp[0] + i*delta
                xmax = xtemp[0] + (i+1)*delta 
                #these_labels = [label 
                                #for label in line_labels 
                                #if label[1] <= xtemp[0] + (i+1)*delta 
                                    #and label[1] >= xtemp[0] + i*delta]
            else: 
                xmin = xtemp[0]
                xmax = xtemp[-1]
                #these_labels = [label 
                                #for label in line_labels 
                                #if label[1] <= xtemp[-1]
                                    #and label[1] >= xtemp[0]]
            y_pos = max([max(yi) for yi in y if list(yi)])*1.3
            setLineLabels(line_labels=line_labels,y_pos=y_pos,\
                          line_label_spectrum=line_label_spectrum,\
                          line_label_color=line_label_color,\
                          line_label_lines=line_label_lines,\
                          baseline_line_labels=baseline_line_labels,\
                          fontsize_label=fontsize_label,\
                          no_line_label_text=no_line_label_text,\
                          short_label_lines=short_label_lines,\
                          line_label_types=line_label_types,\
                          linewidth=line_label_linewidth,xmin=xmin,xmax=xmax,\
                          line_label_dashedvib=line_label_dashedvib)
        pl.ylabel(yaxis,fontsize=fontsize_axis)
        if i == 0:
            pl.title(plot_title,fontsize=fontsize_title)  
            if keytags:
                keytags = [keys for keys,yi in zip(keytags,y) if list(yi)]
            #if keytags:
                #these_tags = [keys for keys,yi in zip(keytags,y) if list(yi)]
                #lg = pl.legend(tuple(these_tags),loc=key_location,prop=\
                     #pl.matplotlib.font_manager.FontProperties(size=fontsize_key))
                #lg.legendPatch.set_alpha(0.0)
        if i == number_subplots-1 or all_xaxislabels:
            pl.xlabel(xaxis,fontsize=fontsize_axis)
        ax = pl.gca()
        if vert_lines:
            for x_coord in vert_lines:
                ax.axvline(x=x_coord,linewidth=2, ls='--', color='k')
        if horiz_lines:
            for y_coord in horiz_lines:
                ax.axhline(y=y_coord, linewidth=2, ls='--',color='k')
        if horiz_rect:
            for y1,y2,col in horiz_rect:
                ax.axhspan(y1,y2,facecolor=str(col),alpha=0.8)
        if vert_rect:
            for x1,x2,col in vert_rect:
                ax.axvspan(x1,x2,facecolor=str(col),alpha=0.8,zorder=-400)
        if labels:
            for s,xpos,ypos in labels:
                pl.text(xpos,ypos,s,transform=ax.transAxes,\
                        fontsize=fontsize_label)
        if localized_labels:
            ll_xmin = min([sorted(dd[0])[0] for dd in these_data])
            ll_xmax = max([sorted(dd[0])[-1] for dd in these_data])
            for s,xpos,ypos,col in localized_labels:
                if xpos < ll_xmax and xpos > ll_xmin:
                    pl.text((xpos-ll_xmin)/(ll_xmax-ll_xmin),ypos,s,\
                            transform=ax.transAxes,\
                            fontsize=fontsize_localized_label,\
                            color=col)
        if arrows:
            ll_xmin = min([sorted(dd[0])[0] for dd in these_data])
            ll_xmax = max([sorted(dd[0])[-1] for dd in these_data])
            for x0,y0,dx,dy,w,fc,zo in arrows:
                if x0 < ll_xmax and x0 > ll_xmin:
                    arr = pl.Arrow(x0,y0,dx,dy,width=w,zorder=zo)
                    arr.set_facecolor(fc)
                    ax.add_patch(arr)            
                    
        if removeYvalues:
            ax.set_yticks([])
        if removeXvalues:
            ax.set_xticks([])
        for label in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
            label.set_fontsize(fontsize_ticklabels)
            if bold_ticklabels:
                label.set_fontweight('bold')
        for tl in sub.get_xticklines() + sub.get_yticklines():
            tl.set_markersize(size_ticklines)
            tl.set_markeredgewidth(1.2)
        for tl in sub.yaxis.get_minorticklines() + sub.xaxis.get_minorticklines():
            tl.set_markersize(size_ticklines/2.)
            tl.set_markeredgewidth(1.2)
        pl.axes(ax) 
        if transparent:
            ax.patch.set_alpha(0.5)
        if ymin <> None:
            pl.ylim(ymin=ymin) # min([min(xi) for xi in x])
        if ymax <> None:
            pl.ylim(ymax=ymax)
        if twinyaxis <> None:
            twiny_legends = []
            sub2 = sub.twinx()
            twindata = []
            [twindata.extend([xi,yi,lp])
             for xi,yi,lp in zip(twiny_x,twiny_y,twiny_line_types)
             if list(yi)]
            twiny_legends.append(sub2.plot(linewidth=linewidth,*twindata))
            sub2.autoscale_view(tight=True,scaley=False)
            if twiny_ymin is None:
                twiny_ymin = 0.9*min([min(yi) for yi in twiny_y])
            pl.ylim(ymin=twiny_ymin) # min([min(xi) for xi in x])
            if twiny_ymax is None:
                twiny_ymax = 1.1*max([max(yi) for yi in twiny_y])
            pl.ylim(ymax=twiny_ymax)
            if twiny_keytags:
                twiny_keytags = [keys 
                                 for keys,yi in zip(twiny_keytags,twiny_y) 
                                 if list(yi)]
            sub2.set_ylabel(twinyaxis,fontsize=fontsize_axis)
            for label in sub2.xaxis.get_ticklabels() + sub2.yaxis.get_ticklabels():
                label.set_fontsize(fontsize_ticklabels)
                if bold_ticklabels:
                    label.set_fontweight('bold')
            for tl in sub2.get_xticklines() + sub2.get_yticklines():
                tl.set_markersize(size_ticklines)
                tl.set_markeredgewidth(1.2)
            for tl in sub2.yaxis.get_minorticklines() + sub2.xaxis.get_minorticklines():
                tl.set_markersize(size_ticklines/2.)
                tl.set_markeredgewidth(1.2)
        if xmin <> None:
            pl.xlim(xmin=xmin) # min([min(xi) for xi in x])
        if xmax <> None:
            pl.xlim(xmax=xmax)
        if i == 0 and keytags:
            these_legs = legends[:len(keytags)]
            if twiny_keytags:
                these_legs += twiny_legends[:len(twiny_keytags)]
                keytags += twiny_keytags
            prop = pl.matplotlib.font_manager.FontProperties(size=fontsize_key)
            lg = pl.legend(these_legs,keytags,loc=key_location,\
                           numpoints=legend_numpoints,prop=prop)
            lg.legendPatch.set_alpha(0.8)
            lg.set_zorder(max(zorder)+1)
    if filename <> None:
        if extension is None:
            extension = ['.pdf','.png','.eps']
        elif type(extension) is types.StringType:
            extension = [extension]
        for iext,ext in enumerate(extension):
            if ext[0] != '.': 
                extension[iext] = '.' + ext
        [pl.savefig(filename+ext,\
                    orientation=(landscape and 'landscape' or 'portrait'))
         for ext in extension]
    if show_plot or filename is None:
        pl.subplots_adjust(bottom=0.45)
        pl.subplots_adjust(top=0.95)
        pl.subplots_adjust(right=0.7)
        pl.subplots_adjust(left=0.05)
        pl.show()
    pl.close('all')
    return filename <> None and filename+extension[0] or None
    


def setLineTypes(x,line_types,extra_line_types):
    
    '''
    Set the line types for this plot. 
    
    If they are given as input, nothing is changed.
    
    If the input is wrong, they are given by default values. 
    
    Zeroes in the input are also replaced by defaults.
    
    @param x: The datasets to be plotted
    @type x: list[array]
    @param line_types: The input value for the line_types
    @type line_types: list[string]
    @param extra_line_types: The pool of extra line_types, which is destroyed
                             as defaults are used.
    @type extra_line_types: list[string]
    
    @return: The set line types and default pool of line types (possibly 
             smaller than the input value)
    @rtype: (list[string],list[string])
    
    '''
    
    if type(line_types) is types.StringType:
        line_types=[line_types]
    elif not type(line_types) is types.ListType:
        line_types=[]
    if len(line_types) == 1:
        line_types *= len(x)
    elif not line_types or len(line_types) != len(x):
        if len(x) > len(extra_line_types):
            raise IOError('Too many datasets for plotting requested with respect to the amount of extra line types available.')
        line_types = [lp for xi,lp in zip(x,extra_line_types)]
    extra_line_types = [extra  for extra in extra_line_types 
                               if extra not in line_types]
    line_types = [lp and str(lp) or extra_line_types.pop(0) 
                  for lp in line_types]
    return (line_types,extra_line_types)
    
                  
    
def setLineLabels(line_labels,line_label_spectrum,fontsize_label,y_pos,\
                  line_label_color,line_label_lines,baseline_line_labels,\
                  no_line_label_text,short_label_lines,line_label_types,\
                  linewidth,xmin,xmax,line_label_dashedvib):
    
    '''
    Set the line labels in a plot. 
    
    @param line_labels: The line labels to be set including 
                        (labels,x_pos,mol_index,vibrational)
    @type line_labels: list(tuple)
    @param line_label_spectrum: linelabels are set at the top and bottom of 
                                the spectrum, as for PACS spectra. If 2, all 
                                labels are set at the top. 
    @type line_label_spectrum: bool
    @param fontsize_label: fontsize of the labels
    @type fontsize_label: int
    @param y_pos: The y position of the label if line_label_spectrum is False
    @type y_pos: float
    @param line_label_color: Color the line labels and lines
    @type line_label_color: bool
    @param line_label_lines: Draw vertical lines at the label
    @type line_label_lines: bool
    @param baseline_line_labels: Put the line labels at the baseline or the 
                                 bottom of the plot
    @type baseline_line_labels: bool
    @param no_line_label_text: Don't show text for line labels, only lines are
                               shown if requested.
    @type no_line_label_text: bool
    @param short_label_lines: The label lines are short and at the top of plot
    @type short_label_lines: bool
    @param line_label_types: The line types for the line label lines.
    @type line_label_types: list
    @param linewidth: line width of the line label lines.
    @type linewidth: int
    @param xmin: The minimum allowed x value for the line labels
    @type xmin: float
    @param xmax: The maximum allowed x value for the line labels
    @type xmax: float
    @param line_label_dashedvib: Use dashed lines for vibrational transitions 
                                 in line label lines. Only applied if the 
                                 ground state label line is a full line.
    @type line_label_dashedvib: bool
        
    ''' 
    
    if line_label_color:
        colors = ['r','b','k','g','m','c','y']
    else:
        colors = ['k']*7
    linestyles = ['-','-','-','-','-','-','-',\
                  '--','--','--','--','--','--','--',\
                  '.-','.-','.-','.-','.-','.-','.-',\
                  '.','.','.','.','.','.','.',\
                  'x','x','x','x','x','x','x',\
                  'o','o','o','o','o','o','o']
    ltd = dict()
    all_indices = list(set([index for label,x_pos,index,vib in line_labels]))
    for ii,this_i in enumerate(all_indices): 
        if len(all_indices) != len(line_label_types):
            ltd[this_i] = '%s%s'%(linestyles[ii],colors[ii%7])
        else:
            ltd[this_i] = line_label_types[ii]
    
    allowed_labels = [label 
                      for label in line_labels 
                      if label[1] <= xmax and label[1] >= xmin]
    if not no_line_label_text:
        for l,(label,x_pos,index,vib) in enumerate(allowed_labels):
            if line_label_spectrum == 2:
                pl.text(x_pos,pl.axis()[3],label,\
                        ha=line_label_lines and 'right' or 'center',\
                        rotation='vertical',\
                        fontsize=fontsize_label,\
                        color=splitLineStyle(ltd[index])[1])
            elif line_label_spectrum:
                pl.text(x_pos,l%2==0 and pl.axis()[3] or pl.axis()[2],label,\
                        va=baseline_line_labels \
                                    and (l%2==0 and 'top' or 'baseline') \
                                    or (l%2==0 and 'top' or 'bottom'),\
                        ha=line_label_lines and 'right' or 'center',\
                        rotation='vertical',\
                        fontsize=fontsize_label,\
                        color=splitLineStyle(ltd[index])[1])
            else:
                pl.text(x_pos,l%2==0 and pl.axis()[2] or y_pos,label,\
                        va=baseline_line_labels \
                                    and (l%2==0 and 'top' or 'baseline') \
                                    or (l%2==0 and 'top' or 'bottom'),\
                        ha=line_label_lines and 'right' or 'center',\
                        rotation='vertical',\
                        fontsize=fontsize_label,\
                        color=splitLineStyle(ltd[index])[1])
    if line_label_lines:
        for label,x_pos,index,vib in allowed_labels:
            if line_label_dashedvib and vib \
                    and splitLineStyle(ltd[index])[0] == '-':
                ls = '--'
            else: 
                ls = splitLineStyle(ltd[index])[0]
            pl.axvline(x=x_pos,\
                       ymin=short_label_lines and 0.7 or 0,ls=ls,\
                       c=splitLineStyle(ltd[index])[1],\
                       linewidth=linewidth,zorder=-100)
                       
    
    
def splitLineStyle(lp):
    
    '''
    Split a line style string in its color and line type components.
    
    Note that grayscale values should always be given as a four-character 
    string that represents a float between 0 and 1.
    
    @param lp: the line type string, eg '.-k'
    @type lp: string
    @return: two string, the first giving the line type, the second the color
    @rtype: string, string
    
    '''
    
    linetype = ''
    color = ''
    for char in lp:
        if char in ['r','b','k','g','m','y','c','w','0']:
            if char == '0':
                color += lp[lp.index(char):lp.index(char)+4]
            else:
                color += char
            break
    linetype = lp.replace(color,'',1)
    return linetype,color
    
    

def makeHistoPlot(x,y,indices=[]):
    
    '''
    [DEPRECATED]  -- not used by Plotting2 anymore.
    
    Convert x and y data into input for an unbinned 'histogram' plot.
    
    (by Pieter de Groote)
    
    @param x: input x-values
    @type x: list[array/list]
    @param y: input y-values
    @type y: list[array/list]
    @keyword indices:  Make a histogram plot if list is not empty.
                       List holds indices of input data lists for the ones that
                       will be shaped as a histogram.
                       The others will be plotted as normal.
                       
                       default: []
    @type indices: list[int]
    @return: the converted x and y data -- xo, yo
    @rtype: list[array], list[array]
    
    '''
    
    xo, yo = [], []
    indices = [int(i) for i in indices]
    for i,(xi,yi) in enumerate(zip(x,y)):
        if i in indices:
            xoi = pl.zeros(2*len(xi))
            yoi = pl.zeros(2*len(xi))
            xoi[1:-1:2] = xi[:-1]+pl.diff(xi)/2.
            xoi[:-2:2] = xi[:-1]-pl.diff(xi)/2.
            xoi[-1] = xi[-1]+(xi[-1]-xi[-2])/2.
            xoi[-2] = xi[-1]-(xi[-1]-xi[-2])/2.
            yoi[::2] = yi
            yoi[1::2] = yi
            xo.append(xoi)
            yo.append(yoi)
        else:
            xo.append(xi)
            yo.append(yi)
    return xo,yo



def readCfg(cfg):
        
    '''
    Read a cfg file. 
    
    @param cfg: path to the config file. If default,
                the hard-coded default plotting options are used.
    @type cfg: string
            
    @return: The contents of the cfg file are returned.
    @rtype: dict
    
    '''
    
    if cfg:
        cfg_dict = DataIO.readDict(cfg,convert_lists=1,convert_floats=1)
    else:
        cfg_dict = dict()
    
    return cfg_dict

    