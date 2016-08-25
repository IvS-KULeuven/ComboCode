# -*- coding: utf-8 -*-

"""
Several tools to write complex tables.

Author: R. Lombaert

"""

import os

import cc.path
from cc.tools.io import DataIO
from cc.data.instruments.Spire import Spire
from cc.data.instruments.Pacs import Pacs



def writeIntIntTable(filename,stars,trans,instrument='PACS',ddict=dict(),\
                     searchstring='os2_us3',\
                     mark_trans=[],extra_marker=r'\tablefootmark{f}',\
                     blend_mark=r'\tablefootmark{$\star$}',print_summary=0,\
                     sort_freq=1):

    '''
    Write a table with integrated line intensities and their uncertainties for
    multiple stars. 
    
    The table - for now - is only available in LaTeX format.
    
    Can be used for unresolved line strengths from PACS and SPIRE measurements.
    
    @param filename: The filename of the to be written table.
    @type filename: string
    @param stars: The stars for which the table is created.
    @type stars: list[string]
    @param trans: The transitions for which the integrated intensities are 
                  selected from the PACS line fit results of each star. 
    @type trans: list[Transition()]
    
    @keyword instrument: The name of the instrument for which this is done. For
                         now only 'PACS' and 'SPIRE'.
                         
                         (default: 'PACS')
    @type instrument: str
    @keyword ddict: The data objects for PACS or SPIRE for each star. If not 
                    given, they are generated automatically with path_linefit 
                    'lineFit', oversampling 6 (for PACS) or 4 (for SPIRE), and 
                    resolution of 0.04 (for SPIRE).
                    
                    (default: None)
    @type ddict: dict(Instrument())
    @keyword blend_mark: The marker used for blended lines.
    
                         (default: \tablefootmark{$\star$})
    @type blend_mark: string
    @keyword mark_trans: If a subset of transitions has to be marked with an 
                         extra mark, they are included in this list.
                         
                         (default: [])
    @type mark_trans: list
    @keyword extra_marker: The marker used for the subset mark_trans
    
                           (default: \tablefootmark{f})
    @type extra_marker: string
    @keyword searchstring: the searchstring conditional for the auto-search, if 
                           data have not been read into given Pacs objects or 
                           if new Pacs objects are generated.
        
                           (default: 'os2_us3')
    @type searchstring: string
    @keyword sort_freq: Sort the transitions on frequency. Otherwise sort on 
                        wavelength.
    
                        (default: 1)
    @type sort_freq: bool
    @keyword print_summary: Print a summary at the end. 
    
                            (default: 0)
    @type print_summary: bool
    
    '''
    
    #-- Check what instrument is requested
    if ddict:
        all_instr = [dd.instrument for dd in ddict.values()]
        if len(set(all_instr)) > 1: 
            print "Too many instruments defined. Make sure only PACS or SPIRE"+\
                  " are requested."
            return
        instrument = all_instr[0].upper()
        
    if isinstance(stars,str):
        stars = [stars]
    
    #-- Note that keeping track of transitions per star is not necessary. Here,
    #   we look at data, and these are kept track of by the filename which 
    #   includes the star name. There is no confusion possible.
    trans = sorted(trans,\
                   key=lambda x: sort_freq and x.frequency or x.wavelength)
    if set([t.vup for t in trans]) == set([0]):
        no_vib = 1
    else:
        no_vib = 0
               
    #-- Reset the integrated line strength info in the transitions, in case it 
    #   was already set previously. There's no way to be sure for every trans
    #   individually if the match-up has been done before. And in addition, it
    #   needs to be done for every star separately. So play safe, and reset in 
    #   the sample transition list.
    for tr in trans:
        tr.unreso = dict()
        tr.unreso_err = dict()
        tr.unreso_blends = dict()
    for star in stars:
        if not ddict.has_key(star):
            if instrument == 'PACS':
                ddict[star] = Pacs(star,6,path_linefit='lineFit')
            elif instrument == 'SPIRE':
                ddict[star] = Spire(star,resolution=0.04,oversampling=4,\
                                    path_linefit='lineFit')
        ddict[star].setData(searchstring=searchstring)
        for ifn in range(len(ddict[star].data_filenames)):
            ddict[star].intIntMatch(trans,ifn)
    
    istars = [DataIO.getInputData().index(star) for star in stars]
    pstars = [DataIO.getInputData(keyword='STAR_NAME_PLOTS',rindex=istar)
              for istar in istars]
    pstars = [s.replace('_',' ') for s in pstars]
    all_molecs = DataIO.getInputData(keyword='TYPE_SHORT',make_float=0,\
                                     filename='Molecule.dat')
    all_pmolecs = DataIO.getInputData(keyword='NAME_PLOT',make_float=0,\
                                     filename='Molecule.dat')
    inlines = []
    inlines.append('&'.join(['']*(no_vib and 4 or 5)+pstars[:-1]+\
                            [r'%s \\\hline'%pstars[-1]]))
    line_els = [instrument,'Molecule']
    if not no_vib: 
        line_els.append('Vibrational')
    line_els.extend(['Rotational','$\lambda_0$',\
                     r'\multicolumn{%i}{c}{$F_\mathrm{int}$} \\'%len(pstars)])
    inlines.append('&'.join(line_els))
    line_els = ['band','']
    if not no_vib:
        line_els.append('state')
    line_els.extend(['transition',r'$\mu$m',\
                     r'\multicolumn{%i}{c}{(W m$^-2$))} \\\hline'%len(pstars)])
    inlines.append('&'.join(line_els))
    all_bands = ['SLW','SSW','R1B','R1A','B2B','B2A','B3A'] 
    bands = set([ib for v in ddict.values() for ib in v.data_ordernames])
    bands = [ib for ib in all_bands if ib in bands]
    if not sort_freq: bands.reverse()
    line_counter = dict()
    for s in stars:
        line_counter[s] = 0
    for band in bands:
        #inlines.append(r'\multicolumn{4}{c}{PACS Band: %s} & \multicolumn{%i}{c}{} \\\hline'\
                       #%(band,len(pstars)))
        new_band = 1
        for it,t in enumerate(trans):
            #-- Check if there's any actual line strength result for any star
            #   for this particular band; in any filename in this band.
            #   Otherwise, exclude the line from the table. If there's no 
            #   spectrum for this band at all, the line is excluded as well.
            #   In this case, the band will not be added to the table at all.
            #   Just look at the keys of t.unreso. If there's none of the 
            #   given band, then exclude the transition. All stars are 
            #   included in the same Transition() objects.
            all_keys = [k 
                        for k in t.unreso.keys()
                        if band in os.path.split(k)[1].split('_')]
            if not all_keys:
                continue
            
            #-- There's at least one star with a measured line strength in this 
            #   band, so add a column. 
            col1 = all_pmolecs[all_molecs.index(t.molecule.molecule)]
            
            #-- If the band has not had a line added before, add the band now
            if new_band:
                col0 = band
                new_band = 0
            else:
                col0 = ''
                
            #-- Define the input for the table.
            parts = [col0,\
                     col1]
            if not no_vib:
                parts.append(t.makeLabel(return_vib=1))
            parts.extend([t.makeLabel(inc_vib=0),'%.2f'%(t.wavelength*10**4)])
            #-- For every star, add the measured line strength of the 
            #   transition, if available. For now, it is assumed only one line
            #   strength measurement is available per star per band. (This is 
            #   not strictly true, for instance for overlapping line scans, but
            #   the line integration tool makes it impossible to discern 
            #   between multiple measurements of the same line in the same 
            #   band)
            for s in stars:
                #-- Grab the filename available in the transition object for 
                #   the measured line strength. If multiple filenames with the
                #   correct band are available, only the first is taken.
                all_fn = [sfn for sfn in t.unreso.keys()
                          if band in os.path.split(sfn)[1].split('_')\
                              and s in os.path.split(sfn)[1].split('_')]
                
                #-- Possibly, for this star, no datafiles of given band are
                #   present. Then just add no flux measurement and continue to
                #   the next star.
                if not all_fn: 
                    parts.append(r'/')
                    continue
                else:
                    fn = all_fn[0]
                
                #-- The filename is present, and thus should have a line 
                #   strength indicated, or be flagged as a blend. If not, an 
                #   error will be raised, but this should not happen!
                fint,finterr,fintblend = t.getIntIntUnresolved(fn)
                if fint == 'inblend':
                    parts.append('Blended')
                else:
                    line_counter[s] += 1
                    parts.append('%s%s%.2e (%.1f%s)'\
                                 %(t in mark_trans and extra_marker or r'',\
                                   fint<0 and blend_mark or r'',abs(fint),\
                                   finterr*100,r'\%'))
            parts[-1] = parts[-1] + r'\\'
            inlines.append('&'.join(parts))   
        if not new_band and band != bands[-1]: 
            inlines[-1] = inlines[-1] + r'\hline'
    DataIO.writeFile(filename,input_lines=inlines)
    if print_summary:
        print('Summary')
        for s in stars:
            print('%s: %i lines measured'%(s,len(ddict[s].linefit.wave_fit))+\
                  ', of which %i lines have been identified.'%line_counter[s])
