# -*- coding: utf-8 -*-

"""
Producing Pacs-related output.

Author: R. Lombaert

"""

import types
import os
import subprocess
import cPickle
from glob import glob
from time import gmtime
from scipy import array,argsort,interpolate
import numpy as np

import cc.path
from cc.tools.io import DataIO
from cc.data.instruments.Instrument import Instrument
from cc.tools.io import Database
from cc.data import Data


def compareInts(pp1,pp2):
    
    '''
    Compare integrated line fluxes for an identical target with the same 
    input Hipe line list.
    
    @param pp1: The first Pacs object with the PACS data and integrated line
                strengths.
    @type pp1: Pacs()
    @param pp2: The second Pacs object with the PACS data and integrated line
                strengths.
    @type pp2: Pacs()
    
    @return: The two one-to-one comparable arrays of line fluxes.
    @rtype: (arr,arr)
    
    #pp1 = Pacs.Pacs('waql',6,'codeJun2013',path_linefit='/home/robinl/Data/PACS/waql/lineFit/')
    #pp2 = Pacs.Pacs('waql',6,'codeJun2013',path_linefit='/home/robinl/Data/PACS/waql/lineFit_normalization/')
    #wl, fl1, fl2, fl1err, fl2err = Pacs.compareInts(pp1,pp2)
    #Plotting2.plotCols(x=[wl],y=[fl1/fl2],line_types=['or'],ylogscale=1,ymin=0.75,ymax=1.25,horiz_lines=[0.8,0.9,1.1,1.2])
    '''
    
    fl1 = []
    fl2 = []
    fl1err = []
    fl2err = []
    lf1 = pp1.linefit
    lf2 = pp2.linefit       
    wl = []
    for gi1 in lf1['groupid']:
        for w1 in lf1['wave_in'][lf1['groupid']==gi1]:
            if w1 in lf2['wave_in'][lf2['groupid']==gi1]:
                f1 = lf1['line_flux'][lf1['groupid']==gi1]
                f2 = lf2['line_flux'][lf2['groupid']==gi1]
                f1err = lf1['line_flux_err'][lf1['groupid']==gi1]
                f2err = lf2['line_flux_err'][lf2['groupid']==gi1]
                i1 = list(lf1['wave_in'][lf1['groupid']==gi1]).index(w1)
                i2 = list(lf2['wave_in'][lf2['groupid']==gi1]).index(w1)
                fl1.append(f1[i1])
                fl2.append(f2[i2])
                fl1err.append(f1err[i1])
                fl2err.append(f2err[i2])
                wl.append(w1)
    wl, fl1, fl2, fl1err, fl2err = array(wl), array(fl1), array(fl2), array(fl1err), array(fl2err)
    asor = argsort(wl)
    return (wl[asor],fl1[asor],fl2[asor],fl1err[asor],fl2err[asor])
    


def writeIntIntTable(filename,stars,trans,dpacs=dict(),searchstring='os2_us3',\
                     mark_trans=[],extra_marker=r'\tablefootmark{f}',\
                     blend_mark=r'\tablefootmark{$\star$}',print_summary=0,\
                     sort_freq=1):

    '''
    Write a table with integrated line intensities and their uncertainties for
    multiple stars. 
    
    The table - for now - is only available in LaTeX format.
    
    @param filename: The filename of the to be written table.
    @type filename: string
    @param stars: The stars for which the table is created.
    @type stars: list[string]
    @param trans: The transitions for which the integrated intensities are 
                  selected from the PACS line fit results of each star. 
    @type trans: list[Transition()]
    
    @keyword path_pacs: full path to PACS data folder, excluding star_name
                        Only needed when PACS data objects are created from 
                        scratch.
                        
                        (default: ~/Data/PACS/)
    @type path_pacs: string
    @keyword dpacs: The data objects for PACS for each star. If not given, they
                    are generated automatically with path_linefit 'lineFit', 
                    oversampling 6 and given path_pacs.
                    
                    (default: None)
    @type dpacs: dict(Pacs())
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
    
    if type(stars) is types.StringType:
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
        if not dpacs.has_key(star):
            dpacs[star] = Pacs(star,6,path_linefit='lineFit')
        dpacs[star].setData(searchstring=searchstring)
        for ifn in range(len(dpacs[star].data_filenames)):
            dpacs[star].intIntMatch(trans,ifn)
    
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
    line_els = ['PACS','Molecule']
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
    bands = ['R1B','R1A','B2B','B2A','B3A']
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
            print('%s: %i lines measured'%(s,len(dpacs[s].linefit.wave_fit))+\
                  ', of which %i lines have been identified.'%line_counter[s])



class Pacs(Instrument):
    
    """
    Environment with several tools for Pacs molecular spectroscopy managing.
    
    """
    
    def __init__(self,star_name,oversampling,path=None,redo_convolution=0,\
                 intrinsic=1,path_linefit=''):
        
        '''
        Initializing an instance of Pacs().
        
        @param star_name: Name of the star from Star.dat
        @type star_name: string
        @param oversampling: The PACS instrumental oversampling, for correct
                             convolution of the Sphinx output]
        @type oversampling: int
        
        @keyword path: Output folder in the code's home folder. Used to locate 
                       the PACS database. If None, it is not used (eg for line
                       measurement matching/identification)
                       
                       (default: None)
        @type path: string
        @keyword intrinsic: Use the intrinsic Sphinx line profiles for 
                            convolving with the spectral resolution? Otherwise
                            the beam convolved line profiles are used.
                            
                            (default: 1)
        @type intrinsic: bool
        @keyword redo_convolution: if you want to do the convolution of the 
                                 sphinx models regardless of what's already in 
                                 the database. The pacs id will not change, nor
                                 the entry in the db, and the old convolved 
                                 model will be copied to a backup
                                 
                                 (default: 0)
        @type redo_convolution: bool
        @keyword path_linefit: The folder name for linefit results from Hipe
                               (created by Pierre, assuming his syntax). The 
                               folder is located in $path_pacs$/$star_name$/.
                               If no folder is given, no linefits are 
                               processed.
                               
                               (default: '')
        @type path_linefit: string
        
        '''
        
        super(Pacs,self).__init__(star_name=star_name,code='GASTRoNOoM',\
                                  path_linefit=path_linefit,path=path,\
                                  oversampling=oversampling,\
                                  instrument_name='PACS',intrinsic=intrinsic)
        self.data_wave_list = []
        self.data_flux_list = []
        self.data_ordernames = []
        self.data_orders = []
        self.data_delta_list = []
        self.redo_convolution = redo_convolution

        if self.path <> None:
            #-- Convenience path
            cc.path.gout = os.path.join(cc.path.gastronoom,self.path)
            #-- Check the path for the PACS database if a model folder is known
            db_path = os.path.join(cc.path.gout,'stars',self.star_name)
            if os.path.isdir(db_path):
                self.db_path = os.path.join(db_path,'GASTRoNOoM_pacs_models.db')
                self.db = Database.Database(self.db_path)
                DataIO.testFolderExistence(os.path.join(cc.path.gout,'stars',\
                                           self.star_name,'PACS_results'))
            else:
                self.db = None
                self.db_path = None
        
        self.sphinx_prep_done = 0
        self.readLineFit()


    def readLineFit(self):
        
        '''
        Read the data from the line fit procedure done with Pierres Hipe 
        script.
        
        Assumes structure and syntax as given in the example file
        /home/robinl/Data/PACS/v1362aql/lineFit/lineFitResults
        
        The line fit results are saved in self.linefit as a np.recarray.
        
        The columns include (with unit indicated): 
        groupid
        band
        wave_in (micron), 
        wave_fit (micron), 
        line_flux (W/m2),
        line_flux_err (W/m2), 
        line_flux_rel (%), 
        continuum (W/m2 m), 
        line_peak (W/m2 m), 
        fwhm_fit (micron),
        fwhm_pacs (micron),
        fwhm_rel
        
        '''
        
        kwargs = dict([('start_from_keyword','GROUPID'),\
                       ('start_row',1)])
        dd = super(Pacs,self).readLineFit(**kwargs)
        if dd is None: return
        
        #-- If no % symbol in 6th last column, new way of giving int ints:
        #   get rid of last 2 columns which do not contain relevant information
        #   as well as the wave ratio column
        if type(dd[-6][0]) is not types.StringType:
            del dd[-1]
            del dd[-1]
            del dd[-10]
        #-- Convert relative error in % to ratio
        dd[-6] = [float(val.strip('%'))/100. for val in dd[-6]]
        del dd[-8]
        #bands = []
        #obs_bands = set([(band,obsid) 
                         #for band,obsid in zip(dd[-11],dd[-13])
                         #if 'R1' not in band])
        #for ival,val in enumerate(dd[-11]):
            #if 'R1' in val:
                #thisid = dd[-13][ival]
                #b2b = ('B2B',thisid)
                #b2a = ('B2A',thisid)
            #else:
                #bands.append(val)
        del dd[1:-11]
        names = ('groupid','band','wave_in','wave_fit','line_flux',\
                 'line_flux_err','line_flux_rel','continuum','line_peak',\
                 'fwhm_fit','fwhm_pacs','fwhm_rel')
        self.linefit = np.recarray(shape=(len(dd[-1]),),\
                                   dtype=zip(names,[int]+['|S3']+[float]*10))
        for n,d in zip(names,dd):
            self.linefit[n] = d

                
                
    def addStarPars(self,star_grid):
        
        '''
        Set parameters taken from the PACS database into empty parameter sets.
        
        Typically used in conjunction with making Star() templates through 
        Star.makeStars().
        
        @param star_grid: The parameter sets that will updated in situ
        @type star_grid: list[Star()]

        '''

        for star in star_grid:    
            model = star['LAST_PACS_MODEL']
            pacs_model = self.db[model].copy()
            star.update({'LAST_GASTRONOOM_MODEL': pacs_model['cooling_id'],\
                         'MOLECULE': list(set(\
                                [(db_entry[2]['TRANSITION'].split()[0],\
                                  db_entry[1]) 
                                 for db_entry in pacs_model['trans_list']])),\
                         'TRANSITION': \
                                [db_entry[2]['TRANSITION'].split() + \
                                 [db_entry[2]['N_QUAD']] 
                                 for db_entry in pacs_model['trans_list']]})
            [trans.setModelId(db_entry[0]) 
             for trans,db_entry in zip(star['GAS_LINES'],\
                                       pacs_model['trans_list'])]



    def setData(self,data_filenames=[],searchstring=''):
        
        '''
        Read data from given filenames and remember them.
        
        If no filenames given, an auto search is performed.
        
        The data are saved as wavelength and flux lists in the object.
        
        @keyword data_filenames: The data filenames. If empty, auto-search is 
                                 done
                                 
                                 (default: [])
        @type data_filenames: list[string]
        @keyword searchstring: the searchstring conditional for the auto-search
        
                               (default: '')
        @type searchstring: string
        
        '''
        
        super(Pacs,self).setData(data_filenames, searchstring)
        if self.data_filenames:
            self.data_orders = [int(o[1]) for o in self.data_ordernames]
            


    def prepareSphinx(self,star_grid,redo_sphinx_prep=0):
        
        '''
        Prepare Sphinx PACS models by checking if the Star() instance already 
        has a PACS id associated with it, and if not calling convolveSphinx.
        
        @param star_grid: The parameter sets
        @type star_grid: list[Star()]
        @keyword redo_sphinx_prep: redo the sphinx prep, regardless of whether 
                                   this Pacs instance already did it once (for 
                                   instance in case the star_grid changed)
                                   
                                   (default: 0)
        @type redo_sphinx_prep: bool
        
        '''
        
        if not self.sphinx_prep_done or redo_sphinx_prep:
            print '** Loading from database, or convolving with ' + \
                  'Gaussian and rebinning to data wavelength grid.'  
            for i,star in enumerate(star_grid):
                print '* Sphinx model %i out of %i.' %(i+1,len(star_grid))
                if not star['LAST_GASTRONOOM_MODEL']: 
                    print '* No cooling model found.'
                else:
                    self.__convolveSphinx(star=star)
                    if star['LAST_PACS_MODEL']:
                        print '* %s is done!'%star['LAST_PACS_MODEL']
            self.sphinx_prep_done = 1
            self.db.sync()
                      


    def getSphinxConvolution(self,star,fn):
        
        '''
        Read the sphinx convolution and return if it has already been done. 
        
        Returns None if the convolution is not available. 
        
        @param star: The Star() object
        @type star: Star()
        @param fn: The filename of the dataset (band) for which the convolution
                   is to be returned.
        @type fn: str
        
        @return: The sphinx convolution result. (wavelength, flux)
        @rtype: array
        
        '''
        
        this_id = star['LAST_PACS_MODEL']
        if not this_id:
            return ([],[])
        fn = os.path.split(fn)[1]
        sphinx_file = os.path.join(cc.path.gout,'stars',self.star_name,\
                                  'PACS_results',this_id,'%s_%s'%('sphinx',fn))
        return DataIO.readCols(sphinx_file)



    def __convolveSphinx(self,star):
        
        '''
        Check if sphinx output has already been convolved (pacs db) and do so 
        if not.
        
        @param star: The parameter set
        @type star: Star()
        
        '''     
        
        #- check for which filenames the convolution has already been done
        finished_conv_filenames = self.checkStarDb(star)
        finished_conv_filenames = [this_f 
                                   for this_f in finished_conv_filenames 
                                   if this_f in [os.path.split(f)[1] 
                                                 for f in self.data_filenames]]
        
        #- if after db check pacs id is undefined, then there is no pacs model,
        #- and the convolution will be done anyway
        if self.redo_convolution and star['LAST_PACS_MODEL']:    
            for filename in finished_conv_filenames:
                ori = os.path.join(cc.path.gout,'stars',self.star_name,\
                                   'PACS_results',star['LAST_PACS_MODEL'],\
                                   '_'.join(['sphinx',filename]))
                backup = os.path.join(cc.path.gout,'stars',self.star_name,\
                                      'PACS_results',star['LAST_PACS_MODEL'],\
                                      '_'.join(['backup','sphinx',filename]))
                subprocess.call(['mv %s %s'%(ori,backup)],shell=True)
            self.db[star['LAST_PACS_MODEL']]['filenames'] = []
            finished_conv_filenames = []
        filenames_to_do = [this_f 
                           for this_f in [os.path.split(f)[1] 
                                          for f in self.data_filenames]
                           if this_f not in finished_conv_filenames]     
        
        #-Get sphinx model output and merge, for all star models in star_grid
        if filenames_to_do:        
            #- if list is not empty, some filenames still need a convolution    
            print '* Reading Sphinx model and merging.'
            merged = star['LAST_GASTRONOOM_MODEL'] \
                        and self.mergeSphinx(star) \
                        or [[],[]]
            sphinx_wave = merged[0]
            sphinx_flux = merged[1]
            if not sphinx_wave: 
                print '* No Sphinx data found.'
                return
  
            #- convolve the model fluxes with a gaussian at central wavelength 
            #- from data_wave_list for every star, and appropriate sigma
            print '* Convolving Sphinx model, after correction for v_lsr.'
            if not self.data_delta_list:
                self.setDataResolution()
            for filename in filenames_to_do:
                i_file = [os.path.split(f)[1] 
                          for f in self.data_filenames].index(filename)
                if not star['LAST_PACS_MODEL']:
                    star['LAST_PACS_MODEL'] = \
                        'pacs_%.4i-%.2i-%.2ih%.2i-%.2i-%.2i'\
                        %(gmtime()[0],gmtime()[1],gmtime()[2],\
                          gmtime()[3],gmtime()[4],gmtime()[5])
                    self.db[star['LAST_PACS_MODEL']] = \
                        dict([('filenames',[]),\
                              ('trans_list',star.getTransList(dtype='PACS')),\
                              ('cooling_id',star['LAST_GASTRONOOM_MODEL'])])
                    DataIO.testFolderExistence(\
                        os.path.join(cc.path.gout,'stars',self.star_name,\
                                     'PACS_results',star['LAST_PACS_MODEL']))
                #-- Correct for the v_lsr of the central source
                sphinx_wave_corr = array(sphinx_wave)*(1./(1-self.vlsr/self.c))
                sph_conv = Data.doConvolution(\
                                        x_in=sphinx_wave_corr,\
                                        y_in=sphinx_flux,\
                                        x_out=self.data_wave_list[i_file],\
                                        widths=self.data_delta_list[i_file],\
                                        oversampling=self.oversampling)
                sph_fn = os.path.join(cc.path.gout,'stars',self.star_name,\
                                      'PACS_results',star['LAST_PACS_MODEL'],\
                                      '_'.join(['sphinx',filename])) 
                DataIO.writeCols(filename=sph_fn,\
                                 cols=[self.data_wave_list[i_file],sph_conv])
                self.db[star['LAST_PACS_MODEL']]['filenames'].append(filename)
                self.db.addChangedKey(star['LAST_PACS_MODEL'])        
        
        #- Idea:
        #- 1) grab model flux and wavelength, data flux and wavelengthm, 
        #-    instrument intrinsic resolution and wavelength
        #- 2) scan model fluxes with gaussian window of variable sigma, 
        #-    which is equal to intrinsic resolution at gaussian center0
        #-    and gaussian center at wavelengths in the data arrays, 
        #-    which may or may not be spaced by the intrinsic resolution;
        #-    take sigma = intrinsic_resolution/delta(lambda) * len(lambda)
        #-    where lambda is [abs(l_model-l_instr)<intrinsic_resolution]
        #- 3) the list found has the same wavelengths as the data input and can
        #-    be directly compared
        


    def getPacsResolution(self):
        
        '''
        Get the Pacs resolution from cc aux file for all orders.
        
        @return: The resolution as a function of wavelength for the 3 orders
        @rtype: (list[list],list[list])
        
        '''
        
        reso_wave_list = [DataIO.getInputData(path=cc.path.aux,\
                                              filename='Pacs_Resolution.dat',\
                                              keyword='WAVE='+str(order))
                          for order in range(1,4)]
        reso_delta_list = [DataIO.getInputData(path=cc.path.aux,\
                                               filename='Pacs_Resolution.dat',\
                                               keyword='ORDER='+str(order))
                           for order in range(1,4)]
        return reso_wave_list,reso_delta_list



    def setDataResolution(self):
        
        '''
        Read and remember the PACS resolution, based on the wavelength grid of 
        the data. 
        
        Returns a list of deltas for every filename, within a list.
        
        '''
        
        print '** Reading PACS native resolution.'
        reso_wave_list,reso_delta_list = self.getPacsResolution()
        
        #- Make interpolators for every order
        print '** Interpolating PACS native resolution.'
        interpolator_list = [interpolate.interp1d(x,y) 
                             for x,y in zip(reso_wave_list,reso_delta_list)]
        
        #- Interpolate the PACS resolution for appropriate order to get sigma's
        #- for data wavelengths
        self.data_delta_list = [interpolator_list[order-1](x_new) 
                                for x_new,order in zip(self.data_wave_list,\
                                                       self.data_orders)]        



    def checkStarDb(self,star):
        
        '''
        Check if for this star the convolution has already been done.
        
        @param star: The parameter set
        @type star: Star()
        
        @return: The filenames of PACS convolutions found for this parameter
                 set. Can be multiple filenames if different bands associated 
                 with this set.
        @rtype: list[string]
        
        '''
        
        if star['LAST_PACS_MODEL']:
            print star['LAST_PACS_MODEL']
            return self.db[star['LAST_PACS_MODEL']]['filenames']
        else:
            #selection = [(k,v) for k,v in self.db.items() 
            #             if v['cooling_id'] == star['LAST_GASTRONOOM_MODEL'] \
            #                and self.compareTransLists(starlist=star.getTransList(),dblist=v['trans_list'])]
            selection = [(k,v) 
                         for k,v in self.db.items() 
                         if v['cooling_id'] == star['LAST_GASTRONOOM_MODEL'] \
                            and star.getTransList(dtype='PACS') == v['trans_list']]
            if len(selection) == 1:
                star['LAST_PACS_MODEL'] = selection[0][0]
                return selection[0][1]['filenames']
            elif len(selection) > 1:
                raise IOError('Multiple PACS ids found for a single star. ' + \
                              'Something very fishy is going on... Contact ' +\
                              'Robin with this error message!')
            else:
                return []
                