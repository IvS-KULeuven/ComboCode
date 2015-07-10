# -*- coding: utf-8 -*-

"""
Running MCMax and managing output from MCMax.

Author: R. Lombaert

"""

import os
import subprocess
from glob import glob

import cc.path
from cc.tools.io import DataIO, Database
from cc.modeling.ModelingSession import ModelingSession



def readModelSpectrum(dpath,rt_sed=1,fn_spec='spectrum45.0.dat'):
     
    '''
    Read the model output spectrum.
     
    If no ray-tracing is requested or no ray-tracing output is found, the 
    average of the MC spectra is taken.
     
    @param dpath: folder that contains the MCMax outputfiles
    @type dpath: string
    
    @keyword rt_sed: If a ray-traced spectrum is requested
     
                     (default: 1)
    @type rt_sed: bool
    @keyword fn_spec: The filename of the ray-traced spectrum. Typically this 
                      is the default name, but can be different depending on 
                      the ray-tracing angle that is used. 
                      Not used if MCSpec are used.
                      
                      (default: spectrum45.0.dat)
    @type fn_spec: str
    
    @return: The wavelength and flux grids (micron,Jy)
    @rtype: (array,array)
     
    '''
     
    rt_sed = int(rt_sed)
    try:    
        if rt_sed:  
            dfile = os.path.join(dpath,fn_spec)
            this_data = DataIO.readCols(dfile)
            #- if the lists are not empty
            if list(this_data[0]) and list(this_data[1]):    
                w = this_data[0]
                f = this_data[1]
            else: raise IOError
        else: raise IOError                                
    except IOError:
        print 'No spectrum was found or ray-tracing is off for ' + \
              'this model. Taking average of theta-grid MCSpectra.'
        dfiles = glob(os.path.join(dpath,'MCSpec*.dat'))
        w = DataIO.readCols(filename=dfiles[0])[0]
        mcy_list = [DataIO.readCols(f)[1] for f in dfiles]
        f = sum(mcy_list)/len(mcy_list)
    return (w,f)



def readVisibilities(dpath,fn_vis='visibility01.0.dat'):
    
    '''
    Read the model output visibilities, either as function of wavelength or
    baseline. 
     
    @param dpath: folder that contains the MCMax outputfiles
    @type dpath: string
    
    @keyword fn_spec: The filename of the ray-traced visibilities. Typically 
                      this is the default name, but can be different depending 
                      on the inclination (or baseline) that is used. 
                      
                      (default: visibility01.0.dat)
    @type fn_spec: str
    
    @return: A dictionary containing either wavelength or baseline, the flux, 
             and the visibilities for either given baselines or wavelengths
    @rtype: dict
     
    '''
    
    #-- Read file and 
    dfile = os.path.join(dpath,fn_vis)
    if not os.path.isfile(dfile):
        return dict()
    cols, comments = DataIO.readCols(dfile,return_comments=1)
    comments = [comment for comment in comments if comment]
    
    if 'visibility' in fn_vis: 
        xtype = 'wavelength'
        seltype = 'baseline'
    elif 'basevis' in fn_vis: 
        xtype = 'baseline'
        seltype = 'wavelength'
    model = dict()
    model[xtype] = cols[0]
    model['flux'] = cols[1]
    model[seltype] = dict()
    for i,comment in enumerate(comments[2:]):
        val = float(comment.partition(seltype)[2].partition(',')[0])
        model[seltype][val] = cols[2+i]
    return model



def rayTraceSpectrum(model_id,path_mcmax='runTestDec09',inputfilename='',\
                     redo_rt=0):
    
    '''
    Do the ray-tracing of the spectrum according to cc.path.mobs/Spec.out, but 
    only if spectrum45.0.dat does not exist yet.
    
    @param model_id: the model_id of the requested model
    @type model_id: string
    
    @keyword path_mcmax: modeling folder in MCMax home
    @type path_mcmax: string
    @keyword inputfilename: the inputfilename of the model. if '': filename is 
                            inputMCMax_model_YYYY-MM-DDhHH-mm-ss.
                            
                            (default: '')
    @type inputfilename: string          
    @keyword redo_rt: redo the ray tracing of the spectrum regardless of the 
                      spectrum already existing or not
                      
                      (default: 0)
    @type redo_rt: bool
    
    '''
    
    redo_rt = int(redo_rt)
    spectrum_file = os.path.join(cc.path.mcmax,path_mcmax,'models',model_id,\
                                 'spectrum45.0.dat')
    if not os.path.isfile(spectrum_file) or redo_rt:
        print '** Ray-tracing spectrum now...'
        if not inputfilename:
            inputfilename=os.path.join(cc.path.mcmax,path_mcmax,'models',\
                                       'inputMCMax_%s.dat'%model_id)
        output_folder = os.path.join(cc.path.mcmax,path_mcmax,'models',\
                                     model_id)
        spec_file = os.path.join(cc.path.mobs,'Spec.out')
        subprocess.call([' '.join(['MCMax ' + inputfilename,'0','-o',\
                                   output_folder,spec_file])],shell=True)
    else:
        print '** Spectrum ray-tracing is already finished.'
        


def rayTraceImage(model_id,path_mcmax='runTestDec09',inputfilename='',\
                  remove_source=0,output_folder=''):
    
    '''
    Do the ray-tracing of images, according to cc.path.mobs/Image.out.
    
    @param model_id: the model_id of the requested model
    @type model_id: string
        
    @keyword path_mcmax: modeling folder in MCMax home
    @type path_mcmax: string
    @keyword inputfilename: the inputfilename of the model. if '': filename is 
                            inputMCMax_model_YYYY-MM-DDhHH-mm-ss.
                            
                            (default: '')
    @type inputfilename: string                        
    @keyword remove_source: remove the central source from the image
    
                            (default: 0)
    @type remove_source: bool
    @keyword output_folder: The location of the output folder. By default, set 
                            at the model output folder.
    
                            (default: '')
    @type output_folder: str
    
    '''
    
    remove_source = int(remove_source)
    print '** Ray-tracing images now...'
    if not inputfilename:
        inputfilename=os.path.join(cc.path.mcmax,path_mcmax,'models',\
                                   'inputMCMax_%s.dat'%model_id)
    model_folder = os.path.join(cc.path.mcmax,path_mcmax,'models',model_id)
    image_file = os.path.join(cc.path.mobs,'Image.out')
    if remove_source:
        subprocess.call([' '.join(['MCMax ' + inputfilename,'0',\
                                   '-s tracestar=.false.','-o',model_folder,\
                                   image_file])],shell=True)
    else:
        subprocess.call([' '.join(['MCMax ' + inputfilename,'0','-o',\
                                   model_folder,image_file])],shell=True)
    
    if not output_folder:
        output_folder = model_folder
    else: 
        subprocess.call(['mv %s %s'%(os.path.join(model_folder,'Image*'),\
                                     os.path.join(output_folder,'.'))],\
                        shell=True)
    print '** Your images can be found at:'
    print output_folder
                                


def rayTraceVisibilities(model_id,path_mcmax='runTestDec09',inputfilename=''):
    
    '''
    Do the ray-tracing of visibilities, according to 
    cc.path.mobs/Visibilities.out.
    
    @param model_id: the model_id of the requested model
    @type model_id: string
        
    @keyword path_mcmax: modeling folder in MCMax home
    @type path_mcmax: string
    @keyword inputfilename: the inputfilename of the model. if '': filename is 
                            inputMCMax_model_YYYY-MM-DDhHH-mm-ss.
                            
                            (default: '')
    @type inputfilename: string                        
    
    '''
    
    print '** Ray-tracing visibilities now...'
    if not inputfilename:
        inputfilename=os.path.join(cc.path.mcmax,path_mcmax,'models',\
                                   'inputMCMax_%s.dat'%model_id)
    output_folder = os.path.join(cc.path.mcmax,path_mcmax,'models',model_id)
    visibilities_file = os.path.join(cc.path.mobs,'Visibilities.out')
    subprocess.call([' '.join(['MCMax ' + inputfilename,'0','-o',\
                               output_folder,visibilities_file])],shell=True)
    print '** Your visibilities can be found at:'
    print output_folder
                               
                               

class MCMax(ModelingSession):
    
    """ 
    Class that includes all methods required for creating an MCMax model. 
    
    """
    
    def __init__(self,path_mcmax='runTest',replace_db_entry=0,db=None,\
                 new_entries=[]):
        
        """ 
        Initializing an instance of ModelingSession.
        
        @keyword db: the MCMax database
        
                          (default: None)
        @type db: Database()
        @keyword replace_db_entry: replace an entry in the MCMax database with 
                                   a newly calculated model with a new model id 
                                   (for instance if some general data not 
                                   included in the inputfiles is changed)
                                   
                                   (default: 0)
        @type replace_db_entry: bool
        @keyword path_mcmax: modeling folder in MCMax home
        
                             (default: 'runTest')
        @type path_mcmax: string
        @keyword new_entries: The new model_ids when replace_db_entry is 1
                                   of other models in the grid. These are not 
                                   replaced!
                                   
                                   (default: [])
        @type new_entries: list[str]     
        
        """
        
        super(MCMax, self).__init__(code='MCMax',path=path_mcmax,\
                                    replace_db_entry=replace_db_entry,\
                                    new_entries=new_entries)
        #-- Convenience path
        cc.path.mout = os.path.join(cc.path.mcmax,self.path)
        DataIO.testFolderExistence(os.path.join(cc.path.mout,\
                                                'data_for_gastronoom'))
        self.db = db
        self.mcmax_done = False
        
        #- Read standard input file with all parameters that should be included
        #- as well as some dust specific information
        inputfilename = os.path.join(cc.path.aux,'inputMCMax.dat')
        self.standard_inputfile = DataIO.readDict(inputfilename,\
                                                  convert_floats=1,\
                                                  convert_ints=1,\
                                                  comment_chars=['#','*'])
                

    def rayTrace(self,star):

        '''
        Ray trace the spectrum and images, if requested and not yet finished.

        @param star: The parameter set for this session
        @type star: Star()
        
        '''

        if self.model_id and int(star['RT_SED']):
            rayTraceSpectrum(model_id=self.model_id,path_mcmax=self.path,\
                             redo_rt=star['REDO_OBS'])
        if self.model_id and int(star['IMAGE']):
            rayTraceImage(model_id=self.model_id,path_mcmax=self.path,\
                          remove_source=star['IMAGE_NOSOURCE'])
        if self.model_id and int(star['VISIBILITIES']):
            rayTraceVisibilities(model_id=self.model_id,path_mcmax=self.path)


            
    def setCommandKey(self,comm_key,star,star_key=None,alternative=None):
        
        '''
        Try setting a key in the command_list from a star instance. 
        
        If the key is unknown, it is left open and will be filled in from the 
        standard gastronoom inputfile.
        
        @param comm_key: the name of the keyword in the command list
        @type comm_key: string
        @param star: Parameter set for this session
        @type star: Star()
        
        @keyword star_key: the name of the keyword in Star() (minus '_%s'
                           %key_type (DUST or GAS), which is added as well in a 
                           second attempt if the first without the addition is 
                           not found). If None, it's equal to comm_key
                           
                           (default: None)
        @type star_key: string
        @keyword alternative: a default value passed from the standard 
                              inputfile that is used if the keyword or the 
                              keyword + '_%s'%key_type is not found in Star().
                              
                              (default: None)
        @type alternative: string
        
        @return: Success? 
        @rtype: bool
        
        '''    
        
        keyword_int_list = ['nrad','ntheta','nspan','nlev','ntspan','ntlev',\
                            'nlam','nzlam','tmaxiter','nbw']
        if comm_key.lower() in keyword_int_list: make_int = 1
        else: make_int = 0     
        return super(MCMax, self).setCommandKey(comm_key,star,'DUST',\
                                                star_key,alternative,\
                                                make_int)



    def compareCommandLists(self,this_list,modellist):
        
        """
        Comparing a command_list with a database entry.
        
        @param this_list: parameters in this modeling session
        @type this_list: dict
        @param modellist: parameters from database model
        @type modellist: dict
        
        @return: Comparison between the two parameter sets
        @rtype: bool
        
        """
        
        #- Run the modelingsession entry of this method, then if True check
        #- the dust species entry in the command list
        keybools = super(MCMax, self).compareCommandLists(this_list=this_list,\
                                                          modellist=modellist,\
                                                          code='mcmax')
        if keybools:
            #- Check the dust species entries 
            this_dust = this_list['dust_species']
            model_dust = modellist['dust_species']
            if sorted(this_dust.keys()) != sorted(model_dust.keys()):
                return False
            else:
                dustbools = []
                for species in this_dust.keys():
                    dustbools.append(super(MCMax,self)\
                        .compareCommandLists(this_list=this_dust[species],\
                                             modellist=model_dust[species],\
                                             code='mcmax'))
                if False in dustbools:
                    return False
                else:
                    return True
        else:
            return False
        
        
        
    def checkDatabase(self):

        """
        Checking cooling database.
        
        @return: The presence of the MCMax model in the database
        @rtype: bool
        
        """
        
        for model_id,cool_dict in sorted(self.db.items()):
            model_bool = self.compareCommandLists(self.command_list.copy(),\
                                                  cool_dict)
            if model_bool:
                if self.replace_db_entry \
                        and model_id not in self.new_entries: 
                    print 'Replacing MCMax database entry for old ID %s.'\
                          %model_id
                    del self.db[model_id]
                    return False
                else:
                    print 'MCMax model has been calculated ' + \
                          'before with ID %s.'%model_id
                    self.model_id = model_id
                    return True
        print 'No match found in MCMax database. Calculating new model.'
        return False
        
        
            
    def doMCMax(self,star):
        
        """
        Running MCMax.
        
        @param star: The parameter set for this session
        @type star: Star()
        
        """

        print '***********************************'                                       
        #- Create the input dictionary for this MCMax run
        print '** Making input file for MCMax'
        #-- Add the previous model_id to the list of new entries, so it does 
        #   not get deleted if replace_db_entry == 1. 
        if self.model_id: 
            self.new_entries.append(self.model_id)
        self.model_id = ''
        self.command_list = dict()
        self.command_list['photon_count'] = star['PHOTON_COUNT']
        if star['STARFILE']:
            self.command_list['startype'] = "'FILE'"
            self.command_list['Lstar'] = star['L_STAR']
            self.command_list['Rstar'] = star['R_STAR']
            self.command_list['starfile'] = "'%s'"%star['STARFILE']
        else:
            self.command_list['Tstar'] = float(star['T_STAR'])
            self.command_list['Rstar'] = float(star['R_STAR'])            
        self.command_list['tcontact'] = star['T_CONTACT'] \
                                            and '.true.' \
                                            or '.false.'
        for boolkey in ['FLD','iter','randomwalk','tdesiter','storescatt',\
                        'multiwav']:
            try:
                self.command_list[boolkey] = int(star[boolkey.upper()]) \
                                                and '.true.' \
                                                or '.false.'
            except KeyError:
                self.setCommandKey(boolkey,star,star_key=boolkey.upper(),\
                                   alternative=self.standard_inputfile[boolkey])
            
        if not int(star['TDESITER']):
            self.command_list['Rin'] = float(star['R_INNER_DUST'])\
                                        *float(star['R_STAR'])\
                                        *star.Rsun/star.au
        self.command_list['Rout'] = float(star['R_OUTER_DUST'])\
                                        *float(star['R_STAR'])\
                                        *star.Rsun/star.au
        
        self.command_list['denstype'] = "'%s'"%star['DENSTYPE']
        if star['DENSTYPE'] == 'MASSLOSS':
            self.command_list['Mdot'] = float(star['MDOT_DUST'])*100.
            self.command_list['vexp'] = float(star['V_EXP_DUST'])
        elif star['DENSTYPE'] == 'SHELLFILE' or star['DENSTYPE'] == 'FILE':
            self.command_list['densfile'] = "'%s'"%star['DENSFILE']
        elif star['DENSTYPE'] == 'POW':
            self.command_list['denspow'] = star['DENSPOW']
            self.command_list['mdust'] = star['M_DUST']
        if int(star['MRN_DUST']):
            self.command_list['mrn'] = '.true.'
            self.command_list['mrn_index'] = star['MRN_INDEX']
            self.command_list['mrn_ngrains'] = star['MRN_NGRAINS']
            self.command_list['mrn_rmax'] = star['MRN_RMAX']
            self.command_list['mrn_rmin'] = star['MRN_RMIN']
        if int(star['SCSET']):
            self.command_list['scset'] = '.true.'
            self.command_list['scseteq'] = int(star['SCSETEQ']) \
                                                and '.true' or '.false.'
            self.command_list['alphaturb'] = star['ALPHATURB']
        self.setCommandKey('Mstar',star,star_key='M_STAR',\
                           alternative=self.standard_inputfile['Mstar'])       
        add_keys = [k  
                    for k in self.standard_inputfile.keys() 
                    if not self.command_list.has_key(k)]
        [self.setCommandKey(k,star,star_key=k.upper(),\
                            alternative=self.standard_inputfile[k])
         for k in add_keys]
        
        #- Following keywords are only added if corresponding abundance is
        #- present, ie not zero.
        #- Note that no indices are associated with each species yet. This is 
        #- done in the creation of the inputfile for MCMax. This way, confusion
        #- is avoided when comparing to older models, in the case the order in
        #- Dust.dat is changed. However, this does mean some species specific
        #- outputfiles may not match anymore when e.g. plotting dust opacities
        dust_dict = dict()
        for species in star.getDustList():
            species_dict = dict()
            if star['TDESITER']:
                species_dict['TdesA'] = star['T_DESA_' + species]
                species_dict['TdesB'] = star['T_DESB_' + species]
            if star.has_key('R_MIN_%s'%species) and star['R_MIN_%s'%species]:
                species_dict['minrad'] = star['R_MIN_%s'%species]\
                                           *star['R_STAR']*star.Rsun/star.au
            #- R_MAX is always created by Star(), even if not requested.
            #- Will be empty string if not available; no maxrad is given
            if star['R_MAX_%s'%species]:
                species_dict['maxrad'] = star['R_MAX_%s'%species]\
                                           *star['R_STAR']*star.Rsun/star.au
            if int(star['MRN_DUST']) and star.has_key(['RGRAIN_%s'%species]):
                species_dict['rgrain'] = star['RGRAIN_%s'%species]
            else:
                species_dict['abun'] = star['A_%s'%species]
            if not os.path.split(star.dust[species]['fn'])[0]:
                print('WARNING! %s has an old opacity file. Should replace for reproducibility.'%species)
            dust_dict[star.dust[species]['fn']] = species_dict
        self.command_list['dust_species'] = dust_dict
        print '** DONE!'
        print '***********************************'
        
        #-- Check the MCMax database if the model was calculated before
        modelbool = self.checkDatabase()
                
        #-- if no match found in database, calculate new model with new model id 
        #-- if the calculation did not fail, add entry to database for new model
        if not modelbool:
            self.model_id = self.makeNewId()
            input_dict = self.command_list.copy()
            del input_dict['photon_count']
            del input_dict['dust_species']
            #-- dust_list in star is already sorted. rgrains species first, 
            #   then the rest, according to the order of appearance in Dust.dat
            for index,species in enumerate(star.getDustList()):
                speciesfile = star.dust[species]['fn']
                speciesdict = self.command_list['dust_species'][speciesfile]
                for k,v in speciesdict.items():
                    input_dict['%s%.2i'%(k,index+1)] = v
                #-- If speciesfile is .topac, they are T-dependent opacities
                #   and should always be given as topac##. The file then points
                #   to the .particle files of the T-dependent opacities.
                if speciesfile.find('.topac') != -1:
                    ftype = 'topac'
                #-- When full scattering is requested, always use .particle 
                #   files if they are available. If not, use whatever is in 
                #   Dust.dat, but then the species will not be properly 
                #   included for full scattering (requires scattering matrix)
                #   It is OK to have .opac files in Dust.dat, as long as 
                #   .particle files exist in the same location
                elif star['SCATTYPE'] == 'FULL':
                    partfile = os.path.splitext(speciesfile)[0] + '.particle'
                    if os.isfile(os.path.join(cc.path.mopac,partfile)):
                        ftype = 'part'
                        speciesfile = partfile
                    else:
                        ftype = 'opac'
                #-- If not full scattering, opacity files are fine. So, use 
                #   whatever is in Dust.dat. Dust.dat should preferentially 
                #   include .opac (or .opacity) files, but can be .particle 
                #   files if opacity files are not available, in which case
                #   ftype should still be 'part'.
                else:
                    if speciesfile.find('.particle') != -1: ftype = 'part'
                    else: ftype = 'opac'
                #-- Add the opacities home folder (not saved in db)
                input_dict['%s%.2i'%(ftype,index+1)] = "'%s'"\
                            %(os.path.join(cc.path.mopac,speciesfile))       
            input_filename = os.path.join(cc.path.mout,'models',\
                                          'inputMCMax_%s.dat'%self.model_id)
            output_folder = os.path.join(cc.path.mout,'models',self.model_id)
            input_lines = ["%s=%s"%(k,str(v)) 
                           for k,v in sorted(input_dict.items())]
            DataIO.writeFile(filename=input_filename,input_lines=input_lines)
            subprocess.call(' '.join(['MCMax',input_filename,\
                                      str(self.command_list['photon_count']),\
                                      '-o',output_folder]),shell=True)
            self.mcmax_done = True
            testf1 = os.path.join(output_folder,'denstemp.dat')
            testf2 = os.path.join(output_folder,'kappas.dat')
            if os.path.exists(testf1) and os.path.exists(testf2) and \
                    os.path.isfile(testf1) and os.path.isfile(testf2):
                self.db[self.model_id] = self.command_list
                self.db.sync()
            else:
                print '** Model calculation failed. No entry is added to ' + \
                      'the database and LAST_MCMAX_MODEL in STAR dictionary '+\
                      'is not updated.'
                self.model_id = ''
                
        #- add/change 'LAST_MCMAX_MODEL' entry
        if self.model_id:
            star['LAST_MCMAX_MODEL'] = self.model_id
        #- Note that the model manager now adds/changes MUTABLE input keys, 
        #- which MAY be overwritten by the input file inputComboCode.dat
        print '***********************************'
        
