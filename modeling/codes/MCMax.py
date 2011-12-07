# -*- coding: utf-8 -*-

"""
Running MCMax and managing output from MCMax.

Author: R. Lombaert

"""

import os
from time import gmtime
import cPickle
import subprocess
from glob import glob

from cc.tools.io import DataIO
from cc.modeling.ModelingSession import ModelingSession
from cc.tools.io import Database



def readModelSpectrum(path_mcmax,model_id,rt_sed=1):
     
    '''
    Read the model output spectrum.
     
    If no ray-tracing is requested or no ray-tracing output is found, the 
    average of the MC spectra is taken.
     
    @param path_mcmax: modeling folder in MCMax home
    @type path_mcmax: string
    @param model_id: Model id of the spectrum to be read
    @type model_id: string
    
    @keyword rt_sed: If a ray-traced spectrum is requested
     
                     (default: 1)
    @type rt_sed: bool
    
    @return: The wavelength and flux grids
    @rtype: (array,array)
     
    '''
     
    rt_sed = int(rt_sed)
    try:    
        if rt_sed:  
            dfile = os.path.join(os.path.expanduser('~'),'MCMax',path_mcmax,\
                                 'models',model_id,'spectrum45.0.dat')
            this_data = DataIO.readCols(dfile)
            #- if the lists are not empty
            if list(this_data[0]) and list(this_data[1]):    
                w = this_data[0]
                f = this_data[1]
            else: raise IOError
        else: raise IOError                                
    except IOError:
        print 'No spectrum was found or ray-tracing is off for ' + \
              '%s. Taking average of theta-grid MCSpectra.'%model_id
        dfiles = glob(os.path.join(os.path.expanduser('~'),'MCMax',path_mcmax,\
                      'models',model_id,'MCSpec*.dat'))
        w = DataIO.readCols(filename=dfiles[0])[0]
        mcy_list = [DataIO.readCols(f)[1] for f in dfiles]
        f = sum(mcy_list)/len(mcy_list)
    return (w,f)



def rayTraceSpectrum(model_id,path_mcmax='runTestDec09',inputfilename='',\
                     redo_rt=0):
    
    '''
    Do the ray-tracing of the spectrum according to ~/MCMax/src/Spec.out, 
    but only if spectrum45.0.dat does not exist yet.
    
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
    spectrum_file = os.path.join(os.path.expanduser("~"),'MCMax',path_mcmax,\
                                 'models',model_id,'spectrum45.0.dat')
    if not os.path.isfile(spectrum_file) or redo_rt:
        print '** Ray-tracing spectrum now...'
        if not inputfilename:
            inputfilename=os.path.join(os.path.expanduser('~'),'MCMax',\
                                       path_mcmax,'models',\
                                       'inputMCMax_%s.dat'%model_id)
        output_folder = os.path.join(os.path.expanduser("~"),'MCMax',\
                                     path_mcmax,'models',model_id)
        spec_file = os.path.join(os.path.expanduser('~'),'MCMax','src',\
                                 'Spec.out')
        subprocess.call([' '.join(['MCMax ' + inputfilename,'0','-o',\
                                   output_folder,spec_file])],shell=True)
    else:
        print '** Spectrum ray-tracing is already finished.'
        


def rayTraceImage(model_id,path_mcmax='runTestDec09',inputfilename='',\
                  remove_source=0):
    
    '''
    Do the ray-tracing of images, according to ~/MCMax/src/Image.out.
    
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
    
    '''
    
    remove_source = int(remove_source)
    print '** Ray-tracing images now...'
    if not inputfilename:
        inputfilename=os.path.join(os.path.expanduser('~'),'MCMax',path_mcmax,\
                                   'models','inputMCMax_%s.dat'%model_id)
    output_folder = os.path.join(os.path.expanduser("~"),'MCMax',path_mcmax,\
                                 'models',model_id)
    image_file = os.path.join(os.path.expanduser('~'),'MCMax','src',\
                              'Image.out')
    if remove_source:
        subprocess.call([' '.join(['MCMax ' + inputfilename,'0',\
                                   '-s tracestar=.false.','-o',output_folder,\
                                   image_file])],shell=True)
    else:
        subprocess.call([' '.join(['MCMax ' + inputfilename,'0','-o',\
                                   output_folder,image_file])],shell=True)
    print '** Your images can be found at:'
    print output_folder
                                


class MCMax(ModelingSession):
    
    """ 
    Class that includes all methods required for creating an MCMax model. 
    
    """
    
    def __init__(self,path_mcmax='runTest',replace_db_entry=0,
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode')):
        
        """ 
        Initializing an instance of ModelingSession.
        
        @keyword replace_db_entry: replace an entry in the MCMax database with 
                                   a newly calculated model with a new model id 
                                   (for instance if some general data not 
                                   included in the inputfiles is changed)
                                   
                                   (default: 0)
        @type replace_db_entry: bool
        @keyword path_mcmax: modeling folder in MCMax home
        
                             (default: 'runTest')
        @type path_mcmax: string
        @keyword path_combocode: CC home folder
        
                                 (default: /home/robinl/ComboCode/')
        @type path_combocode: string
        
        """
        
        super(MCMax, self).__init__(code='MCMax',path=path_mcmax,\
                                    path_combocode=path_combocode,\
                                    replace_db_entry=replace_db_entry)
        DataIO.testFolderExistence(os.path.join(os.path.expanduser("~"),\
                                   'MCMax',self.path,'data_for_gastronoom'))
        self.done_mcmax = False
        
        #- Read standard input file with all parameters that should be included
        inputfilename = os.path.join(os.path.expanduser("~"),'MCMax','src',\
                                     'inputMCMax.dat')
        self.standard_inputfile = [line
                                   for line in DataIO.readFile(inputfilename)]
        self.standard_inputfile = [line.split('=') 
                                   for line in self.standard_inputfile 
                                   if line[0] != '*']
        
        

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



    def doMCMax(self,star):
        
        """
        Running MCMax.
        
        @param star: The parameter set for this session
        @type star: Star()
        
        """

        print '***********************************'
        
        #Create command for running mcmax for every STAR instance
        print '** Making input file for MCMax'
        self.command_list = []
        for line in self.standard_inputfile:                                                      ## R_SOLAR = 0.004652 AU
            if line[0].upper()=='STARTYPE':
                if star['STARTYPE'] == 'BB':
                    self.command_list.append('%s=%s' %('Tstar',str(float(star['T_STAR']))))
                    self.command_list.append('%s=%s' %('Rstar',str(float(star['R_STAR']))))  
                elif star['STARTYPE'] == 'MARCS':
                    #self.command_list.append("%s='%s'" %(line[0],star['STARTYPE']))
                    #self.command_list.append("%s='%s'" %('starfile',os.path.join(os.path.expanduser('~'),self.path_combocode,\
                    #                                                                'StarFiles',star['STARFILE'])))
                    #self.command_list.append("%s=%s" %('Lstar',star['L_STAR']))
                    ##if self.prepMarcs(starParsList[i][0],marcstypestr,vlsr,MARCS_kernel,useMARCS):      #Returns a boolean,makes a file with data (see self.command_list for format)
                    #self.command_list = [string for string in self.command_list if string.find('Tstar') == -1 and string.find('Rstar') == 1]
                    #self.command_list.extend([" startype='FILE'",'_'.join("starfile='marcs",star['MARCS_TYPE'],star['MARCS_KERNEL'],star['T_STAR']) + ".dat'",'Lstar=' + str(star['L_STAR']) + 'd0'])
                    raise IOError('STARTYPE=MARCS not yet supported.')
                elif star['STARTYPE'] == 'FILE':
                    self.command_list.append("%s='%s'" %(line[0],star['STARTYPE']))
                    self.command_list.append("%s='%s'" %('starfile',os.path.join(os.path.expanduser('~'),self.path_combocode,\
                                                                                 'StarFiles',star['STARFILE'])))
                    self.command_list.append("%s=%s" %('Lstar',star['L_STAR']))
            elif line[0].upper()=='RIN':
                if int(star['TDESITER']): 
                    pass
                    #self.command_list.append('%s=%s' %(line[0],str(float(star['R_STAR'])*0.004652)))
                else:
                    self.command_list.append('%s=%s' %(line[0],str(float(star['R_INNER_DUST'])*float(star['R_STAR'])*0.004652)))
            elif line[0].upper()=='ROUT':
                self.command_list.append('%s=%s' %(line[0],str(float(star['R_OUTER_DUST'])*float(star['R_STAR'])*0.004652)))
            elif line[0].upper()=='MDOT':
                if star['DENSTYPE'] == "MASSLOSS":
                    self.command_list.append('%s=%s' %(line[0],str(float(star['MDOT_DUST'])*100)))                
            elif line[0].upper()=='VEXP':
                if star['DENSTYPE'] == "MASSLOSS":
                    self.command_list.append('%s=%s' %(line[0],str(float(star['V_EXP_DUST']))))                
            elif line[0].upper()=='TDESITER':
                self.command_list.append('%s=%s' %(line[0],int(star['TDESITER']) and '.true.' or '.false.'))     
            elif line[0].upper()=='FLD':
                self.command_list.append('%s=%s' %(line[0],int(star['FLD']) and '.true.' or '.false.'))     
            elif line[0].upper()=='TCONTACT':
                self.command_list.append('%s=%s' %(line[0],int(star['T_CONTACT']) and '.true.' or '.false.'))        
            elif line[0].upper()=='TMAX': 
                self.command_list.append('%s=%s' %(line[0],str(float(star['T_INNER_DUST'])))) 
            elif line[0].upper()=='DENSTYPE':
                self.command_list.append('%s=%s' %(line[0],"'%s'"%star['DENSTYPE']))
                if star['DENSTYPE'] == "SHELLFILE" or star['DENSTYPE'] == "FILE":
                    self.command_list.append('%s=%s' %('densfile',star['DENSFILE']))
                if star['DENSTYPE'] == "POW":
                    self.command_list.append('%s=%s' %('denspow',star['DENSPOW']))
                    self.command_list.append('%s=%s' %('mdust',star['M_DUST']))
            else:
                try:
                    self.command_list.append('%s=%s' %(line[0],star[line[0].upper()]))
                except KeyError:
                    try:
                        self.command_list.append('%s=%s' %(line[0],star[line[0].upper().replace(line[0].upper()[0],line[0].upper()[0] + '_',1)]))
                    except KeyError:
                        try:
                            self.command_list.append('%s=%s' %(line[0],star[line[0].upper().replace(line[0].upper()[0],line[0].upper()[0] + '_',1) + '_DUST']))
                        except KeyError:
                            try:
                                self.command_list.append('%s=%s' %(line[0],star[line[0].upper()+ '_DUST']))
                            except KeyError:
                                if line[0][:4] != 'abun' and line[0][:4] != 'part':
                                    self.command_list.append('%s=%s' %(line[0],line[1]))
                                    #print line[0].upper() + ' not found in inputComboCode.dat. Taking standard value from inputMCMax.dat: ' + line[1] + '.'
        
        #Following five keywords are only added if the corresponding abundance is present, ie not zero
        dust_list = DataIO.getInputData(path=os.path.join(self.path_combocode,'Data'),keyword='SPECIES_SHORT',filename='Dust.dat')
        kappas = DataIO.getInputData(path=os.path.join(self.path_combocode,'Data'),keyword='PART_FILE',filename='Dust.dat')
        
        self.command_list.extend(['abun%.2i=%s'%(index+1,star['A_' + species]) 
                                  for index,species in enumerate(star['DUST_LIST'])])
        self.command_list.extend(['%s%.2i=%s' %(kappas[dust_list.index(species)].find('.opacity') != -1 and 'opac' or 'part',\
                                                index+1,"'" + os.path.join(os.path.expanduser('~'),'MCMax','src',\
                                                kappas[dust_list.index(species)]) + "'")
                                  for index,species in enumerate(star['DUST_LIST'])])
        #This is always present in the inputfile 
        self.command_list.extend(['TdesA%.2i=%s' %(index+1,str(float(star['T_DESA_' + species])))
                                  for index,species in enumerate(star['DUST_LIST'])])
        #This is always present in the inputfile 
        self.command_list.extend(['TdesB%.2i=%s' %(index+1,str(float(star['T_DESB_' + species])))
                                  for index,species in enumerate(star['DUST_LIST'])])
        #R_MIN can only be added if it is present in inputComboCode.dat, hence only if star has the key, it is NEVER made
        self.command_list.extend(['minrad%.2i=%s' %(index+1,str(float(star['R_MIN_' + species])*float(star['R_STAR'])*0.004652))
                                  for index,species in enumerate(star['DUST_LIST'])
                                  if star.has_key('R_MIN_' + species)])
        #R_MAX can be derived from T_MIN, hence it is always created and if T_MIN is not present, an empty string is added
        self.command_list.extend(['maxrad%.2i=%s' %(index+1,str(float(star['R_MAX_' + species])*float(star['R_STAR'])*0.004652))
                                  for index,species in enumerate(star['DUST_LIST'])
                                  if star['R_MAX_' + species]])
        self.command_list[0:0] = [str(star['PHOTON_COUNT'])]
        print '** DONE!'
        print '***********************************'
        #command line doesnt include output directory, such that all the rest can be used to compare with database entries
        #model id is also included with command list in db, an id is based on time of calculation in UTC in seconds siunce 1970
        self.model_id = ''
        try:
            mcmax_db = open(os.path.join(os.path.expanduser('~'),'MCMax',self.path,'MCMax_models.db'),'r')
            try:
                while True:
                    model = cPickle.load(mcmax_db)
                    if model[0] == self.command_list:
                        if self.replace_db_entry and 0:
                            mcmax_db.close()
                            Database.deleteModel(model_id=str(model[1]),db_path=os.path.join(os.path.expanduser('~'),'MCMax',self.path,'MCMax_models.db'))
                            raise EOFError
                        else:
                            print '** MCMax model has been calculated before with ID ' + str(model[1]) + '.'
                            self.model_id = model[1]
                        break
            except EOFError:
                print '** No match found in MCMax database. Calculating new model.'
            finally:
                mcmax_db.close()
        except IOError:
            print '** No database present at ' + os.path.join(os.path.expanduser('~'),'MCMax',self.path,'MCMax_models.db') + '. Creating new one.'
        #if no match is found between command lines, calculate new model with new model id 
        #if the calculation did not fail, append an entry to the database of the new model, with its id     
        if not self.model_id:
            self.model_id = 'model_%.4i-%.2i-%.2ih%.2i-%.2i-%.2i' %(gmtime()[0],gmtime()[1],gmtime()[2],gmtime()[3],gmtime()[4],gmtime()[5])
            input_filename = os.path.join(os.path.expanduser("~"),'MCMax',self.path,'models','inputMCMax_%s.dat'%self.model_id)
            DataIO.writeFile(filename=input_filename,input_lines=self.command_list[1:])
            subprocess.call([' '.join(['MCMax',input_filename,self.command_list[0],'-o',\
                            os.path.join(os.path.expanduser("~"),'MCMax',self.path,'models',self.model_id)])],shell=True)
            mcmax_db = open(os.path.join(os.path.expanduser('~'),'MCMax',self.path,'MCMax_models.db'),'a')
            try:
                new_model = open(os.path.join(os.path.expanduser("~"),'MCMax',self.path,'models',self.model_id,'MCScattered.dat'), 'r')
                new_model2 = open(os.path.join(os.path.expanduser("~"),'MCMax',self.path,'models',self.model_id,'denstemp.dat'), 'r')
                new_model.close()
                new_model2.close()
                cPickle.dump([self.command_list,self.model_id],mcmax_db)
                self.done_mcmax = True
            except IOError:
                print '** Model calculation failed. No entry is added to the database and LAST_MCMAX_MODEL in STAR dictionary is not updated.'
                self.model_id = ''
            finally:
                mcmax_db.close()
        #clean up: add/change 'LAST_MCMAX_MODEL' entry, check T_MIN, delete some variables that should be changed
        if self.model_id:
            #add/change MUTABLE input keys, which will not be overwritten by the input file inputComboCode.dat
            star['LAST_MCMAX_MODEL'] = self.model_id
            #add/change MUTABLE input keys, which MAY be overwritten by the input file inputComboCode.dat
            #Done in ModelingManager now, after the star instance has been copied to a backup in the case of iterative              
        print '***********************************'
        
