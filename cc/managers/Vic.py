# -*- coding: utf-8 -*-

"""
Interface for communicating with and managing modeling through VIC3.

Author: R. Lombaert

"""

import os
from scipy import log10
import subprocess
from time import gmtime, sleep

from cc.tools.io import DataIO
from cc.modeling.codes import Gastronoom



class Vic():
    
    """ 
    Creating a vic manager which communicates with the Vic3 supercomputer and
    updates modeling results on the home disk.
    
    """
    
    def __init__(self,account,path='code23-01-2010',credits_acc=None,\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode'),\
                 time_per_sphinx=30,recover_sphinxfiles=0):
        
        """ 
        Initializing a Vic instance.
        
        @param account: the name of the user account at VIC
                          
                          (default: 'vsc30226')
        @type account: string
        
        @keyword path: The output folder in the GASTRoNOoM home folder
        
                            (default: 'codeSep2010')
        @type path: string
        @keyword credits_acc: the name of the credits account to be charged. If
                              one, the standard user account is charged.
                              
                              (default: None)
        @type credits_acc: string
        @keyword time_per_sphinx: the expected calculation time for one sphinx 
                                  model in minutes
                                  
                                  (default: 30)
        @type time_per_sphinx: int
        @keyword recover_sphinxfiles: Try to recover sphinx files from the Vic
                                      disk in case they were correctly 
                                      calculated, but not saved to the database
                                      for one reason or another. 
                                      
                                      (default: 0) 
        @type recover_sphinxfiles: bool
        @keyword path_combocode: CC home folder
        
                                 (default: '~/ComboCode/')
        @type path_combocode: string 

        """

        self.finished = dict()
        self.failed = dict()
        self.models = dict()
        self.command_lists = dict()
        self.transitions = dict()
        self.sphinx_model_ids = dict()
        self.inputfiles = dict()
        self.trans_in_progress = [] 
        self.path = path
        self.path_combocode = path_combocode
        self.account = account
        self.disk = account[3:6]
        self.uname = os.path.split(os.path.expanduser('~')+'/'.rstrip('/'))[-1]
        self.time_per_sphinx = float(time_per_sphinx)
        self.recover_sphinxfiles = recover_sphinxfiles
        self.current_model = 0
        if not credits_acc:
            self.credits_acc = None
        else:
            self.credits_acc = credits_acc
        main_vic_folder = os.path.join('/data','leuven',self.disk,\
                                       self.account,'COCode')
        subprocess.call('ssh %s@login.vic3.cc.kuleuven.be mkdir %s/output/'\
                        %(self.account,main_vic_folder),shell=True)
        subprocess.call('ssh %s@login.vic3.cc.kuleuven.be mkdir %s'\
                        %(self.account,main_vic_folder) + '/dust_files/',\
                        shell=True)
        subprocess.call('ssh %s@login.vic3.cc.kuleuven.be mkdir %s'\
                        %(self.account,main_vic_folder) +'/CustomAbundances/',\
                        shell=True)
        subprocess.call('ssh %s@login.vic3.cc.kuleuven.be mkdir %s'\
                        %(self.account,main_vic_folder) +'/StarFiles/',\
                        shell=True)


    def setSphinxDb(self,sph_db):
        
        '''
        Set the Sphinx db for this Vic instance.
        
        @param sph_db: The sphinx database
        @type sph_db: Database()
        
        '''
        
        self.sph_db = sph_db    

   
    
    def updateLineSpec(self):
        
        '''
        Update telescope.spec files on VIC3.
        
        '''
        
        homespec = os.path.join(os.path.expanduser('~'),'GASTRoNOoM','src',\
                                'data','*spec')
        vicspec = os.path.join('/user','leuven',self.disk,self.account,\
                               'COCode','data','.')
        subprocess.call(['scp %s %s@login.vic3.cc.kuleuven.be:%s'\
                         %(homespec,self.account,vicspec)],shell=True)
        


    def addModel(self,model_id,command_list):
        
        '''
        Add model to the list of to be processed models.
        
        Here, dictionaries are initiated to contain information about the
        modeling session. 
        
        Every entry in the dictionaries have an index number associated with 
        them to uniquely identify a modeling session across the vic manager.
        
        Dictionaries are kept for model ids, parameter sets, transitions, 
        failed calculations and finished calculations, and vic inputfiles for 
        every transition.
        
        The current_model index is the same between a call to the addModel 
        method and the queueModel() OR the reset methods. Anything between uses
        the same index. queueModel() will move to the next index value, while
        reset() will reset the current index number.
        
        @param model_id: The cooling model_id
        @type model_id: string
        @param command_list: The parameters for this GASTRoNOoM model
        @type command_list: dict()
         
        '''
        
        if self.models.has_key(self.current_model):
            raise Error('Vic().addModel() is trying to add a model_id ' + \
                        'to a session that was already assigned an id. ' + \
                        'Reset or queue the previous model first.')
        self.models[self.current_model] = model_id
        self.command_lists[self.current_model] = command_list
        self.transitions[self.current_model] = []
        self.failed[self.current_model] = []
        self.finished[self.current_model] = []
        self.inputfiles[self.current_model] = []



    def addTransInProgress(self,trans):
        
        '''
        Add a transition to the list of transitions in progress. They will be
        checked at the end of the VIC run to see if they have been correctly 
        calculated.
        
        This concerns transitions that are requested, but are already present 
        in the sphinx database with an "IN_PROGRESS" keyword included in the 
        transition dictionary. These will not be calculated again, instead they 
        are remembered and checked at the end of the full modeling run.
        
        @param trans: The transition calculating on Vic3
        @type trans: Transition()
        
        '''
        
        self.trans_in_progress.append(trans)
        


    def addTrans(self,trans):
        
        '''
        Add a transition to the vic transition list.
        
        These will be calculated on Vic3 and inputfiles will be prepared for 
        them.
        
        The entries in the transitions dictionary are deleted when the entry 
        has been completely finished. 
        
        @param trans: The transition to be calculated on Vic3
        @type trans: Transition()
        
        '''
        
        self.transitions[self.current_model].append(trans)
                


    def queueModel(self):
        
        '''
        Queue the current model on VIC3, which will cp mline and cooling output 
        file to VIC, run the shell script there, and create the necessary input 
        files.
        
        Once everything has been started up, the current model index number is 
        increased by one to allow for a new model to be added. You cannot add a 
        model and then add another model unless you reset or queue the previous
        model.
        
        '''
        
        self.sphinx_model_ids[self.current_model] \
            = list(set([trans.getModelId() 
                        for trans in self.transitions[self.current_model]]))
        if not self.recover_sphinxfiles:
            printing = self.makeJobFile()
            self.makeInputFiles()
            model_id = self.models[self.current_model]
            jobfile = os.path.join('/user','leuven',self.disk,self.account,\
                                   'COCode','vic_run_jobs_%s_%i.sh'\
                                   %(model_id,self.current_model))
            subprocess.Popen('ssh %s@login.vic3.cc.kuleuven.be %s'\
                             %(self.account,jobfile),shell=True,\
                             stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            if printing: print '\n'.join(printing)
        self.current_model += 1


  
    def reset(self):
        
        '''
        If a model has been added, and no transitions were required to be 
        calculated, remove that model entry here.
        
        This only removes the self.models and self.transitions entries. The 
        other dictionaries are re-initiated anyway when adding a model. 
       
        '''
        
        del self.models[self.current_model]
        del self.transitions[self.current_model]
        


    def makeJobFile(self):
        
        '''
        Make the job file that will run the loop on VIC3 and copy the cooling
        and mline output to VIC3.
        
        @return: to be printed strings once all the copying is done, which 
                 shows how many transitions are being calculated for which 
                 sphinx model id
        @rtype: list[string]
        
        '''
        
        model_id = self.models[self.current_model]
        vic_server = '%s@login.vic3.cc.kuleuven.be'%self.account
        jobfiles = []
        printing = []
        for model_id_sphinx in self.sphinx_model_ids[self.current_model]:
            these_trans = [trans 
                           for trans in self.transitions[self.current_model] 
                           if trans.getModelId() == model_id_sphinx]
            models_in_job = (int(log10(len(these_trans))))**2 
            if not models_in_job: models_in_job = 1
            job_number = len(these_trans)%models_in_job==0.0 \
                            and len(these_trans)/models_in_job \
                            or int(len(these_trans)/models_in_job)+1
            #- job_number: this is the number of jobs you want to queue
            job_number = job_number/8+1  
            time_per_job = self.time_per_sphinx*models_in_job
            walltimestring = '%.2i:00:00'%(int(time_per_job/60)+1)
            
            #- Create job file
            jobfile = DataIO.readFile(os.path.join(os.path.expanduser('~'),\
                                      'GASTRoNOoM','vic_job_example.sh'))
            new_jobfile = []
            for line in jobfile:
                #if line.find('#PBS -l nodes=1:ppn=8') != -1:
                #    new_line = line.replace('ppn=8','ppn=%i'%)
                if line.split('=')[0].find('#PBS -l walltime') != -1:
                    new_line = '='.join([line.split('=')[0],walltimestring])
                elif line.split('=')[0].find('export COCODEHOME') != -1:
                    new_line = line.replace('vsc30226',self.account)\
                                        .replace('/302/','/%s/'%self.disk)
                elif line.split('=')[0].find('export COCODEDATA') != -1:
                    new_line = '='.join([line.split('=')[0],\
                            os.path.join(line.split('=')[1],model_id+'_'+\
                            str(self.current_model)+'/')])\
                            .replace('vsc30226',self.account)\
                            .replace('/302/','/%s/'%self.disk)
                elif line.find('for i in $(seq 1 1)') != -1:
                    new_line = 'for i in $(seq 1 %i)'%models_in_job
                elif line.split('=')[0].find('export MODELNUMBER') != -1:
                    new_line = '='.join([line.split('=')[0],line.split('=')[1]\
                                    .replace('model_id',model_id_sphinx)])
                else:
                    new_line = line
                new_jobfile.append(new_line)
                
            #- Save job file, change permission and copy to VIC
            local_folder = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                        self.path,'models',model_id_sphinx)
            jobfilename_vic = '/user/leuven/'+self.disk+'/' + self.account + \
                              '/COCode/vic_job_' + model_id_sphinx + '.sh'
            jobfilename_local = os.path.join(local_folder,'vic_input',\
                                             'vic_job_%s.sh'%model_id_sphinx)
            DataIO.writeFile(jobfilename_local,new_jobfile)
            subprocess.call(['chmod +x %s'%jobfilename_local],shell=True)
            subprocess.call(['scp %s %s:%s'%(jobfilename_local,vic_server,\
                                             jobfilename_vic)],shell=True)
            jobfiles.append((jobfilename_vic,job_number))  
            
            #- Make the output folder on VIC
            vic_folder = '/data/leuven/%s/%s/COCode/output/%s/'\
                         %(self.disk,self.account,model_id_sphinx)
            subprocess.call('ssh %s mkdir %s'%(vic_server,vic_folder),\
                            shell=True)
            
            #-copy required GASTRoNOoM files, molecule specific.
            these_molecules = set(['sampling'] + \
                                  [trans.molecule.molecule 
                                   for trans in these_trans])
            to_be_copied = ['coolfgr*','input%s.dat'%model_id_sphinx]
            to_be_copied.extend(['cool*_%s.dat'%molec 
                                 for molec in these_molecules])
            to_be_copied.extend(['ml*_%s.dat'%molec 
                                 for molec in these_molecules 
                                 if molec != 'sampling'])
            for filecopy in to_be_copied:
                subprocess.call(['scp %s %s:%s.'\
                                 %(os.path.join(local_folder,filecopy),\
                                   vic_server,vic_folder)], shell=True)
            
            #- number of nodes*number of cpus=amount of times to queue it
            printing.append('Running %i jobs with %i models each for ID %s.' \
                            %(job_number*7, models_in_job,model_id_sphinx))
        
        #- Create the run-jobs file.
        runjobsfile = DataIO.readFile(os.path.join(os.path.expanduser('~'),\
                                      'GASTRoNOoM','vic_run_jobs_example.sh'))
        new_runjobsfile = []
        for i,job in enumerate(jobfiles):
            for line in runjobsfile:
                if line.find('#!/bin/bash -l') != -1 and i == 0:
                    new_line = line
                elif line.find('cd COCode') != -1 and i == 0:
                    new_line = 'cd /user/leuven/'+self.disk+'/'+ self.account+\
                               '/COCode/'
                elif line.find('module load worker') != -1 and i == 0:
                    new_line = line
                elif line.find('for i in $(seq 1 1) ;') != -1:
                    new_line = 'for i in $(seq 1 %i) ;' %(job[1])    
                elif line.find('do wsub -A dummy -t 1-8 -batch') != -1:
                    new_line = line.replace('vic_job_example.sh',\
                                            os.path.split(job[0])[1])   
                    if self.credits_acc <> None:
                        new_line = new_line.replace('-A dummy',\
                                                    '-A %s'%self.credits_acc)
                    else:
                        new_line = new_line.replace('-A dummy ','')
                elif line.find('done') != -1:
                    new_line = line
                else:
                    new_line = ''
                new_runjobsfile.append(new_line)
        
        #- Run-jobs file: Write, change permission and copy to VIC
        runjobsfilename_local = os.path.join(os.path.expanduser('~'),\
                                             'GASTRoNOoM',self.path,'models',\
                                             model_id,'vic_input',\
                                             'vic_run_jobs_%s_%s.sh'\
                                             %(model_id,\
                                               str(self.current_model)))
        runjobsfilename_vic = '/user/leuven/%s/%s/COCode/vic_run_jobs_%s_%s.sh'\
                              %(self.disk,self.account,model_id,\
                                str(self.current_model))
        DataIO.writeFile(runjobsfilename_local,new_runjobsfile)
        subprocess.call(['chmod +x %s'%runjobsfilename_local],shell=True)
        subprocess.call(['scp %s %s:%s'\
                         %(runjobsfilename_local,vic_server,\
                           runjobsfilename_vic)],shell=True)
        return printing
                
  

    def makeInputFiles(self):
        
        '''
        Make the input files with just one line request in each.
        
        These inputfiles are then converted to format appropriate for Vic3 and
        subsequently copied to Vic3.
        
        '''
        
        model_id = self.models[self.current_model] 
        vic_model_folder = os.path.join('/data','leuven',self.disk,\
                                        self.account,'COCode','%s_%i/'\
                                        %(model_id,self.current_model))
        subprocess.call('ssh %s@login.vic3.cc.kuleuven.be mkdir %s'\
                        %(self.account,vic_model_folder),shell=True)
        will_calculate_stuff = 0
        custom_files = []
        opacity_files = []
        starfiles = []
        for model_id_sphinx in self.sphinx_model_ids[self.current_model]:
            these_trans = [trans 
                           for trans in self.transitions[self.current_model]
                           if trans.getModelId() == model_id_sphinx]
            for i,trans in enumerate(these_trans):
                will_calculate_stuff = 1
                actual_command_list \
                    = self.command_lists[self.current_model].copy()
                actual_command_list['DATA_DIRECTORY'] \
                    = '"/user/leuven/%s/%s/COCode/data/"'\
                      %(self.disk,self.account)
                actual_command_list['OUTPUT_DIRECTORY'] \
                    = '"/data/leuven/%s/%s/COCode/output/%s/"'\
                      %(self.disk,self.account,trans.getModelId())
                actual_command_list['PARAMETER_FILE'] \
                    = '"/data/leuven/%s/%s/COCode/output/%s/%s"'\
                      %(self.disk,self.account,trans.getModelId(),\
                        'parameter_file_%s.dat'%trans.getModelId())
                actual_command_list['OUTPUT_SUFFIX'] = trans.getModelId()
                opacity_files.append(actual_command_list['TEMDUST_FILENAME'])
                if int(actual_command_list['KEYWORD_DUST_TEMPERATURE_TABLE']):
                    actual_command_list['DUST_TEMPERATURE_FILENAME'] \
                        = '"/data/leuven/%s/%s/COCode/dust_files/%s"'\
                          %(self.disk,self.account,\
                            os.path.split(actual_command_list\
                                          ['DUST_TEMPERATURE_FILENAME']\
                                          .strip('"'))[1])
                path = '/data/leuven/%s/%s/COCode/CustomFiles/'\
                       %(self.disk,self.account)
                molec_dict = trans.molecule.makeDict(path)
                starfiles.append(molec_dict.pop('STARFILE',''))
                commandfile = \
                     ['%s=%s'%(k,v) 
                      for k,v in sorted(actual_command_list.items()) 
                      if k != 'R_POINTS_MASS_LOSS'] + ['####'] + \
                     ['%s=%s'%(k,v) 
                      for k,v in sorted(molec_dict.items())] + ['####'] + \
                     ['%s=%s'%(k,v) 
                      for k,v in sorted(trans.makeDict().items())] + \
                     ['######']
                for key,fkey in zip(['ENHANCE_ABUNDANCE_FACTOR',     
                                     'SET_KEYWORD_CHANGE_ABUNDANCE',
                                     'SET_KEYWORD_CHANGE_TEMPERATURE'],\
                                    ['ABUNDANCE_FILENAME',\
                                     'CHANGE_FRACTION_FILENAME',\
                                     'NEW_TEMPERATURE_FILENAME']):
                     if getattr(trans.molecule,key.lower()):
                          custom_files.append((getattr(trans.molecule,\
                                                       fkey.lower()),\
                                               molec_dict[fkey].strip('"')))
                if actual_command_list.has_key('R_POINTS_MASS_LOSS'):
                    commandfile.extend(['%s=%s'%('R_POINTS_MASS_LOSS',v) 
                                        for v in actual_command_list\
                                                 ['R_POINTS_MASS_LOSS']] + \
                                       ['####'])
                infile = '_'.join(['gastronoom',trans.getModelId(),\
                                   '%i.inp'%(i+1)])
                DataIO.writeFile(os.path.join(os.path.expanduser('~'),\
                                              'GASTRoNOoM',self.path,'models',\
                                              model_id,'vic_input',\
                                              infile),commandfile)
                self.inputfiles[self.current_model]\
                    .append(os.path.join('/data','leuven',self.disk,\
                                         self.account,'COCode','%s_%i'\
                                         %(model_id,self.current_model),\
                                         infile.replace('.inp','.*')))
            #- There is no overlap between filenames: All filenames with the 
            #- same trans model id get an increasing number i
            #- Then they are copied to VIC3 and removed afterwards, so there is 
            #- never a problem to store them in the cooling model_id folder
        if not will_calculate_stuff:
            return
        else:
            starfiles = list(set([f for f in starfiles if f]))
            if len(starfiles) > 1: 
                print('WARNING! Multiple starfiles detected in grid in Vic.py!')
            if starfiles:
                path = os.path.join('/data','leuven',self.disk,self.account,\
                                 'COCode','StarFiles','starfile_tablestar.dat')
                subprocess.call(['scp ' + starfiles[0] + ' ' + self.account + \
                                 '@login.vic3.cc.kuleuven.be:' + \
                                 path],shell=True)
            opacity_files = [f 
                             for f in set(opacity_files)
                             if f != 'temdust.kappa']
            for filename in opacity_files:
                full_path = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                         'src','data',filename)
                full_path_vic = '/user/leuven/'+self.disk+'/'+self.account + \
                                '/COCode/data/.'
                subprocess.call(['scp ' + full_path + ' ' + self.account + \
                                 '@login.vic3.cc.kuleuven.be:'+full_path_vic],\
                                shell=True)
            for filename,vicfile in set(custom_files):
                subprocess.call(['scp ' + filename + ' ' + self.account + \
                                 '@login.vic3.cc.kuleuven.be:' + vicfile],\
                                shell=True)
            if int(self.command_lists[self.current_model]\
                                     ['KEYWORD_DUST_TEMPERATURE_TABLE']):
                homefile = self.command_lists[self.current_model]\
                                             ['DUST_TEMPERATURE_FILENAME']\
                                             .strip('"')
                vicfile = actual_command_list['DUST_TEMPERATURE_FILENAME']\
                                             .strip('"')
                subprocess.call(['scp %s %s@login.vic3.cc.kuleuven.be:%s'\
                                 %(homefile,self.account,vicfile)],\
                                shell=True)
            homeinput = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                     self.path,'models',model_id,'vic_input',\
                                     'gastronoom_*.inp')
            vicinput = os.path.join('/data','leuven',self.disk,self.account,\
                                    'COCode',\
                                    '%s_%i/.'%(model_id,self.current_model))
            subprocess.call(['scp %s %s@login.vic3.cc.kuleuven.be:%s'\
                             %(homeinput,self.account,vicinput)],shell=True)
            subprocess.call(['rm %s'%homeinput],shell=True)
        
   
     
    def finalizeVic(self):
        
        '''
        Finalize a modeling procedure on VIC: successful and failed results 
        are printed to a file, including the transitions.
        
        This log file can be used as input for ComboCode again by putting 
        LINE_LISTS=2.
        
        '''
        
        for trans in self.trans_in_progress:
            filename = os.path.join(os.path.expanduser('~'),'GASTRoNOoM',\
                                    self.path,'models',trans.getModelId(),\
                                    trans.makeSphinxFilename(2))
            if not os.path.isfile(filename):              
                trans.setModelId('') 
        if self.models.keys():
            time_stamp = '%.4i-%.2i-%.2ih%.2i:%.2i:%.2i' \
                         %(gmtime()[0],gmtime()[1],gmtime()[2],\
                           gmtime()[3],gmtime()[4],gmtime()[5])
            results = ['# Successfully calculated models:'] \
                    + [self.models[current_model] 
                       for current_model in self.models.keys() 
                       if current_model not in self.failed.keys()] \
                    + ['# Unsuccessfully calculated models (see 3 logfiles '+ \
                       'for these models):'] \
                    + [self.models[current_model] 
                       for current_model in self.models.keys() 
                       if current_model in self.failed.keys()]
            DataIO.writeFile(os.path.join(os.path.expanduser('~'),\
                                          'GASTRoNOoM',self.path,\
                                          'vic_results','log_' + time_stamp),\
                             results)
            for current_model,model_id in self.models.items():
                model_results = ['# Successfully calculated transitions:'] + \
                    ['Sphinx %s: %s' %(trans.getModelId(),str(trans)) 
                     for trans in self.finished[current_model]] + \
                    ['# Unsuccessfully calculated transitions (see 2 other ' + \
                     'logfiles for these transitions):'] + \
                    ['Sphinx %s: %s' %(trans.getModelId(),str(trans)) 
                     for trans in self.failed[current_model]]
                DataIO.writeFile(os.path.join(os.path.expanduser('~'), \
                                              'GASTRoNOoM',self.path,\
                                              'vic_results','log_results%s_%i'\
                                              %(time_stamp,current_model)),\
                                 model_results)
                for this_id in self.sphinx_model_ids[current_model]:
                    sphinx_files = os.path.join(os.path.expanduser('~'),\
                                                'GASTRoNOoM',self.path,\
                                                'models',this_id,'sph*')
                    subprocess.call(['chmod a+r %s'%sphinx_files],shell=True)
        


    def checkProgress(self,wait_qstat=0):
        
        '''
        Checks progress on all queued model_ids by checking if the output has 
        been copied to the local disk from VIC3.
        
        In this method, the self.failed and self.finished keywords are updated,
        and the self.inputfiles and self.transitions is cleaned up on the go.
        
        @keyword wait_qstat: wait 10 seconds before checking the qstat query on
                             vic, in order to make sure that the qeueu command 
                             on vic is finished queueing up the new models, as 
                             this may be slower than the python script. Use 
                             this for progress check after just queueing a new 
                             model
                             
                             (default: 0)
        @type wait_qstat: bool
        
        @return: Have all models been finished on Vic3? 
        @rtype: bool
        
        '''
        
        #- iteration op self.transitions, not self.models since the latter is 
        #- not updated as progress is checked...
        for current_model in self.transitions.keys():     
            model_id = self.models[current_model]
            print 'Currently checking %s...' %model_id
            lsprocess = subprocess.Popen('ssh %s@login.vic3.cc.kuleuven.be '\
                                         %self.account + \
                                         'ls /data/leuven/%s/%s/COCode/%s_%i/'\
                                         %(self.disk,self.account,model_id,\
                                           current_model),\
                                         shell=True,stdout=subprocess.PIPE)
            lsfile = lsprocess.communicate()[0].split('\n')
            lsfile = [f 
                      for f in lsfile 
                      if f[-9:] != '.inp.done' and f != 'jobs.log' and f]
            if wait_qstat:
                sleep(10)    
            qstatprocess = subprocess.Popen(\
                    'ssh %s@login.vic3.cc.kuleuven.be qstat | grep %s'\
                    %(self.account,self.account),shell=True,\
                    stdout=subprocess.PIPE)
            qstatfile = qstatprocess.communicate()[0].split('\n')
            if not lsfile or qstatfile == ['']:
                for model_id_sphinx in self.sphinx_model_ids[current_model]:
                    trans_strings = [trans.makeSphinxFilename() 
                                     for trans in self.transitions\
                                                        [current_model]
                                     if trans.getModelId() == model_id_sphinx]
                    heresph = os.path.join(os.path.expanduser('~'),\
                                           'GASTRoNOoM',self.path,\
                                           'models',model_id_sphinx)
                    for trans_string in trans_strings:
                        vicsph = '/data/leuven/%s/%s/COCode/output/%s/%s'\
                                 %(self.disk,self.account,model_id_sphinx,\
                                   trans_string)
                        subprocess.call('scp %s@login.vic3.cc.kuleuven.be:%s '\
                                        %(self.account,vicsph)+'%s/.'%heresph,\
                                        shell=True,stdout=subprocess.PIPE,\
                                        stderr=subprocess.PIPE)
                    vicjob = '/data/leuven/%s/%s/COCode/%s_%i/jobs.log'\
                             %(self.disk,self.account,model_id,current_model)
                    subprocess.call('scp %s@login.vic3.cc.kuleuven.be:%s '\
                                    %(self.account,vicjob) + \
                                    '%s/log_vic_jobs'%(heresph),\
                                    shell=True,stdout=subprocess.PIPE)
                    gas_session = Gastronoom.Gastronoom(\
                                        path_combocode=self.path_combocode,\
                                        path_gastronoom=self.path,\
                                        sph_db=self.sph_db)
                    gas_session.model_id = model_id
                    gas_session.trans_list \
                        = [trans 
                           for trans in self.transitions[current_model] 
                           if trans.getModelId() == model_id_sphinx]
                    for trans in gas_session.trans_list:
                        gas_session.checkSphinxOutput(trans)
                    gas_session.finalizeSphinx()
                    self.finished[current_model] \
                        = [trans 
                           for trans in gas_session.trans_list 
                           if trans.getModelId()]
                    self.failed[current_model] \
                        = [trans 
                           for trans in gas_session.trans_list 
                           if not trans.getModelId()]
                    vicjob = '/user/leuven/%s/%s/COCode/vic_job_%s.sh*'\
                             %(self.disk,self.account,model_id_sphinx)
                    subprocess.call(['ssh %s@login.vic3.cc.kuleuven.be rm %s'\
                                     %(self.account,vicjob)],shell=True)
                del self.transitions[current_model]
                while self.inputfiles[current_model]:
                    subprocess.call([' '.join(\
                        ['ssh','%s@login.vic3.cc.kuleuven.be'%self.account,\
                         'rm'] + \
                        [self.inputfiles[current_model].pop() 
                         for i in range(len(self.inputfiles[current_model]) \
                                                < 100\
                                            and len(self.inputfiles\
                                                        [current_model])\
                                            or 100)])],\
                                    shell=True,stdout=subprocess.PIPE)
                vicjobfile = '/user/leuven/%s/%s/COCode/vic_run_jobs_%s_%i.sh'\
                             %(self.disk,self.account,model_id,current_model)
                subprocess.call(['ssh %s@login.vic3.cc.kuleuven.be rm %s'\
                                 %(self.account,vicjobfile)],shell=True)
        self.sph_db.sync()
        if self.transitions:
            return True
        else:
            return False

    
      
    def getQueue(self):
        
        '''
        Get a list of unique queue number + cooling model_id for those models
        that are still in progress on VIC.
        
        @return: The queue numbers and model_ids still in progress.
        @rtype: list[(int,string)]
        
        '''
        
        return [(k,v) 
                for k,v in self.models.items() 
                if self.transitions.has_key(k)]
                