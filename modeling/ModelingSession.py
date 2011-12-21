# -*- coding: utf-8 -*-

"""
Interface for creating modeling environments.

Author: R. Lombaert

"""

import os
from time import gmtime

from cc.tools.io import DataIO



class ModelingSession(object):
    
    """
    The basic modeling environment. Inherited by MCMax() and Gastronoom().
    
    """
      
    def __init__(self,code,path,replace_db_entry=0,\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode')):
        
        """ 
        Initializing an instance of ModelingSession.
        
        @param code: code for which the modelingsession is created
        @type code: string
        @param path: modeling output folder in the code's home folder
        @type path: string
        
        @keyword replace_db_entry: 
        @keyword path_combocode: CC home folder
          
                                 (default: '/home/robinl/ComboCode')
        @type path_combocode: string
        @keyword replace_db_entry: replace an entry in the database with a 
                                   newly calculated model with a new model id 
                                   (eg if some general data not included in 
                                   the inputfiles is changed)
                                              
                                   (default: 0)
        @type replace_db_entry: bool
          
        """
        
        self.path_combocode = path_combocode
        self.path = path
        self.code = code
        self.replace_db_entry = replace_db_entry
        mutablefile = os.path.join(self.path_combocode,'CC',\
                                   'Mutable_Parameters_' + code + '.dat')
        self.mutable = [line[0] 
                        for line in DataIO.readFile(mutablefile,delimiter=' ')
                        if ' '.join(line)]
        self.mutable = [line for line in self.mutable if line[0] != '#']
        DataIO.testFolderExistence(os.path.join(os.path.expanduser('~'),\
                                                self.code,self.path,'models'))
        


    def makeNewId(self):
        
        '''
        Make a new model_id based on the current UTC in seconds since 1970.
        
        '''
        
        return 'model_%.4i-%.2i-%.2ih%.2i-%.2i-%.2i' \
                %(gmtime()[0],gmtime()[1],gmtime()[2],\
                  gmtime()[3],gmtime()[4],gmtime()[5])
                  
                  
                  
    def setCommandKey(self,comm_key,star,key_type,star_key=None,\
                      alternative=None,make_int=0):
        
        '''
        Try setting a key in the command_list from a star instance. 
        
        If the key is unknown, it is left open and will be filled in from the 
        standard gastronoom inputfile.
        
        @param comm_key: the name of the keyword in the command list
        @type comm_key: string
        @param star: The parameter set
        @type star: Star()
        @param key_type: the type of the keyword, either 'DUST' or 'GAS' 
        @type key_type: string
        
        @keyword star_key: the name of the keyword in the star instance 
                           (minus '_%s'%key_type, which is added as well in a 
                           second attempt if the first without the addition is 
                           not found), if None, it is equal to comm_key
                           
                           (default: None)
        @type star_key: string
        @keyword alternative:  a default value passed from the standard 
                               inputfile that is used if the keyword or the 
                               keyword + '_%s'%key_type is not found in Star()
                               
                               (default: None)
        @type alternative: string
        @keyword make_int: make an integer before converting to string for this
                           keyword.
                           
                           (default: 0)
        @type make_int: boolean
        
        @return: True if successful, otherwise False.
        @rtype: bool
        
        ''' 
        
        if star_key is None: star_key = comm_key
        try:
            self.command_list[comm_key] = \
                    make_int \
                        and str(int(float(star[star_key]))) \
                        or str(star[star_key])
            return True
        except KeyError: 
            try:
                self.command_list[comm_key] = \
                    make_int \
                        and str(int(float(star[star_key+ '_%s'%key_type]))) \
                        or str(star[star_key+ '_%s'%key_type])
                return True
            except KeyError:
                if alternative <> None: 
                    self.command_list[comm_key] = \
                        make_int \
                            and str(int(float(alternative))) \
                            or str(alternative)
                    return True
                else:
                    return False



    def compareCommandLists(self,this_list,modellist,code,ignoreAbun=0,\
                            extra_dict=None):
        
        """
        Comparing a command_list with a database entry.
        
        @param this_list: parameters in this modeling session
        @type this_list: dict
        @param modellist: parameters from database model
        @type modellist: dict
        @param code: The GASTRoNOoM subcode
        @type code: string
        
        @keyword ignoreAbun: only relevant for mline: ignore the 4 abundance 
                             parameters (such as for co)
                             
                             (default: 0)
        @type ignoreAbun: bool
        @keyword extra_dict: if not None this gives extra dictionary entries 
                             to be used in the comparison on top of this_list.
                             The extra entries are assumed present in modellist
                             otherwise the comparison will return False.
                             
                             (default: None)
        @type extra_dict: dict
        
        @return: Comparison between the two parameter sets
        @rtype: bool
        
        """
        
        model_bool_list = []
        if extra_dict <> None: this_list.update(extra_dict)
        if code == 'mcmax':
            keywords = set(this_list.keys()+modellist.keys())
            if 'dust_species' in keywords:
                keywords.remove('dust_species')
        else:
            keywords = getattr(self,code + '_keywords')
        if code == 'mline' and ignoreAbun:
            keywords = [key 
                        for key in keywords 
                        if key not in ['ABUN_MOLEC','ABUN_MOLEC_RINNER',\
                                       'ABUN_MOLEC_RE','RMAX_MOLEC']]
        for keyword in keywords:
            try:
                try:
                    try:
                        val = float(this_list[keyword])
                    except TypeError:
                        raise ValueError
                    delta = not val and 1e-10 or 0.001*val
                    model_bool_list.append\
                        (val - delta < float(modellist[keyword]) < val + delta)
                except ValueError:
                    model_bool_list.append\
                        (this_list[keyword]==modellist[keyword])
            except KeyError:
                if keyword not in this_list.keys() \
                        and keyword not in modellist.keys():
                    model_bool_list.append(True)
                else:
                    model_bool_list.append(False)
        if False not in model_bool_list: 
            return True
        else: 
            return False