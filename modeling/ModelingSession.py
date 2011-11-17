# -*- coding: utf-8 -*-

"""
Interface for creating modeling environments.

Author: R. Lombaert

"""

import os

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
                    