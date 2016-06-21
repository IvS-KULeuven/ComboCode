# -*- coding: utf-8 -*-

import os

#-- Method only accessible when importing cc.data. 
#   Should not be called otherwise. Method is here to avoid unnecessary globals
def __readPaths():

    '''
    Read the relevant paths from the usr/Path.dat file.
    
    Default values typically return an expected location for model related 
    folders, and an empty string for data locations (in which case no data are
    used).
    
    @return: The pathnames with associated path are returned.
    @rtype: dict
    
    '''
    
    dd = dict([('gdata','GASTRoNOoM/src/data'),('gastronoom','GASTRoNOoM'),\
               ('chemistry','Chemistry'),\
               ('mcmax','MCMax'),('mobs','MCMax/Observation_Files'),\
               ('mopac','MCMax/Opacities'),('data',''),('dradio',''),\
               ('dpacs',''),('dspire',''),('dsed',''),('dphot',''),\
               ('dcflux',''),('ll','LineLists'),('ivsdata',''),('atm',''),\
               ('starf',''),('densf','')])
    filename = os.path.join(usr,'Path.dat')
    FILE = open(filename,'r')
    lines = FILE.readlines()
    FILE.close()
    lines = [line.strip('\n').strip() for line in lines if line.strip('\n')]
    dd.update(dict([line.split('=') for line in lines if line[0] != '#']))
    for k,v in dd.items():
        if not v:
            dd[k] = ''
        elif not v[0] == '/':
            dd[k] = os.path.join(os.path.expanduser('~'),v)
    return dd



#-- Define the home, usr and aux paths for ComboCode
home = os.path.dirname(__file__).replace('cc/path','')
usr = os.path.join(home,'usr')
aux = os.path.join(home,'aux')

#-- Insert the path names and values taken from usr/Path.dat into the globals 
#   of the module's namespace. These paths are only available through, e.g., 
#   cc.paths.pathname
globals().update(__readPaths())

#-- Note that this module can be used to set additional paths in any given 
#   python session. The paths are then always remembered until they are changed
#   within that session. e.g. mout and gout give the mcmax and gastronoom 
#   output folders, but they are only created when ComboCode is ran. 
