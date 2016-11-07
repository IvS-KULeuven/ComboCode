# -*- coding: utf-8 -*-

"""
Examination of the Chemistry analysis routine output.

Author: M. Van de Sande

"""

import os
import scipy
from scipy import argmin,ones
from scipy import array
from scipy import sum
import operator
import types
import numpy as np

import cc.path
from cc.tools.io import DataIO
#from cc.modeling.tools import CodeIO
from cc.statistics.Statistics import Statistics
import matplotlib.pyplot as plt
from cc.modeling.objects import Transition
from cc.statistics import BasicStats    
from time import gmtime

class ChemStats(object):
    
    """
    Environment with several tools to perform statistics on the peak and 
    integrated intensities of resolved molecular lines.
    
    """
        
    def __init__(self,star_name,path_code='OutputClumpy',\
                 star_grid=[],molecules = []):        
        
        """ 
        
        @param star_name: Star name from Star.dat
        @type star_name: string
        
        @keyword code: the code used for producing your output 
        
                       (default: 'GASTRoNOoM')
        @type code: string
        @keyword path_code: Output folder in the code's home folder
                       
                            (default: 'codeSep2010')
        @type path_code: string
        @keyword molecules: Molecules for which the analysis output needs to be
                            retrieved
                            (default: [])
        @type molecules: list
                
        """
        
        #super(ChemStats,self).__init__(star_name=star_name,\
                                       #code=code,path_code=path_code)

        self.path_code = path_code
        self.molecules = molecules
        
        #-- Keep track of each Star() model that is valid.
        self.star_grid = dict()
        self.star_grid = [star for star in star_grid
                          if star['LAST_CHEMISTRY_MODEL']
                          and star['PERFORM_ROUTINE']==1]
        
        self.molecules = molecules
        #-- Convenience path
        cc.path.cout = os.path.join(cc.path.chemistry,self.path_code)
        
        
    def analyseMolecule(self, star, molec, print_main=1, print_sec=0):
        
        analyse = self.readAnalyse(star)
        
        if print_main == 1:
            self.printMainPathways(analyse,molec)
        
        if print_sec == 1:
            sec = self.getSecondaryPathways(analyse,molec,print_sec=1)
        
        
        
    def readAnalyse(self, star):
        
        '''
        Reads the output of the analyse routine performed.
        
        @keyword star: Star object containing the chemistry model
        @type star: Star()
        
        @return: Dictionary containing the analyse output for every
                         molecule, eg analyse['SiO']. The analyse radius is 
                         also included (analyse['radius']).
        @rtype: dict()
        
        '''
        
        #- Read in analyse.out
        filename = os.path.join(cc.path.cout,'models',\
                                star['LAST_CHEMISTRY_MODEL'],'analyse.out')
        data = DataIO.readFile(filename)
        data = [x.split() for x in data]
        
        #- Initialise output dictionary
        analyse = dict()
        
        #- Determine and save routine radius
        radius = float(data[0][-2])
        index = np.where([len(i)==2 for i in data])[0]
        analyse['radius'] = radius
        
        #- Save analyse routine output in dictionary per molecule
        for ii in range(len(index[:-1])):
            name = data[index[ii]][1][1:]
            analyse[name] = dict()
            analyse[name]['MAIN'] = data[index[ii]+1:index[ii+1]-1]
            analyse[name]['DRATE'] = data[index[ii+1]-1][1]
            analyse[name]['PRATE'] = data[index[ii+1]-1][3]
        
        return analyse

    
    def getSecondaryPathways(self,analyse,molec,print_sec=1,\
                             destruction=1,formation=1):
        
        others = analyse.keys()
        others.remove('radius')
        others.remove(molec)
        
        sec = dict()
        
        for x in others:
            sec[x] = dict()
            if destruction == 1 and formation == 1:
                dummy = [st for st in analyse[x]['MAIN'] if molec in st]
            elif destruction == 1 and formation == 0:
                dummy = [st for st in analyse[x]['MAIN'] if molec in st[:3]]
            elif destruction == 0 and formation == 1:
                dummy = [st for st in analyse[x]['MAIN'] if molec in st[3:]]
            if dummy != []:
                sec[x]['MAIN'] = dummy
                sec[x]['PRATE'] = analyse[x]['PRATE']
                sec[x]['DRATE'] = analyse[x]['DRATE']
            else:
                sec.pop(x)
        
        if print_sec == 1:
            print '***   Minor pathways at radius '+str(analyse['radius'])+', '+molec+'  ***'         
            for x in sec.keys():
                print '\t ***      '+x+'       ***'
                for st in sec[x]['MAIN']:
                    print '\t'.join(st)
                print '\t '+'DRATE = '+sec[x]['DRATE']+\
                    '\t'+'PRATE = '+sec[x]['PRATE']+'\n'

        return sec
        

    def printMainPathways(self,analyse,molec):
        
        print '***   Main pathways at radius '+str(analyse['radius'])+', '+molec+'  ***'         
        for m in analyse[molec]['MAIN']:
            print '\t'.join(m)
        print '\t '+'DRATE = '+analyse[molec]['DRATE']+\
                '\t'+'PRATE = '+analyse[molec]['PRATE']+'\n'


    def printSecondaryPathways(self,analyse,molec,\
                             destruction=1,formation=1):
        
        others = analyse.keys()
        others.remove('radius')
        others.remove(molec)
        
        sec = dict()
        
        for x in others:
            sec[x] = dict()
            if destruction == 1 and formation == 1:
                dummy = [st for st in analyse[x]['MAIN'] if molec in st]
            elif destruction == 1 and formation == 0:
                dummy = [st for st in analyse[x]['MAIN'] if molec in st[:3]]
            elif destruction == 0 and formation == 1:
                dummy = [st for st in analyse[x]['MAIN'] if molec in st[3:]]
            if dummy != []:
                sec[x]['MAIN'] = dummy
                sec[x]['PRATE'] = analyse[x]['PRATE']
                sec[x]['DRATE'] = analyse[x]['DRATE']
            else:
                sec.pop(x)
        
        print '***   Minor pathways at radius '+str(analyse['radius'])+', '+molec+'  ***'         
        for x in sec.keys():
            print '\t ***      '+x+'       ***'
            for st in sec[x]['MAIN']:
                print '\t'.join(st)
            print '\t '+'DRATE = '+sec[x]['DRATE']+\
                '\t'+'PRATE = '+sec[x]['PRATE']+'\n'

        
    #def printMainPathways(self,star):
        
        ##for star in self.star_grid:
        #analyse = self.readAnalyse(star)
    
        #print '\n'
        #print '***************************************************'
        #print '***   Main formation and destruction pathways   ***'
        #print '*** \t \t at radius '+str(analyse['radius'])+' \t  \t***'
        #print '***************************************************'
        
        #for molec in self.molecules:
            #print '\t ***      '+molec+'       ***'
            #for m in analyse[molec]['MAIN']:
                #print '\t'.join(m)
            #print '\t '+'DRATE = '+analyse[molec]['DRATE']+\
                    #'\t'+'PRATE = '+analyse[molec]['PRATE']+'\n'


        
        
        
        
        















