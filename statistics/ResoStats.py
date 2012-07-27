# -*- coding: utf-8 -*-

"""
Performing statistics on the peak and integrated intensity of resolved 
molecular lines.

Author: R. Lombaert

"""

import os
import scipy
from scipy import argmin
from scipy import array as arr
import operator

from cc.tools.io import DataIO
from cc.statistics.Statistics import Statistics



class ResoStats(Statistics):
    
    """
    Environment with several tools to perform statistics on the peak and 
    integrated intensities of resolved molecular lines.
    
    """
        
    def __init__(self,star_name,path_code='codeSep2010',stat_print=0,\
                 path_combocode=os.path.join(os.path.expanduser('~'),\
                                             'ComboCode')):        
        
        """ 
        Initializing an instance of IntIntStats.
        
        @param star_name: Star name from Star.dat
        @type star_name: string
        
        @keyword print_stats: Print statistics for individual transitions and 
                              models.
                              
                              (default: 0)
        @type print_stats: bool
        @keyword path_code: Output folder in the code's home folder
                       
                            (default: 'codeSep2010')
        @type path_code: string
        @keyword path_combocode: CC home folder
        
                                 (default: '~/ComboCode/')
        @type path_combocode: string
        
        """
        
        super(ResoStats,self).__init__(star_name=star_name,\
                                       path_combocode=path_combocode,\
                                       code='GASTRoNOoM',path_code=path_code)
        #- Dicts keeping integrated and peak intensities of datafiles. 
        #- key: the sample Transition(), 
        #- value: list of values for first datafile included
        #- Note that for every model the data value is recalculated, based on 
        #- the best vlsr guess which depends on the model.
        #- As a result, every element in the value list corresponds to the 
        #- Star() object with the same index in self.star_grid.
        self.dinttmb = dict()
        self.dpeaktmb = dict()
        #- Dicts keeping integrated and peak intensities of sphinx models, as 
        #- well as the loglikelihood measure of fitness between model and data
        #- key: the sample Transition(), 
        #- value: list of values for every model included
        #- Every element in the value list corresponds to the Star() object 
        #- with the same index in self.star_grid.
        self.minttmb = dict()
        self.mpeaktmb = dict()
        self.loglikelihood = dict()
        self.stat_print = stat_print
        
    
    
    def setInstrument(self,sample_transitions):
       
        '''
        Set and read the data objects for this statistics module. 
        
        In this case, a list of sample transitions in which the data will be 
        read is the only input required.
    
        @param sample_transitions: Sample transitions used as reference for 
                                   the data files. 
        @type sample_transitions: list[Transition()]
        
        '''
        
        super(ResoStats,self).setInstrument('FREQ_RESO',sample_transitions)
            
            
            
    def setModels(self,star_grid):
       
        '''
        Set and read the models for this statistics module. 
        
        In this case, a list of Star() objects in which the sphinx files will 
        be read is the only input required.
    
        @param star_grid: The parameter sets, if not given: model ids needed
        
                          (default: [])
        @type star_grid: list[Star()]
        
        '''
        
        super(ResoStats,self).setModels('FREQ_RESO',star_grid)        
        
        
        
    def getIntensities(self):
        
        """
        The data intensities are stored in the dictionary 
        self.dintint/self.dpeakint, the model intensities in 
        self.mintint/self.mpeakint. 
        
        They keys in the dictionaries are the Transition to which the sphinx 
        or datafiles belong. 
        
        The model values are lists in accordance with self.star_grid, the data
        values are single floats for the first dataset for the transition.

        """
        
        tnodata = [t for t in self.instruments['FREQ_RESO'] 
                     if not t.lpdata and t.getModelId()]
        tnomodels = [t for t in self.instruments['FREQ_RESO'] 
                       if not t.getModelId() and t.lpdata]
        translist = [t for t in self.instruments['FREQ_RESO']
                       if t.lpdata and t.getModelId()]
                   
        #- For every sample transition (st), collect the equivalent transitions
        #- in the model grid. Then retrieve all integrated and peak tmb values,
        #- for both data and model. Data as well due to the use of the best 
        #- vlsr results, which is based on sphinx. In theory these should be identical!
        for st in translist:
            mtrans = [star.getTransition(st)
                      for star in self.star_grid]
            for mt in mtrans: 
                mt.setData(st)
                mt.getBestVlsr()
            self.minttmb[st] = arr([mt.getIntTmbSphinx() 
                                    for mt in mtrans])
            self.mpeaktmb[st] = arr([mt.getPeakTmbSphinx() 
                                     for mt in mtrans])
            self.dinttmb[st] = arr([mt.getIntTmbData() 
                                    for mt in mtrans])
            self.dpeaktmb[st] = arr([mt.getPeakTmbData() 
                                     for mt in mtrans])
            self.loglikelihood[st] = arr([mt.getLoglikelihood() 
                                          for mt in mtrans])
            ratioint = self.minttmb[st]/self.dinttmb[st]
            ratiopeak = self.mpeaktmb[st]/self.dpeaktmb[st]
            
            if self.stat_print:
                print '*******************************************************'
                print 'Statistics for %s:'%str(st)
                print '-------------------------------------'
                print 'Data intensities:'
                print 'MEAN (int): %f K, STD (int): %f K'\
                      %(mean(self.dinttmb[st]), std(self.dinttmb[st]))
                print 'MEAN (peak): %f K, STD (peak): %f K'\
                      %(mean(self.dpeaktmb[st]), std(self.dpeaktmb[st])) 
                print '-------------------------------------'
                print 'Model/Data intensities [integrated --- peak --- lll]:'
                lines = ['- %s: \t %.1f \t---\t %.1f \t---\t %.4f\
                         '%(str(st),ri,rp,lll)
                        for ri,rp,lll in zip(ratioint,ratiopeak,\
                                             self.loglikelihood[st])]
                print '\n'.join(lines)