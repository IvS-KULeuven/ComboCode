# -*- coding: utf-8 -*-

"""
A class for reading and managing Mline output.

Author: M. van de Sande and R. Lombaert

"""

import os
import numpy as np

from cc.tools.io.Reader import Reader
from cc.tools.io import DataIO



class MlineReader(PopReader,MolReader):
    
    '''
    A Reader for Mline output files. 
    
    This class inherits from both PopReader and MolReader. 
    
    NYI - Work in progress.
    
    '''
    
    def __init__(self,fn):
        
        '''
        Creating an Mline object ready for reading Mline output.
        
        Reading mline output through full filename, including filepath.
        
        The mlinefile number is given as *, which is then replaced as needed.
        
        @param filename: The mline filename, including filepath. The mline
                         file number is given by *.
        @type filename: string        
        
        
        '''

        self.molecule = os.path.splitext(fn)[0].split('_')[-1]
        fn = fn.replace('ml1','ml*').replace('ml2','ml*').replace('ml3','ml*')
        super(MlineReader, self).__init__(fn=fn)
        self.molecule = os.path.splitext(self.fn)[0].split('_')[-1]
        self.readML1()
        self.readMlineParameters()
        self.readML3()



    def readMoleculeSettings(self):
    
        '''
        Read molecule spectroscopy settings from the ml1 file.
        
        This includes ny_up, ny_low, nline, n_impact and n_impact_extra.
        
        '''
        
        fn = self.fn.replace('ml*','ml1')
        dd = DataIO.readDict(filename=fn,start_index=8,end_index=38,\
                             convert_floats=1)
        self.contents['pars'] = dd
        
        trans = DataIO.readCols(filename=fn,)
        data[39:39+dd['NLINE']]
        trans = [t.split(' ') for t in trans]
        for ii in range(NLINE):
            trans[ii] = [t for t in trans[ii] if t]
        transnumber = [int(trans[x][0]) for x in range(NLINE)]
        lup = [int(trans[x][1]) for x in range(NLINE)]
        llow = [int(trans[x][2]) for x in range(NLINE)]
        
        


    def readMlineParameters(self):
        
        '''
        Read the extra mline parameters from the mline log file written by CC.
        
        Only when this file is available!
        
        Nothing is read from the database. This method is standalone, similar to
        SphinxReader. 
        
        '''
        
        
    def readML3(self):
    
        '''
        Read the information from the ML3 file. 
        
        This includes:
            - Scattering integral
            - Source function
            - Line opacities
            - Number densities 
            
        '''
        
        fn = self.fn.replace()
        
    
    
    
    
    def readMline3(modeldir, model, molecule,NY,NLINE,N_IMPACT,N_FREQ):
    f = open(modeldir+'ml3'+model[:-1]+'_'+molecule+'.dat')
    data = f.read().split('\n')
    f.close()
    data = [d.strip() for d in data]

    #- Split data into the six main blocks
    delimiters = np.where(array(data) == '--------------------------------------------------          --------------------------------------------------                --------------------------------')[0]
    data = [d.split(' ') for d in data]
    blok1 = data[:delimiters[0]]
    blok2 = data[delimiters[0]+1:delimiters[1]]
    blok3 = data[delimiters[1]+1:delimiters[2]]
    blok4 = data[delimiters[2]+1:delimiters[3]]
    blok5 = data[delimiters[3]+1:delimiters[4]]
    blok6 = data[delimiters[4]+1:]
    #- Remove empty arrays 
    blok1 = [x for x in blok1 if x != ['']]
    blok2 = [x for x in blok2 if x != ['']]
    blok3 = [x for x in blok3 if x != ['']]
    blok4 = [x for x in blok4 if x != ['']]
    blok5 = [x for x in blok5 if x != ['']]
    blok6 = [x for x in blok6 if x != ['']]
    
    ###- Read in scattering integral SI = blok1
    #- Convert list[list[string]] to list[string]
    blep = [item for sublist in blok1 for item in sublist]

    ##- Negative numbers are joined to the previous number,
    ##- so sometimes string = several numbers.
    #- Check for problems 
    prob = [i for i,x in enumerate(blep) if len(x) > 14]
    for i in range(len(blep)):
        if i in prob:
            x = blep[i]
            x = x.split('E')
            y = []
            for j in range(len(x)-1):
                if j == 0:
                    y.append(x[j]+'E'+x[j+1][:3])
                else:
                    y.append(x[j][3:]+'E'+x[j+1][:3])
            print y
            #- String of numbers stuck together 
            #- list of strings
            blep[i] = y
        else:
            #- Convert ok strings also to list[string]
            blep[i] = [blep[i]]
    #- Convert list[list[string]] to list[string]
    blok1 = [item for sublist in blep for item in sublist]
    blok1 = [float(x) for x in blok1 if x]

    SI = []
    for i in range(N_IMPACT):
        SI.append(blok1[i*NLINE:(i+1)*NLINE])



    ###- Read in source function SF = blok2
    ##- Same reasoning as before
    blep = [item for sublist in blok2 for item in sublist]
    prob = [i for i,x in enumerate(blep) if len(x) > 14]
    for i in range(len(blep)):
        if i in prob:
            x = blep[i]
            x = x.split('E')
            y = []
            for j in range(len(x)-1):
                if j == 0:
                    y.append(x[j]+'E'+x[j+1][:3])
                else:
                    y.append(x[j][3:]+'E'+x[j+1][:3])
            blep[i] = y
        else:
            blep[i] = [blep[i]]
    blok2 = [item for sublist in blep for item in sublist]
    blok2 = [float(x) for x in blok2 if x]

    SF = []
    for i in range(N_IMPACT):
        SF.append(blok2[i*NLINE:(i+1)*NLINE])



    ###- Read in line opacities LO = blok3 
    ##- Same reasoning as before
    blep = [item for sublist in blok3 for item in sublist]
    prob = [i for i,x in enumerate(blep) if len(x) > 14]
    for i in range(len(blep)):
        if i in prob:
            x = blep[i]
            x = x.split('E')
            y = []
            for j in range(len(x)-1):
                if j == 0:
                    y.append(x[j]+'E'+x[j+1][:3])
                else:
                    y.append(x[j][3:]+'E'+x[j+1][:3])
            blep[i] = y
        else:
            blep[i] = [blep[i]]
    blok3 = [item for sublist in blep for item in sublist]
    blok3 = [float(x) for x in blok3 if x]

    LO = []
    for i in range(N_IMPACT):
        LO.append(blok3[i*NLINE:(i+1)*NLINE])


    ###- Read in number densities ND = blok4
    blep = [item for sublist in blok4 for item in sublist]
    prob = [i for i,x in enumerate(blep) if len(x) > 14]
    for i in range(len(blep)):
        if i in prob:
            x = blep[i]
            x = x.split('E')
            y = []
            for j in range(len(x)-1):
                if j == 0:
                    y.append(x[j]+'E'+x[j+1][:3])
                else:
                    y.append(x[j][3:]+'E'+x[j+1][:3])
            blep[i] = y
        else:
            blep[i] = [blep[i]]
    blok4 = [item for sublist in blep for item in sublist]
    blok4 = [float(x) for x in blok4 if x]

    ND = []
    for i in range(N_IMPACT):
        ND.append(blok4[i*NY:(i+1)*NY])

    
    return SI, SF, LO, ND
    
