# -*- coding: utf-8 -*-

"""
Toolbox for creating and maintaining a radio-data database. 

Author: R. Lombaert

"""

import os, pyfits
from glob import glob

import cc.path
from cc.tools.io import DataIO
from cc.tools.io.Database import Database
from cc.data import LPTools


class Radio(Database):
    
    """
    Environment with several tools to load and save a database for Radio data.
    
    """
        
    def __init__(self,auto_parse=0,db_name='radio_data.db'):        
        
        """ 
        Initializing an instance of Radio.
        
        Inherits from the Database class and adds a bit of functionality for 
        the particular case of radio data.
        
        Initializing this object returns a Database() class that can contain 
        radio data in the given path. If no db is already present, a new one is
        created from scratch at given location.
        
        The name of the database is radio_data.db by default. Can be changed if 
        needed, but ComboCode does not work with this other name.
        
        If nothing more is requested, a new database is essentially empty. 
        
        The option to parse the content of the target folder is available: This 
        will add entries to the database for the stars in the folder and the 
        'simple' molecules and transitions. 
                
        Once a db is created, you can perform several methods on the database 
        specific to radio data management.
        
        Examples:
        For data located in /home/robinl/Data/Molecular/
        >>> from cc.data import Radio
        >>> db = Radio.Radio('/home/robinl/Data/Molecular/',auto_parse=1)
        >>> db.sync()

        This will create a database, automatically search your data folder for 
        data files that the script can recognize, and save them in the 
        database. The final command ( .sync() ) will then save the database to
        your hard disk in the same folder, under 'radio_data.db'. Later when 
        you run Radio.Radio('/home/robinl/Data/Molecular/') again, it will load
        that same database that still contains all those datafiles and 
        TRANSITION references.

        The database is structured per star_name and per transition definition.
        So for instance, you will find the file whya_co32_APEX.fits in a dict 
        under
        db['whya']['TRANSITION=12C16O 0 3 0 0 0 2 0 0 APEX 0.0']
        where you will find
        dict([('whya_co32_APEX.fits',None)])
        
        The None refers to the value of the dictionary entry for that file. 
        This is replaced by the fit results of the line when 
        >>> db.fitLP(filename='whya_co32_Maercker_new_JCMT.fits')
        is ran. The method redoes the fit by default, regardless of there being
        an entry for it already. Any required fitting parameters can be passed
        along to the function (see LPTools.fitLP()). CC loads fit results from
        the database.

        If multiple data files are available for the same star, the same 
        transition and the same telescope, the multiple data files will also be
        contained in this dict. Note that the transition definition is quite 
        specific: 1 space between entries, and 11 entries total. It is 
        identical to the input in your combocode inputfile (excluding the final
        number n_quad, which is model input and has nothing to do with your 
        data).

        Lastly, you can do multiple things with this database.
        >>> db.addData(star_name='whya',trans='TRANSITION=12C16O 0 3 0 0 0 2 0 0 JCMT 0.0',\
                       filename='whya_co32_Maercker_new_JCMT.fits')
        will add a new entry to the database if it is not there. This is useful
        for those molecules/filenames which have not been automatically found 
        by the auto_parse (for instance SO2 data files do not stick to the file
        naming convention in my folder (do ls *so2* in my data folder), because
        with this method you can still add them whichever way you want.

        >>> db.removeData(star_name='whya',trans='TRANSITION=12C16O 0 3 0 0 0 2 0 0 JCMT 0.0')
        removes a whole transition from the database. Alternatively

        >>> db.removeData(star_name='whya',filename='whya_co32_Maercker_new_JCMT.fits')
        removes a single filename from the database. (does not delete the file)

        And if you want to know which datafiles in your data folder are not yet
        in the database (so you know which you have to add manually):
        >>> db.checkFolder()
        prints them out for you.

        Remember, every time you add or remove data or do a parseFolder(), you 
        have to save the contents of your database in python to the hard disk 
        by doing:
        >>> db.sync()

        As long as you don't do this, the hard disk version of the datafile 
        will remain unchanged. 
        
        @param path: Path to the database excluding db name. This is the folder
                     where the radio data are located.
        @type path: string
        
        @keyword auto_parse: Automatically parse the target folder for data to be 
                             included in the database. 
                         
                             (default: 0)
        @type auto_parse: bool
        @keyword db_name: The name of the database. The default only is used by 
                          ComboCode.
                      
                          (default: radio_data.db)
        @type db_name: str      
        
        """
        
        db_path = os.path.join(cc.path.dradio,db_name)
        super(Radio,self).__init__(db_path=db_path)
                
        if auto_parse:
            self.parseFolder()
            
            

    def parseFolder(self):
        
        '''
        Parse the db folder and add filenames to recognized transitions. 
        
        Includes:
            - CO and its isotopologues
            - SiO and its isotopologues
            - SiS and its isotopologues
            - H2O and pH2O
            - CS
            - PO and PN
        
        Particularly excludes: (because no naming convention/too complex)
            - SO and SO2
            - H2O and pH2O isotopologues
            - H2CO
            - CN
            - HCN and its isotopologues
            - HCO+
        
        '''
        
        ggf = sorted(glob(self.folder+'/*.dat') + glob(self.folder+'/*.fits'))
        ggf = [os.path.split(ff)[1] for ff in ggf]
        ggf = [ff for ff in ggf if '_' in ff]
        ggf = [ff for ff in ggf if ff[0] != '_']
        fsplit = [ff.split('_') for ff in ggf]
        stars = [ff[0] for ff in fsplit]
        molec_trans = [(ff[1][0] == 'v' and len(ff[1]) == 2) and ff[2] or ff[1] 
                       for ff in fsplit]
        vibs = [(ff[1][0] == 'v' and len(ff[1]) == 2) and ff[1][1] or 0
               for ff in fsplit]
        telescopes = [os.path.splitext(ff[-1])[0] for ff in fsplit]
        
        defmolecs = DataIO.getInputData(keyword='TYPE_SHORT',\
                                        filename='Molecule.dat')
        defmolecs_short = DataIO.getInputData(keyword='NAME_SHORT',\
                                              filename='Molecule.dat')
        defspec_indices = DataIO.getInputData(keyword='SPEC_INDICES',\
                                              filename='Molecule.dat')
        for ff,s,mt,tel,vib in zip(ggf,stars,molec_trans,telescopes,vibs):
            this_dm,this_dms,this_dsi = None, None, None
            for dm,dms,dsi in zip(defmolecs,defmolecs_short,defspec_indices):
                if mt[:len(dms)] == dms:
                    this_dm = dm
                    this_dms = dms
                    this_dsi = dsi
                    break
            evil = ['so','so2','h217o','h218o','ph217o','ph218o',\
                    'cn','hcn','h13cn','hco+']
            if this_dm == None or this_dms in evil:
                continue
            traw = mt.replace(this_dms,'',1)
            if this_dsi == 2:
                if len(traw) == 4:
                    trans = 'TRANSITION=%s %s %s %s 0 %s %s %s 0 %s 0.0'\
                            %(this_dm,vib,traw[0],traw[1],\
                              vib,traw[2],traw[3],tel)
                else:
                    continue
            elif this_dsi == 0:
                if len(traw) == 4:
                    trans = 'TRANSITION=%s %s %s 0 0 %s %s 0 0 %s 0.0'\
                            %(this_dm,vib,traw[0:2],vib,traw[2:4],tel)
                elif len(traw) == 3:
                    trans = 'TRANSITION=%s %s %s 0 0 %s %s 0 0 %s 0.0'\
                            %(this_dm,vib,traw[0:2],vib,traw[2],tel)
                elif len(traw) == 2:
                    trans = 'TRANSITION=%s %s %s 0 0 %s %s 0 0 %s 0.0'\
                            %(this_dm,vib,traw[0],vib,traw[1],tel)
                else:
                    continue
            elif this_dsi == 1:
                if len(traw) == 6:
                    trans = 'TRANSITION=%s %s %s %s %s %s %s %s %s %s 0.0'\
                            %(this_dm,vib,traw[0],traw[1],traw[2],\
                              vib,traw[3],traw[4],traw[5],tel)
                else:
                    continue
            else:
                continue
            self.addData(star_name=s,trans=trans,filename=ff)
            


    def addStar(self,star_name):
        
        '''
        Add a new star to the database. A check is ran to see if the requested 
        star is known in ~/ComboCode/usr/Star.dat. 
        
        @param star_name: The name of the star (from Star.dat)
        @type star_name: str
               
        '''
        
        if star_name in self.keys():
            return
        known_stars = DataIO.getInputData(keyword='STAR_NAME')
        if star_name not in known_stars:
            print('Requested star %s is unknown in %s/Star.dat.'\
                  %(star_name,cc.path.usr))
            return
        self[star_name] = dict()



    def addData(self,star_name,filename,trans):
        
        '''
        Add a new file to the database with given transition definition. If the 
        transition is already present in the database for this particular star, the
        filename is added to the dictionary for that transition.
        
        A single transition for a star can thus have multiple filenames associated
        with it!
        
        Note that this method does NOT automatically sync (ie save changes to the 
        hard disk) the database. That must be done through an additional flag to 
        avoid excess overhead. 
        
        The transition definition is checked for correctness: single spaces between
        entries in the defintion, and a total of 11 entries: 1 molecule (shorthand)
        8 quantum numbers, 1 offset 
        
        The fit is not done here. Instead the filename key is added to the 
        transition's dictionary with value None.
        
        @param star_name: The name of the star for which to add the data.
        @type star_name: str
        @param filename: The filename of the radio data to be added. Must be in the 
                        db folder! Only the filename is used. Any folder path is 
                        cut.
        @type filename: str
        @param trans: The transition definition. Must be in the correct format!
        @type trans: str
        
        '''
        
        filename = os.path.split(filename)[1]
        if not self.has_key(star_name):
            self.addStar(star_name)
        
        entries = trans.split()
        trans = ' '.join(entries)
        
        defmolecs = DataIO.getInputData(keyword='TYPE_SHORT',\
                                        filename='Molecule.dat')
        this_molec = entries[0].replace('TRANSITION=','',1)
        if this_molec not in defmolecs:
            print('Molecule %s unrecognized.'%this_molec)
            return
        
        if len(entries) != 11: 
            print('Incorrect transition definition:')
            print(trans)
            return
        
        if self[star_name].has_key(trans):
            if filename not in self[star_name][trans].keys():
                self[star_name][trans][filename] = None
                self.addChangedKey(star_name)
        else:
            self[star_name][trans] = dict([(filename,None)])
            self.addChangedKey(star_name)
            


    def removeData(self,star_name,filename='',trans=''):
        
        '''
        Remove filenames or transition definitions from the database for a given
        star_name. Either filename or trans must be given. 
        
        @param star_name: Name of the star for which data or trans is removed
        @type star_name: str
        
        @keyword filename: data filename to be removed. If empty string, a 
                           transition definition must be given (ie remove 
                           everything for a given transition definition)
                           
                           (default: '')
        @type filename: str
        @keyword trans: Transition definition to be removed. Irrelevant if a 
                        filename is also given (latter takes precedence)
                        
                        (default: '')                        
        @type trans: str
        
        '''
        
        if not self.has_key(star_name):
            print('Star %s not found in database.'%star_name)
            return
        
        #-- If a filename is given, remove that, regardless of a given transition.
        if filename:
            trans = ''
            filename = os.path.split(filename)[1]
        elif not trans: 
            print('No filename or transition given. No change is made to the'+\
                  'radio database')
            return
        
        if filename: 
            for k,v in self[star_name].items():
                if filename in v.keys():
                    del self[star_name][k][filename]
                    self.addChangedKey(star_name)
                    break
        else:
            if self[star_name].has_key(trans):
                del self[star_name][trans]
                self.addChangedKey(star_name)
        


    def checkFolder(self):
        
        '''
        Check which radio datafiles are not included in the db present in its 
        folder.

        '''
        
        ggf = sorted(glob(self.folder+'/*.dat') + glob(self.folder+'/*.fits'))
        ggf = [os.path.split(ff)[1] for ff in ggf]
        ggf = [ff for ff in ggf if ff[0] != '_']
        
        dbfiles = []
        for s in self.keys():
            for v in self[s].values():
                dbfiles.extend(v)

        print('These files have not been found in the database:')
        for ff in ggf:
            if ff not in dbfiles: 
                print(ff)
        
    
    def convertOldDb(self):
        
        '''
        Convert an existing Radio database from the old format:
        db[star][trans] = list
        to the new format:
        db[trans][trans][filename] = fit dict
        
        '''
        
        for ss in self.keys():
            for tt in self[ss].keys():
                nd = dict()
                for ff in self[ss][tt]:
                    nd[ff] = None
                self[ss][tt] = nd
            self.addChangedKey(ss)
    
    
    def filterContents(self,star_name=''):
        
        '''
        Filter out the fit results from the contents of the database, eg for
        printing purposes. 
        
        If star_name is given, only the star will be printed. Otherwise the 
        full database is printed. 
        
        @keyword star_name: The requested star. Default is to print all of them
        
                            (default: '')
        @type star_name: str
        
        @return: The dictionary with just the relevant contents (trans and file
                 names)
        @rtype: dict
        
        ''' 
        
        def __removeFits(self,ss):
            
            '''
            Helper method to retrieve only the filenames for every transition
            of given star. 
            
            '''
            
            fdd = dict()
            for tt in self[ss].keys():
                fdd[tt] = [ff for ff in sorted(self[ss][tt].keys())]
            return fdd
         
         
        pdd = dict()
        if star_name: 
            pdd[star_name] = __removeFits(self,star_name)
        else:
            for ss in self.keys():
                pdd[ss] = __removeFits(self,ss)
        
        return pdd
        
        
    
    
    def fitLP(self,star_name='',filename='',trans='',**kwargs):
        
        '''
        Fit the data line profiles with a soft parabola or a Gaussian according
        to LPTools.fitLP() (see that method for further details).
        
        Input keywords require one of these:
            - no keys: all lines in db are fitted
            - star_name: All lines for star are fitted
            - star_name, trans: all lines for transition of star are fitted
            - star_name, filename: only this filename is fitted
        
        The fit is redone by default, regardless of there being an entry in the
        db already. 
        
        Note that this method does NOT automatically sync (ie save changes to 
        the hard disk) the database. That must be done through an additional 
        flag to avoid excess overhead. 
        
        @keyword star_name: The name of the star for which to add the data. If
                            not given, all files in db are fitted.
                            
                            (default: '')
        @type star_name: str
        @keyword filename: The filename of the radio data to be fitted. Must be 
                           in the db folder! Only the filename is used. Any 
                           folder path is cut. 

                           (default: '')
        @type filename: str
        @keyword trans: The transition definition. Must be in the correct 
                        format!
        
                        (default: '')
        @type trans: str
        @keyword kwargs: Any additional keywords that are passed on to 
                         LPTools.fitLP()
        @type kwargs: dict
                         
        '''
        
        def __fitLP(self,ss,tt,ff):
            
            '''
            Helper method giving star_name, transition, and filename for the 
            fitter to do and then remember to save when sync() is ran.
            
            Accessed only by Radio().fitLP().
            
            '''
            
            fn = os.path.join(self.folder,ff)
            try:
                fitr = LPTools.fitLP(filename=fn,**kwargs)
            except ValueError:
                fitr = None
            self[ss][tt][ff] = fitr
            self.addChangedKey(ss)
            
            
        if star_name and not self.has_key(star_name):
            print 'Star not found.'
            return
        
        #-- No star_name given, so run through all stars, transitions and files
        if not star_name:
            for ss in self.keys():
                for tt in self[ss].keys():
                    for ff in self[ss][tt].keys():
                        __fitLP(self,ss,tt,ff)
        
        #-- star_name given. If trans is given, run through all its filenames 
        elif trans:
            if trans not in self[star_name].keys():
                print 'Transition not found.'
                return
            for ff in self[star_name][trans].keys():
                __fitLP(self,star_name,trans,ff)
            
        #-- star_name given. If trans is not given, but filename is, fit it. 
        elif filename: 
            filename = os.path.split(filename)[1]
            for k,v in self[star_name].items():
                if filename in v.keys():
                    trans = k
                    break
            if not trans: 
                print 'Filename not found.'
                return
            __fitLP(self,star_name,trans,filename)
        
        #-- star_name given. No trans/filename given. Fit everything for star
        else:
            for tt in self[star_name].keys():
                for ff in self[star_name][tt].keys():
                    __fitLP(self,star_name,tt,ff)
        