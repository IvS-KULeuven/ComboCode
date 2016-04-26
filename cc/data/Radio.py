# -*- coding: utf-8 -*-

"""
Toolbox for creating and maintaining a radio-data database.

Author: R. Lombaert

"""

import os, re
from glob import glob

import cc.path
from cc.tools.io import DataIO
from cc.tools.io.Database import Database
from cc.data import LPTools


class Radio(Database):

    """
    Environment with several tools to load and save a database for Radio data.

    """

    def __init__(self,db_name='radio_data.db',db_path=None):

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
        For data located in dradio in Path.dat
        >>> from cc.data import Radio
        >>> db = Radio.Radio()
        >>> db.sync()

        This will create a database, automatically search your data folder for
        data files that the script can recognize, and save them in the
        database. The final command ( .sync() ) will then save the database to
        your hard disk in the same folder, under 'radio_data.db'. Later when
        you run Radio.Radio() again, it will load that same database that still
        contains all those datafiles and TRANSITION references.

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
        >>> db.addData(star_name='whya',\
                       trans='TRANSITION=12C16O 0 3 0 0 0 2 0 0 JCMT 0.0',\
                       filename='whya_co32_Maercker_new_JCMT.fits')
        will add a new entry to the database if it is not there. This is useful
        for those molecules/filenames which have not been automatically found
        by the parseFolder method, because
        with this method you can still add them whichever way you want.

        >>> db.removeData(star_name='whya',\
                          trans='TRANSITION=12C16O 0 3 0 0 0 2 0 0 JCMT 0.0')
        removes a whole transition from the database. Alternatively

        >>> db.removeData(star_name='whya',\
                          filename='whya_co32_Maercker_new_JCMT.fits')
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

        @keyword db_name: The name of the database. The default only is used by
                          ComboCode.

                          (default: radio_data.db)
        @type db_name: str
        @keyword db_path: The path to the db if not the default dradio in 
                          Path.dat
        
                          (default: None)
        @type db_path: str
        

        """
        
        if db_path is None:
            fn = os.path.join(cc.path.dradio,db_name)
        else:
            fn = os.path.join(db_path,db_name)
        super(Radio,self).__init__(db_path=fn)



    def parseFolder(self):

        '''
        Parse the db folder and add filenames to recognized transitions.

        Includes:
            - CO and its isotopologues
            - SiO and its isotopologues
            - SiS and its isotopologues
            - H2O and pH2O, and their isotopologues
            - SO2 and SO
            - CS
            - PO and PN
            - H2CO

        Particularly excludes: (because no naming convention/too complex)
            - CN
            - HCN and its isotopologues
            - HCO+

        '''

        ggf = sorted(glob(self.folder+'/*.dat') + glob(self.folder+'/*.fits'))
        ggf = [os.path.split(ff)[1] for ff in ggf]
        ggf = [ff for ff in ggf if '_' in ff]
        ggf = [ff for ff in ggf if ff[0] != '_']
        fsplit = [ff.split('_') for ff in ggf]

        #-- Convention: star names come first
        stars = [ff.pop(0) for ff in fsplit]

        #-- Convention: telescope names are at end of filename before extension
        telescopes = [os.path.splitext(ff.pop(-1))[0] for ff in fsplit]

        #-- Convention: vibrational states != 0 in second place in filename.
        #   Vibrational states v1 are added automatically. Others are not, for
        #   now.
        vibs = [(ff[0][0] == 'v' and len(ff[0]) == 2) and ff.pop(0)[1] or '0'
               for ff in fsplit]

        #-- Convention: molecule names always first in the molecule tag
        molec_trans = ['_'.join(ff) for ff in fsplit]

        defmolecs = DataIO.getInputData(keyword='TYPE_SHORT',\
                                        filename='Molecule.dat')
        defmolecs_short = DataIO.getInputData(keyword='NAME_SHORT',\
                                              filename='Molecule.dat')
        defspec_indices = DataIO.getInputData(keyword='SPEC_INDICES',\
                                              filename='Molecule.dat')
        for ff,s,mt,tel,vib in zip(ggf,stars,molec_trans,telescopes,vibs):
            if vib not in ['1','0']:
                continue
            this_dm,this_dms,this_dsi = None, None, None
            for dm,dms,dsi in zip(defmolecs,defmolecs_short,defspec_indices):
                if mt[:len(dms)] == dms:
                    this_dm = dm
                    this_dms = dms
                    this_dsi = dsi
                    break

            #-- Don't bother with these molecules: not clear how they are used
            #   as input for GASTRoNOoM.
            #-- Note that the parser can already detect decimal quantum numbers
            #   as required by CN
            evil = ['cn','hcn','h13cn','hco+','oh']

            #-- If molecule is not found, or is an evil molecule, or is not of
            #   spectral index type 0, 1, or 2, then skip the entry.
            if this_dm == None or this_dms in evil or this_dsi not in [0,1,2]:
                continue

            mt = mt.replace(this_dms,'',1)
            
            #-- Create the regular expression for the filename convention
            #   Allows co32, h2o123132, h2o-1_2_3-1_3_2
            p = re.compile('\d{2,}|-[\d\.?\d_?]+-[\d\.?\d_?]+')            
            parsed = p.match(mt).group().rstrip('_').split('-')
            #-- Note that if a NoneType AttributeError occurs, it's probably 
            #   because the molecule has not been added to the convention list 
            #   yet, and is therefore considered evil. Added it to the evil list

            #-- if parsed contains one element, several options are possible
            if len(parsed) == 1:
                #-- dsi == 1 can only take 6 quantum numbers, 3 per level
                if len(parsed[0]) == 6 and this_dsi == 1:
                    tupper = parsed[0][0:3]
                    tlower = parsed[0][3:6]
                #-- dsi == 2 can only take 4 quantum numbers, 2 per level
                elif len(parsed[0]) == 4 and this_dsi == 2:
                    tupper = parsed[0][0:2]
                    tlower = parsed[0][2:4]
                #-- dsi == 0 can only take 2 quantum numbers, 1 per level
                #   several combinations are possible for multiple digit levels
                elif len(parsed[0]) == 6 and this_dsi == 0:
                    tupper = [parsed[0][0:3]]
                    tlower = [parsed[0][3:6]]
                elif len(parsed[0]) == 5 and this_dsi == 0:
                    tupper = [parsed[0][0:3]]
                    tlower = [parsed[0][3:5]]
                elif len(parsed[0]) == 4 and this_dsi == 0:
                    tupper = [parsed[0][0:2]]
                    tlower = [parsed[0][2:4]]
                elif len(parsed[0]) == 3 and this_dsi == 0:
                    tupper = [parsed[0][0:2]]
                    tlower = [parsed[0][2:3]]
                elif len(parsed[0]) == 2 and this_dsi == 0:
                    tupper = [parsed[0][0]]
                    tlower = [parsed[0][1]]
                #-- Any other combinations are not possible, and must be
                #   separated with '_', '-' in the filename
                else:
                    continue

            #-- If there are three elements, the second and  third give the
            #   upper and lower level respectively, each quantum number
            #   separated by '_'
            elif len(parsed) == 3:
                #-- parsed[0] is just an empty string
                tupper = parsed[1].split('_')
                tlower = parsed[2].split('_')
                #-- Double check if the correct amount of quantum numbers are
                #   given for each type
                if (this_dsi == 0 and len(tupper) != 1) or \
                   (this_dsi == 1 and len(tupper) != 3) or \
                   (this_dsi == 2 and len(tupper) != 2):
                   continue

            #-- If neither, continue without adding anything to the db
            else:
                continue

            if this_dsi == 2:
                trans = 'TRANSITION=%s %s %s %s 0 %s %s %s 0 %s 0.0'\
                        %(this_dm,vib,tupper[0],tupper[1],\
                          vib,tlower[0],tlower[1],tel)
            elif this_dsi == 1:
                trans = 'TRANSITION=%s 0.0'\
                        %(' '.join([this_dm,vib,tupper[0],tupper[1],tupper[2],\
                                    vib,tlower[0],tlower[1],tlower[2],tel]))
            elif this_dsi == 0:
                trans = 'TRANSITION=%s %s %s 0 0 %s %s 0 0 %s 0.0'\
                        %(this_dm,vib,tupper[0],vib,tlower[0],tel)

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

                    #-- Remove the transition if no filenames are associated
                    #   with it anymore, i.e. the dict is empty.
                    if not self[star_name][k]:
                        del self[star_name][k]
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




    def fitLP(self,star_name='',filename='',trans='',replace=0,**kwargs):

        '''
        Fit the data line profiles with a soft parabola or a Gaussian according
        to LPTools.fitLP() (see that method for further details).

        Input keywords require one of these:
            - no keys: all lines in db are fitted
            - star_name: All lines for star are fitted
            - star_name, trans: all lines for transition of star are fitted
            - star_name, filename: only this filename is fitted

        The fit is NOT redone by default, if there is an entry in db already. 
        You can force a replacement fit by turning replace on.

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
        @keyword replace: Replace existing fits with a new one. If off, only 
                          lines that have not been fitted yet are added. 
                          
                          (default: 0)
        @type replace: bool
        @keyword kwargs: Any additional keywords that are passed on to
                         LPTools.fitLP()
        @type kwargs: dict
        

        '''

        def __fitLP(self,ss,tt,ff,replace):

            '''
            Helper method giving star_name, transition, and filename for the
            fitter to do and then remember to save when sync() is ran.

            Accessed only by Radio().fitLP().

            '''
            
            if not replace and self[ss][tt][ff]: return
            fn = os.path.join(self.folder,ff)
            try:
                fitr = LPTools.fitLP(filename=fn,**kwargs)
            except ValueError:
                print 'Line profile fit in %s failed.'%ff
                fitr = None
            self[ss][tt][ff] = fitr
            self.addChangedKey(ss)


        if filename and not star_name:
            star_name = os.path.split(filename)[1].split('_')[0]

        if trans and not star_name:
            print 'Define star_name to continue.'
            return

        if star_name and not self.has_key(star_name):
            print 'Star not found.'
            return

        #-- No star_name given, so run through all stars, transitions and files
        if not star_name:
            for ss in self.keys():
                for tt in self[ss].keys():
                    for ff in self[ss][tt].keys():
                        __fitLP(self,ss,tt,ff,replace)

        #-- star_name given. If trans is given, run through all its filenames
        elif trans:
            if trans not in self[star_name].keys():
                print 'Transition not found.'
                return
            for ff in self[star_name][trans].keys():
                __fitLP(self,star_name,trans,ff,replace)

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
            __fitLP(self,star_name,trans,filename,replace)

        #-- star_name given. No trans/filename given. Fit everything for star
        else:
            for tt in self[star_name].keys():
                for ff in self[star_name][tt].keys():
                    __fitLP(self,star_name,tt,ff,replace)
