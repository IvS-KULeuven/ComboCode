# -*- coding: utf-8 -*-

"""
Module for reading and writing data files as well as parsing of information.

Author: R. Lombaert

"""

import os
import subprocess
from scipy import array,zeros
import types
from pylab import mlab


def getGastronoomOutput(filename,keyword='RADIUS',begin_index=0,\
                        return_array=0,key_index=0):
    
    """
    Search GASTRoNOoM output for relevant envelope information.

    @param filename: The filename of the relevant output GASTRoNOoM file
    @type filename: string
    
    @keyword keyword: the type of information required, always equal to one of 
                      the keywords present in the outputfiles of GASTRoNOoM
                      
                      (default: 'RADIUS')
    @type keyword: string
    @keyword begin_index: start looking for keyword at row with begin_index
                    
                          (default: 0)
    @type begin_index: int
    @keyword return_array: Return a scipy array rather than a python list
    
                           (default: 0)
    @type return_array: bool
    @keyword key_index: If 0 it is automatically determined, otherwise this is 
                        the column index
                        
                        (default: 0)
    @type key_index: int
    
    @return: The requested data from the GASTRoNOoM output
    @rtype: list/array
   
    """
    
    keyword = keyword.upper()
    data = readFile(filename,' ')
    data_col_1 = [d[0] for d in data]
    key_i = findString(begin_index,data_col_1)
    key_j = findFloat(key_i,data_col_1)
    if not key_index:
        keys = ' '.join([' '.join(d).replace('\n','') 
                         for d in data[key_i:key_j]]).split()
        key_index = [key[:len(keyword)].upper() for key in keys].index(keyword)
    #- Data never start on the first line
    #- Starting from 1st float, all floats into list, until EOF OR end of block
    data_i = key_j
    #- Data may end at EOF or before a new block of data (sphinx fi)
    data_j = findString(data_i,data_col_1)     
    if return_array:
        return array([float(line[key_index]) 
                      for line in data[data_i:data_j]])
    else:   
        return [float(line[key_index]) 
                for line in data[data_i:data_j]]



def getInputData(path=os.path.join(os.path.expanduser('~'),'ComboCode','Data'),\
                 keyword='STAR_NAME',filename='Star.dat',remove_underscore=0,\
                 make_float=1,start_index=1):
    
    """
    Search ComboCode Input Data files for parameters. 
    
    Includes files such as Dust.dat, Star.dat, Indices.dat, Molecule.dat.

    @keyword path: CC data path
    
                   (default: ~/ComboCode/Data/)
    @type path: string
    @keyword keyword: the type of information required, always equal to one of 
                      the keywords present in the "Data" of ComboCode, and 
                      automatically also a Star dict keyword
                      
                      (default: STAR_NAME)
    @type keyword: string
    @keyword filename: filename in PATH_COMBOCODE/Data/ that includes wanted 
                       information
                       
                       (default: 'Star.dat')
    @type filename: string
    @keyword remove_underscore: remove the underscores from the entries and 
                                replace them by spaces.
                                
                                (default: 0)
    @type remove_underscore: bool
    @keyword make_float: set to 0, if no floats are desired at all. If 1, all
                         entries will be converted to floats and on failure,
                         the string is returned instead
                         
                         (default: 1)
    @type make_float: bool
    @keyword start_index: Start search for keyword at this value
     
                          (default: 1)
    @type start_index: int

    @return: The requested data from the CC input
    @rtype: list
     
    """
    
    keyword = keyword.upper()
    data = [line 
            for line in readFile(os.path.join(path,filename),' ')
            if ''.join(line).strip()]
    i = int(start_index)
    while ' '.join(data[i-1]).find(keyword) == -1:
        i += 1
    data_index = [line.strip('#') 
                  for line in data[i-1]
                  if line.strip('#')].index(keyword)
    try:
        end_index = i
        while data[end_index][0][0] != '#':
            end_index += 1
    except IndexError:
        end_index = None
    try:
        if not make_float: 
            raise ValueError
        else: 
            return [float(line[data_index]) 
                    for line in data[i:end_index] 
                    if line[0]]
    except ValueError:
        if remove_underscore: 
            return [line[data_index].replace('_',' ') 
                    for line in data[i:end_index] 
                    if line[0]]
        else: 
            return [line[data_index] 
                    for line in data[i:end_index] 
                    if line[0]]



def readFile(filename,delimiter=None,replace_spaces=1):
    
    """
    Read file, and return content with delimiter of choice.
    
    The delimiter allows for lines to be split into a list of substrings.
    
    @param filename: the full filename of to be read file
    @type filename: string
    
    @keyword delimiter: The delimiter, default if strings don't have to be 
                        split into a list of substrings
                        
                        (default: None)
    @type delimiter: string
    @keyword replace_spaces: Replace any number of spaces or tabs by just one
                             space. If delimiter == ' ', then replace_spaces is
                             always active.
                             
                             (default: 1)
    @type replace_spaces: bool 
    
    @return: The lines in the file are returned, either as the full line or 
             split into substrings, depending on the delimiter
    @rtype: list[string] or list[list[string]]
    
    """
    
    if delimiter == ' ':
        replace_spaces = 1
    FILE = open(filename,'r')
    lines = FILE.readlines()
    FILE.close()
    return [line 
            for line in splitLines(lines,delimiter,replace_spaces) 
            if ' '.join(line)]



def readDict(filename,delimiter='=',comment_chars=['#'],convert_lists=0,\
             convert_floats=0,convert_ints=0,multi_keys=[]):
     
    '''
    Read a file as a dictionary.

    Commented lines and lines without the delimiter are ignored.
     
    If given keywords are present more than once, the entries for those 
    keywords are returned as a list in the dictionary with key equal to the 
    keyword. The given keywords are defined by the list multi_keys

    @param filename: the filename of the file to be read
    @type filename: string
    @keyword delimiter: the delimiter defining the key/value pairs
     
                        (default: '=')
    @type delimiter: string
    @keyword comment_chars: single character strings setting the comment characters
     
                            (default: ['#'])
    @type comment_chars: list of strings
    @keyword convert_lists: convert strings that include lists to real lists
                                     
                            (default: 0)
    @type convert_lists: bool
    @keyword convert_floats: convert input values to floats if possible
                                      
                             (default: 0)
    @type convert_floats: bool
    @keyword convert_ints: convert input values to ints if number%1 is not 0
                           only works if also convert_float == 1
                                    
                           (default: 0)
    @type convert_ints: bool
    @keyword multi_keys: Defines the keywords which may be present multiple 
                         times in the dictionary. They are included in the dict
                         as a list even if they are present only once.
                          
                         (default: [])
    @type multi_keys: list(string)
    
    @return: the dictionary with the info from the file.
    @rtype: dict
     
    '''
     
    lines = readFile(filename)
    if len(comment_chars) > 1:
        for char in comment_chars[1:]:
            lines = [line.replace(char,comment_chars[0]) for line in lines]
    lines = [line.partition(comment_chars[0])[0] for line in lines]
    #- Make sure the final character in a value definition doesn't drop off
    #- when splitting the line, in case there's no comment character on the line.
    newdict = dict()
    all_keys = [line.split(delimiter,1)[0].strip() 
                for line in lines
                if len(line.split(delimiter,1)) == 2]
    all_vals = [line.split(delimiter,1)[1].strip()
                for line in lines 
                if len(line.split(delimiter,1)) == 2]
    for mkey in multi_keys:
        mkey_count = all_keys.count(mkey) 
        if mkey_count > 0:
            newdict[mkey] = [v for k,v in zip(all_keys,all_vals) if k == mkey]
    newdict.update(dict([(k,v)
                         for k,v in zip(all_keys,all_vals)
                         if k not in multi_keys]))
    if convert_floats:
        newdict = dict([(k,convertFloat(v,convert_int=convert_ints)) 
                        for k,v in newdict.items()])
    if convert_lists:
        for k,v in newdict.items():
            if type(v) is types.StringType and v.find('[') != -1:
                v = v.strip('[').strip(']')
                newv = []
                while v.find('(') != -1:
                    tu = v[v.find('(')+1:v.find(')')].split(',')
                    tu = tuple([convertString(t) for t in tu])
                    newv.append(len(tu) != 1 and tu or tu[0])
                    v = v[v.find(')')+1:].lstrip(',')
                if not newv:
                    newv = [convertString(t) for t in v.split(',')]
                if newv == ['']:
                    newv = []
                if convert_floats:
                    converted = []
                    for newvi in newv:
                        if type(newvi) is types.TupleType:
                            converted.append(tuple([convertFloat(tui,\
                                                      convert_int=convert_ints) 
                                                    for tui in newvi]))
                        else:
                            converted.append(convertFloat(newvi,\
                                                     convert_int=convert_ints))
                    newv = converted
                newdict[k] = newv
            elif v == 'None':
                newdict[k] = None
                
    return newdict
                
              
              
def convertInt(number):
    
    '''
    Convert a float to an integer.
    
    Is only done if the float%1 != 0. Otherwise, it remains a float.
    
    @param number: The float/string to be converted.
    @type number: float or string
    
    @return: The converted integer is returned, or the input float if 
             float%1 is 0.
    @rtype: float or int
    
    '''
    
    return float(number)%1 != 0 and float(number) or int(float(number))
       


def convertString(string):
    
    '''
    Convert a string to a string, where 'None' is converted to None.
    
    @param string: The string to be converted.
    @type string: string
    
    @return: The new string (identical to input) or None
    @rtype: string
    
    '''
    
    if string == 'None':
        return None
    else:
        return string



def convertFloat(string,nans=0,convert_int=0):
    
    '''
    Convert a string to a float.
    
    If the string cannot be converted, the string itself or a nan is 
    returned.
    
    It is possible to pass objects other than strings and numbers as 'string'.
    In that case a TypeError for the float conversion is raised, and the object
    is returned untouched. The nan conversion only happens if a ValueError is 
    raised, ie it was indeed a string that does not represent a number. 
    
    @param string: The string to be converted.
    @type string: string
    
    @keyword nans: Convert the string to nan if it cannot be converted to float
    
                   (default: 0)
    @type nans: bool
    @keyword convert_int: Convert the float to integer straight away if 
                          float(string)%1 != 0
                          
                          (default: 0)
    @type convert_int: bool
    
    @return: The converted float is returned, or the input string if 
                conversion failed
    @rtype: float or string
    
    '''
    
    try:
        try:
            if convert_int:
                return convertInt(string)
            else:
                return float(string)
        except TypeError:
            return string
    except ValueError:
        return nans and float('nan') or string
        


def inputToString(val,make_int=0,exp_not=0):
    
    """
    Convert an input value to a string.
    
    @param val: The input value
    @type val: str, int, float
    
    @keyword make_int: Turn the input value into an integer
    
                       (default: 0)
    @type make_int: bool
    @keyword exp_not: Convert to exponential notation in a string
    
                      (default: 0)
    @type exp_not: bool
    
    @return: The converted input value
    @rtype: string

    """
    
    if exp_not:
        return make_int and '%.2e'%(int(float(val))) or '%.2e'%(float(val))
    else:
        return make_int and str(int(float(val))) or str(val)
    


def printRecArray(recarr,precision=8):
    
    """
    Print a record array in a mannerly fashion.
    
    @keyword precision: The precision of the floats shown when printed
    
                        (default: 8)
    @type precision: int
    
    """
    
    mlab.rec2txt(recarr,precision=precision)



def readCols(filename,delimiter=' ',make_float=1,remove_comments=1,\
             start_row=0,make_array=1,nans=0):
    
    '''
    Read columns, remove comments and turn into floats.
    
    @param filename: The full filename and path of the file
    @type filename: string
    
    @keyword delimiter: delimiter between the columns
                        
                        (default: ' ')
    @type delimiter: string
    @keyword make_float: turn everything into floats
                         
                         (default: 1)
    @type make_float: bool
    @keyword remove_comments: remove all comment lines, if False the comment 
                              lines are returned as well as a second argument.
                              In-line comments are removed regardless.
    
                              (default: 1)
    @type remove_comments: bool
    @keyword start_row: read from this row number onward
    
                        (default: 0)
    @type start_row: int
    @keyword make_array: return numpy arrays instead of python lists
    
                         (default: 1)
    @type make_array: bool
    @keyword nans: convert any non-float input value to a numpy.NaN. If False,
                   the strings are returned instead. 
                   
                   (default: 0)
    @type nans: bool
    
    @return: The columns are returned, and if requested, the commentlines
    @rtype: (list or array) or (list or array,list(string))
    
    '''
    
    lines = readFile(filename,delimiter=delimiter)[start_row:]
    if not remove_comments: 
        comment_chars = ['#','!']
        commentlines = [line 
                        for line in lines 
                        if delimiter.join(line)[0] in comment_chars]
    lines = removeComments(lines,delimiter)
    if not lines:
        print 'WARNING from DataIO.readCols! No numerical data ' + \
              'are available in %s. Returning empty list.'%filename
        return not remove_comments and ([],[]) or []
    if make_float:
        lines = [[convertFloat(l,nans=nans) for l in line] for line in lines]
    if make_array and make_float:    
        if not remove_comments: 
            return [array([line[i] for line in lines]) 
                    for i in xrange(lines[0])], commentlines
        else: 
            return [array([line[i] for line in lines]) 
                    for i in xrange(len(lines[0]))]
    else:
        if not remove_comments: 
            return [[line[i] for line in lines] 
                    for i in xrange(lines[0])], commentlines
        else: 
            return [[line[i] for line in lines] for i in xrange(len(lines[0]))]
      
      
      
def removeComments(lines,delimiter=' '):
    
    '''
    Remove all commented lines and comments in the lines.
    
    @param lines: the lines to be checked. They are assumed to be split in sub-
                  strings
    @type lines: list[list[string]]
    
    @keyword delimiter: The delimiter used when reading the file
    
                        (default: ' ')
    @type delimiter: string
    
    '''
    
    comment_chars = ['#','!']
    lines = [line 
             for line in lines 
             if delimiter.join(line)[0] not in comment_chars]
    for char in comment_chars:
        lines = [[entry.find(char) != -1 and entry[:entry.find(char)] or entry 
                  for entry in line] 
                 for line in lines]
    return lines



def splitLines(lines,delimiter=None,replace_spaces=1):

    """
    Split lines based on a delimiter of choice, which can be None.
    
    @param lines: The lines to be split up 
    @type lines: list[string]
    
    @keyword delimiter: The delimiter, default if strings don't have to be 
                        split into a list of substrings
                        
                        (default: None)
    @type delimiter: string
    @keyword replace_spaces: replace any number of spaces/tabs by just 1 space 
        
                             (default: 1)
    @type replace_spaces: bool
    @return: The lines split up in substrings if delimiter is not None. 
             Otherwise, the strings are returned as they are, or with spaces
             replaced, depending on the replace_spaces keyword. Empty lines are
             always removed
    @rtype: list[string] or list[list[string]]
    
    """

    #- Make sure a space/tab appears as one space only, if replace_spaces == 1
    if not delimiter is None and replace_spaces:
        return [" ".join(line.split()).split(delimiter) 
                for line in lines 
                if line]    
    elif delimiter is None and replace_spaces:
        return [" ".join(line.split()) for line in lines if line]
    elif delimiter is None and not replace_spaces:
        return [line for line in lines if line]
    elif not delimiter is None and not replace_spaces:
        return [line.split(delimiter) for line in lines if line]
    
    
    
def writeFile(filename,input_lines,mode='w'):
    
    """
    Write file with a list of strings as input.
    
    @param filename: filename of file, includes data type extension
    @type filename: string
    @param input_lines: The lines to be written
    @type input_lines: list[string]
    
    @keyword mode: writing mode ('w' is new file, 'a' appends to existing file)
    
                   (default: 'w')
    @type mode: string
    
    """
    
    #- Check existence folder only if new file is made
    if mode == 'w':                                
        testFolderExistence(os.path.split(filename)[0])
    FILE = open(filename,mode)
    FILE.write('\n'.join(input_lines))
    FILE.close()
    
    
    
def writeCols(filename,cols,mode='w'):
    
    """
    Write columns of data.
    
    @param filename: filename of target file for writing
    @type filename: string
    @param cols: columns to be written, every column given separately in list
    @type cols: list[list or array]
    
    @keyword mode: writing mode ('w' is new file, 'a' appends to existing file)
    @type mode: string
    
    """
    
    #- Check existence folder only if new file is made
    if mode == 'w':                                
        testFolderExistence(os.path.split(filename)[0])
    FILE = open(filename,mode)
    FILE.write('\n'.join([' '.join(['%.15e' %col[i]  
                                    for col in cols]) 
                          for i in xrange(len(cols[0]))]))
    FILE.close()


        
def findNumber(index,floats):

    """ 
    Starting from index, find the index of the next number different from zero 
    in the list. Can be at index itself!
    
    Cannot work with non-float or non-convertible-to-float input values.
    
    @param index: The starting index for the search
    @type index: int
    @param floats: The floats being searched for a non-zero number
    @type floats: list[floats]
    
    @return: The index of the next number different from zero, can be param 
             index itself.
    @rtype: int
    
    """
    
    while float(floats[index]) == 0.0:
        index += 1
    return index



def findZero(index,floats):

    """ 
    Starting from index, find the index of the next number equal to zero in 
    the list. Can be at index itself!
    
    @param index: The starting index for the search
    @type index: int
    @param floats: The floats being searched for a zero number
    @type floats: list[floats]
    
    @return: The index of the next number equal to zero, can be param index 
             itself.
    @rtype: int
    
    """
    
    while float(floats[index]) != 0.0:
        index += 1
    return index



def findFloat(index,vals):
     
    '''
    Starting from index, find the index of the next float in the list, zero or 
    non-zero. Can take strings as input! Goal is to browse until an input value
    convertible to float is found.
    
    In case that no float is found before the end of the file, the length of 
    the list is returned.
    
    If the given index is larger than the length of the list, the length is 
    also returned.
    
    @param index: The starting index in the list
    @type index: int
    @param vals: the list of strings and floats
    @type vals: list
    @return: the index of the next float in the list
    @rtype: int
    
    '''
    
    while True:
        if index >= len(vals):
            return len(vals)
        try:
            dummy = float(vals[index])
            break
        except ValueError:
            index += 1
    return index



def findString(index,vals):
     
    '''
    Starting from index, find the index of the next String in the list.
   
    In case that no float is found before the end of the file, the length of 
    the list is returned.
   
    If the given index is larger than the length of the list, the length is 
    also returned.
   
    @param index: The starting index in the list
    @type index: int
    @param vals: the list of strings and floats
    @type vals: list
   
    @return: the index of the next float in the list
    @rtype: int
   
    '''
    
    while True:
        if index >= len(vals):
            return len(vals)
        try:
            dummy = float(vals[index])
            index += 1
        except ValueError:
            break
    return index



def testFolderExistence(filepath):
    
    """
    Checking if requested folder exists and creating it if not.
    
    @param filepath: the folder name
    @type filepath: string
    
    """
    
    if filepath and not os.path.isdir(filepath):
        subprocess.call(['mkdir ' + filepath],shell=True)
        print 'Made directory: ' + filepath



def joinPdf(old,new):
    
    '''
    Join .pdf files into a single .pdf and remove the separate ones.
    
    Requires the pdftk software installed on the system.
    
    @param old: The input filenames of the .pdf files to be joined
    @type old: list[string]
    @param new: The filename of the joined .pdf file
    @type new: string
    
    '''
    
    subprocess.call([' '.join(['pdftk']+old+['cat','output',new])],shell=True)
    subprocess.call([' '.join(['rm']+old)],shell=True)



def checkLink(path,ln_path,folder=1):

    '''
    Check if a link exists between two paths, and if not create it.
    
    @param path: The name of the path the link refers to
    @type path: string
    @param ln_path: The name of the link
    @type ln_path: string
    
    @keyword folder: if you're linking a folder, instead of a file. If False,
                     a folder check will not be done. 
                     
                     (default: 1)
    @type folder: bool
    
    '''
    
    if not os.path.islink(os.path.join(ln_path,\
                                       os.path.split(path.rstrip('/'))[1])):
        if folder: testFolderExistence(path)
        subprocess.call([' '.join(['ln','-s',path,ln_path])],shell=True)
        