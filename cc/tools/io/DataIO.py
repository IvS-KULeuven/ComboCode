# -*- coding: utf-8 -*-

"""
Module for reading and writing data files as well as parsing of information.

Author: R. Lombaert

"""

import os, sys
import subprocess
from glob import glob
import numpy as np
from scipy import array,zeros
from PyPDF2 import PdfFileMerger
from matplotlib import mlab

import cc.path


def read(func,module=sys.modules[__name__],return_func=0,*args,**kwargs):

    '''
    Takes any function given as a str (preceding modules separated with a dot 
    '.' ) and passes any additional args and kwargs to the function found in 
    either the DataIO module or the given module. 
    
    Note that numpy is important in DataIO as np, and can also be requested this
    way. 
    
    @param func: The requested function. Can also be given as a function, in 
                 which case no search for the function is done, and this returns
                 that function's results.
    @type func: str/function
    
    @keyword module: The home module of the function. Default is DataIO, and any
                     module therein can be accessed by passing module.func to 
                     the keyword. 
                     
                     (default: DataIO)
    @type module: module
    @keyword return_func: Return the function itself, rather than calling it. 
    
                          (default: 0)
    @type return_func: bool
    
    @return: the function's output is returned. 
    @rtype: anything. Or function. 
    
    '''

    #-- Function not yet found. Search it. 
    if isinstance(func,str):
        #-- Recursively find the function of loaded modules in given
        #   module or a function of the DataIO module itself if no '.' 
        for fstr in func.split('.'):
            if fstr == 'DataIO': 
                module = sys.modules[__name__]
                continue
            module = getattr(module,fstr)
    
    #-- Only return the function, don't call it
    if return_func:
        return module
        
    #-- Call the function and return the results.
    return module(*args,**kwargs)


def findKey(i,data,key):

    '''
    Find the index of the line that contains the first occurrence of a keyword.
    
    The search is case-insensitive.
    
    @param i: The starting index of the search
    @type i: int
    @param data: The data read with readFile. Either split on a delimiter, or 
                 the full line as a str.
    @type data: list[list]/list[str]
    @param key: The keyword that is searched
    @type key: str
    
    @return: The index of the line that contains the key for the first time 
    @rtype: int
    
    '''
    
    #-- The search is case-insensitive.
    key = key.upper()
    
    #-- Determine if lines are given as string or a list of split strings.
    if isinstance(data[0],str): 
        key_absent = lambda x: x.upper().find(key) == -1
    else:
        key_absent = lambda x: ' '.join(x).upper().find(key) == -1
    
    #-- key_absent returns True as long as it cannot find the keyword.
    #   while loop ensures we don't have to format the entire file. Can be 
    #   arduous if the file is large.
    while key_absent(data[i]):
        i += 1
    
    return i
    
    
    
def getKeyData(incr,filename,keyword,single=1):
    
    """
    Search a data file with data in a single (or multiple) columns, separated by
    comment lines containing the type of data. 
    
    The data returned follow the line that contains the key, unless incr is set 
    to 0. In this case, the line that contains the data is returned in its 
    entirety, and single is put to 0.
    
    This method is often used for extracting MCMax output. In that case, for 
    radius and theta incr is usually the grid size (NRAD and NTHETA 
    respectively), and for any other quantity incr is NRAD*NTHETA (fi for 
    denstemp.dat). incr==0 can be used to extract inputvalues from log.dat.

    @param incr: length of the data after key that is required.
                 Put this keyword to zero if you are extracting a number
                 from one line that contains the keyword itself. In that case
                 single is put to 0 so you can take your information from the
                 whole line.
    @type incr: int
    @param filename: name and path of the file searched
    @type filename: string
    @param keyword: the type of information required, always equal to one
                    of the keywords present in the file
    @type keyword: string
    
    @keyword single: return a list of only the first element on every row. 
                     Otherwise the entire line is returned. Off by default if 
                     incr == 0.
    
                     (default: 1)
    @type single: bool
    
    @return: The requested data
    @rtype: list[]
    
    """
    
    data = readFile(filename,' ')
    
    #-- The line with the key is usually not what we want. So add 1.
    i = findKey(0,data,keyword) + 1

    #-- If incr is 0, we need the line itself, and not just the first value.
    if not incr:
        single = 0
        i -= 1
        incr = 1
        
    #-- Return a single value (first of the line) or the entire line.
    if single:
        return [float(line[0]) for line in data[i:i+int(incr)]]
    else:
        return [line for line in data[i:i+int(incr)]]



def readFortranFile(convert_cols,func=np.loadtxt,*args,**kwargs):

    '''
    Reads a fortran data file. 
    
    The method is identical to np.loadtxt, but defines a converters dict for 
    converting double-notation into floats. 
    
    Use as np.loadtxt, but leave out the converters argument.
    
    @param convert_cols: The indices of the columns containing double notation
    @type convert_cols: list
    
    @keyword func: The read function used. Default is np.loadtxt, alternative
                   is np.genfromtxt. Function requires converters keyword. 
                   
                   (default: np.loadtxt)
    @type func: function
    
    @return: The np.loadtxt output
    @rtype: array
    
    '''
    
    converters = {i: lambda x:float(x.replace('D','E')) for i in convert_cols}
    return np.loadtxt(converters=converters,*args,**kwargs)    



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
        dd = array([float(line[key_index].replace('D+','E+').replace('D-','E-')) 
                    for line in data[data_i:data_j]])
        return dd
    else:   
        return [float(line[key_index].replace('D+','E+').replace('D-','E-')) 
                for line in data[data_i:data_j]]
    
    
    
def getInputData(path=cc.path.usr,keyword='STAR_NAME',filename='Star.dat',\
                 remove_underscore=0,make_float=1,start_index=1,rindex=None):
    
    """
    Search ComboCode/usr files for parameters. (Can be applied to other files
    as well)
    
    Includes files such as Dust.dat, Star.dat, Indices.dat, Molecule.dat.

    @keyword path: Location of the input file
    
                   (default: cc.path.usr)
    @type path: string
    @keyword keyword: the type of information required, always equal to one of 
                      the keywords present in the "Data" of ComboCode, and 
                      automatically also a Star dict keyword
                      
                      (default: STAR_NAME)
    @type keyword: string
    @keyword filename: filename in that includes wanted information
                       
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
    @keyword start_index: Start search for keyword on the line before this index
                          (ie data is returned from the first noncommented line 
                          at or after this index)
     
                          (default: 1)
    @type start_index: int
    @keyword rindex: Only return element with this index. Default if full list
                     is to be returned.
    
                     (default: None)
    @type rindex: int                 
   
    @return: Requested data from the usr input (either list or single element)
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
            elements = [float(line[data_index]) 
                        for line in data[i:end_index] 
                        if line[0]]
    except ValueError:
        if remove_underscore: 
            elements = [line[data_index].replace('_',' ') 
                        for line in data[i:end_index] 
                        if line[0]]
        else: 
            elements = [line[data_index] 
                        for line in data[i:end_index] 
                        if line[0]]
    if not rindex is None:
        return elements[rindex]
    else:
        return elements


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



def readDict(filename=None,lines=None,delimiter='=',comment_chars=['#'],\
             convert_lists=0,convert_floats=0,convert_ints=0,multi_keys=[],\
             start_row=0,end_row=None,key_modifier=None):
     
    '''
    Read a file as a dictionary.

    Commented lines and lines without the delimiter are ignored.
     
    If given keywords are present more than once, the entries for those 
    keywords are returned as a list in the dictionary with key equal to the 
    keyword. The given keywords are defined by the list multi_keys

    @keyword filename: the filename of the file to be read. In case of default,
                       lines are required. filename takes precedence over lines.
                    
                       (default: None)
    @type filename: string
    @keyword lines: Skip reading a file, and parse these lines instead. 
                    If default, a filename is required. Is a list of strings,
                    each interpreted as a line. Ignored if filename is given.
                    
                    (default: None)
    @type lines: list[str]
    @keyword delimiter: the delimiter defining the key/value pairs
     
                        (default: '=')
    @type delimiter: string
    @keyword comment_chars: single character strings setting the comment 
                            characters
     
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
    @keyword start_row: Limit the text file to lines starting from this index.
    
                          (default: 0)
    @type start_row: int
    @keyword end_row: Limit the text file to lines ending before this index
                        (ie last line is end_index - 1). Default includes up to
                        the last line of the file.
                        
                        (default: None)
    @type end_row: int
    @keyword key_modifier: A function that modifies the key, such as lower(). 
                           Any method that works on a string in principle works.
                           Give the function as a string.
                           
                           (default: None)
    @type key_modifier: str
    
    @return: the dictionary with the info from the file.
    @rtype: dict
     
    '''
    
    if filename:
        lines = readFile(filename)
    lines = lines[start_row:end_row]
    lines, comments = removeComments(lines,comment_chars=comment_chars)
    
    #-- Make sure the final character in a value definition doesn't drop off
    #   when splitting the line in case there's no comment character on the line
    newdict = dict()
    all_keys = [line.split(delimiter,1)[0].strip()
                for line in lines
                if len(line.split(delimiter,1)) == 2]
    all_vals = [line.split(delimiter,1)[1].strip()
                for line in lines 
                if len(line.split(delimiter,1)) == 2]
    
    #-- Apply a key modifier if requested. 
    if key_modifier: 
        all_keys = [getattr(k,key_modifier)() for k in all_keys]
    
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
            if isinstance(v,str) and v.find('[') != -1:
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
                        if isinstance(newvi,tuple):
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
    
    @param recarr: The record array.
    @type recarr: recarray
    
    @keyword precision: The precision of the floats shown when printed
    
                        (default: 8)
    @type precision: int
    
    """
    
    print mlab.rec2txt(recarr,precision=precision)



def removeComments(lines,comment_chars=['#','!',';']):

    '''
    Split input from comments and return both as separate lists.
    
    Takes a list of strings, such as what readFile returns. 
    
    @param lines: The strings in a file
    @type lines: list[str]
    
    @keyword comment_chars: single character strings setting the comment 
                            characters
     
                            (default: ['#','!'])
    @type comment_chars: list[str]
    
    @return: The input and the comment lines in two separate lists
    @rtype: (list,list)
    
    '''

    if len(comment_chars) > 1:
        for char in comment_chars[1:]:
            lines = [line.replace(char,comment_chars[0]) for line in lines]
    data = [line.partition(comment_chars[0])[0] for line in lines]
    comments = [line.partition(comment_chars[0])[2] for line in lines]
    return (data,comments)
    


def readCols(filename,delimiter=' ',make_float=1,start_row=0,make_array=1,\
             nans=0,start_from_keyword='',return_comments=0,\
             comment_chars=['#','!',';'],end_row=None):
    
    '''
    Read columns, remove comments and turn into floats.
    
    Note that number of columns returned is the minimum number of columns from
    all rows, where the columns are split by the chosen delimiter (space-like
    by default). 
    
    @param filename: The full filename and path of the file
    @type filename: string
    
    @keyword delimiter: delimiter between the columns
                        
                        (default: ' ')
    @type delimiter: string
    @keyword make_float: turn everything into floats
                         
                         (default: 1)
    @type make_float: bool
    @keyword start_row: Limit the text file to lines starting from this row. 
                        If start_from_keyword is used, start_row counts from the
                        index where the keyword is first found.
    
                        (default: 0)
    @type start_row: int
    @keyword make_array: return numpy arrays instead of python lists
    
                         (default: 1)
    @type make_array: bool
    @keyword nans: convert any non-float input value to a numpy.NaN. If False,
                   the strings are returned instead. 
                   
                   (default: 0)
    @type nans: bool
    @keyword start_from_keyword: Start returning data from the line that
                                 contains a given keyword. Only used if not 
                                 default. The start_row parameter counts from 
                                 this index onward. (not case sensitive)
                                 
                                 (default: '')
    @type start_from_keyword: string
    @keyword return_comments: Return the comments list in addition to the data
    
                              (default: 0)
    @type return_comments: bool
    @keyword comment_chars: single character strings setting the comment 
                            characters
     
                            (default: ['#','!',';'])
    @type comment_chars: list[str]
    @keyword end_row: Limit the text file to lines ending before this index
                      (ie last line is end_index - 1). Default includes up to
                      the last line of the file. If start_from_keyword is given,
                      end_row counts with respect to the found index.
                        
                      (default: None)
    @type end_row: int
    
    @return: The columns are returned, with in addition the comments if 
             requested
    @rtype: list[list or array] or (list[list or array],list[str])
    
    '''
    
    lines = readFile(filename)
    if str(start_from_keyword):
        #-- Find occurrences of searchstring
        start_from_keyword = start_from_keyword.upper()
        indices = [i for i,line in enumerate(lines) 
                     if line.upper().find(start_from_keyword) != -1]
        #-- If any were found, grab the first and cut the lines above it
        if indices: lines = lines[indices[0]:]
    
    #-- Cut anything above start row, up to the end_row
    lines = lines[start_row:end_row]
    
    #-- Remove the comments and empty lines, then split the lines
    lines,comments = removeComments(lines,comment_chars=comment_chars)
    lines = [line for line in lines if line]
    lines = splitLines(lines,delimiter=delimiter)
    if not lines:
        print 'WARNING from DataIO.readCols! No numerical data ' + \
              'are available in %s. Returning empty list.'%filename
        return []
    
    #-- Apply requests for floatsand arrays and return.
    if make_float:
        lines = [[convertFloat(l,nans=nans) for l in line] for line in lines]
    ndata = min([len(line) for line in lines])
    if make_array and make_float:    
        lines = [array([line[i] for line in lines]) for i in xrange(ndata)]
    else:
        lines = [[line[i] for line in lines] for i in xrange(ndata)]
    return return_comments and (lines,comments) or lines



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
    
    
    
def writeFile(filename,input_lines,mode='w',delimiter='\n'):
    
    """
    Write file with a list of strings as input.
    
    @param filename: filename of file, includes data type extension
    @type filename: string
    @param input_lines: The lines to be written
    @type input_lines: list[string]
    
    @keyword delimiter: The line separator. Typically the new line character, 
                        but can be changed if needed (eg empty string in case
                        new line character is already included in the input 
                        lines)
                        
                        (default: '\\n')
    @type delimiter: string
    @keyword mode: writing mode ('w' is new file, 'a' appends to existing file)
    
                   (default: 'w')
    @type mode: string
    
    """
    
    #- Check existence folder only if new file is made
    if mode == 'w':                                
        testFolderExistence(os.path.split(filename)[0])
    FILE = open(filename,mode)
    FILE.write(delimiter.join(input_lines))
    FILE.close()
    


def replaceString(filename,old_str,new_str):
    
    '''
    Replace a given string with another in requested filename.
    
    The filename can contain whatever works for the glob function, such as the 
    wildcard character, to do multiple files in one go. 

    @param filename: The filename, possibly with a wildcard character
    @type filename: string
    @param old_str: The old string to be replace in all the files
    @type old_str: str
    @param new_str: The new string to be inserted in all the files
    @type new_str: str
    
    '''
    
    gg = glob(filename)
    for gf in gg:
        old_lines = readFile(gf,replace_spaces=0)
        new_lines = [ll.replace(old_str,new_str) for ll in old_lines]
        writeFile(filename=gf,input_lines=new_lines,delimiter='') 
    
    
    
def writeCols(filename,cols,mode='w',delimiter='\t'):
    
    """
    Write columns of data.
    
    @param filename: filename of target file for writing
    @type filename: string
    @param cols: columns to be written, every column given separately in list
    @type cols: list[list or array]
    
    @keyword mode: writing mode ('w' is new file, 'a' appends to existing file)
    
                   (default: 'w')
    @type mode: string
    @keyword delimiter: The delimiter used between columns
    
                        (default: '\t')
    @type delimiter: string
    
    """
    
    #- Check existence folder only if new file is made
    if mode == 'w':                                
        testFolderExistence(os.path.split(filename)[0])
    FILE = open(filename,mode)
    FILE.write('\n'.join([delimiter.join(\
                            [isinstance(col[i],str) and '%s'%col[i] \
                                or '%.3e'%col[i]  
                             for col in cols])
                          for i in xrange(len(cols[0]))]))
    FILE.close()


        
def findNumber(index,floats):

    """ 
    Starting from index, find the index of the next number different from zero 
    in the list. Can be at index itself!
    
    Cannot work with non-float or non-convertible-to-float input values.
    
    If at the end of the list, the index is returned as the length of the list.
    
    @param index: The starting index for the search
    @type index: int
    @param floats: The floats being searched for a non-zero number
    @type floats: list[floats]
    
    @return: The index of the next number different from zero, can be param 
             index itself.
    @rtype: int
    
    """
    
    while index < len(floats) and float(floats[index]) == 0.0:
        index += 1
    return index



def findZero(index,floats):

    """ 
    Starting from index, find the index of the next number equal to zero in 
    the list. Can be at index itself!

    If at the end of the list, the index is returned as the length of the list.

    @param index: The starting index for the search
    @type index: int
    @param floats: The floats being searched for a zero number
    @type floats: list[floats]
    
    @return: The index of the next number equal to zero, can be param index 
             itself.
    @rtype: int
    
    """
    
    while index < len(floats) and float(floats[index]) != 0.0:
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
            dummy = float(vals[index].replace('D+','E+').replace('D-','E-'))
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
            dummy = float(vals[index].replace('D+','E+').replace('D-','E-'))
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



def joinPdf(old,new,del_old=1):
    
    '''
    Join .pdf files into a single .pdf and remove the separate ones.
    
    Requires the pdftk software installed on the system.
    
    @param old: The input filenames of the .pdf files to be joined
    @type old: list[string]
    @param new: The filename of the joined .pdf file
    @type new: string
    
    @keyword del_old: delete the old filenames
    
                      (default: 1)
    @type del_old: bool
    
    '''
    
    pp = PdfFileMerger()
    for ofn in old:
        pp.append(ofn)
    pp.write(new)
    pp.close()
    if bool(del_old):
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
        


def fillOutSpaces(string,nchars):
    
    '''
    Fill out a string to a set number of characters with spaces. 
    
    @param string: The string to be filled out, if nchars is < than len(string)
                   the string itself is returned, unchanged.
    @type string: str
    @param nchars: The number of characters requested in final string
    @type nchars: int
    
    @return: The filled out string.
    @rtype: str
    
    '''
    
    string, nchars = str(string), int(nchars)
    if nchars > len(string):
        string += ''.join([' ']*(nchars-len(string)))
    return string



def checkEntryInfo(input_list,number_of_keys,info_type):
    
    '''
    Specific input keywords for ComboCode (currently: MOLECULE, TRANSITION,
    R_POINTS_MASS_LOSS) require multiple arguments separated by a space. This 
    method sorts those arguments and checks if the format is OK.
    
    The method returns a tuple or a list of these arguments depending on 
    multiplicity of input, and the need to create a grid for the parameter.
    
    @param input_list: argument strings for one of the special CC input keys.
    @type input_list: list[string]
    @param number_of_keys: the number of arguments expected
    @type number_of_keys: int
    @param info_type: MOLECULE, TRANSITION or R_POINTS_MASS_LOSS
    @type info_type: string
    
    @return: the sorted arguments for the requested info_type. The output is 
             given as a list or a tuple, depending on if this gridded parameter
             or not, respectively.
    @rtype: list[tuple(list[string])] or tuple(list[string])
    
    '''
    
    if not input_list:
        return []
    
    if not info_type in ('MOLECULE','R_POINTS_MASS_LOSS','TRANSITION'):
        raise KeyError('The info_type keyword is unknown. Should be ' + \
                       'MOLECULE, TRANSITION or R_POINTS_MASS_LOSS.')
    input_list = [line.split() for line in input_list]

    #-- Make sure identical transitions don't get confused if offset 0 is 
    #   given in different ways, such as more significant numbers. Offset is
    #   always the entry with index 10. 
    def checkOffset(line):
        if float(line[10]) == 0.0:
            line[10] = '0.0'
        return line
    if info_type == 'TRANSITION':
        input_list = [checkOffset(line) for line in input_list]

    if set([len(line) for line in input_list]) != set([number_of_keys]):
        print 'Number of keys should be: %i'%number_of_keys
        print '\n'.join(['%i  for  %s'%(len(line),line) for line in input_list])
        raise IOError('Input for one of the %s lines has wrong '%info_type + \
                      'number of values. Double check, and abort.')
    else:
        #-- if MOLECULE: only molecule string, 
        #   if R_POINTS_MASS_LOSS: only grid id number, 
        #   if TRANSITION: everything except last entry (n_quad)
        entries = [info_type in ('MOLECULE','R_POINTS_MASS_LOSS') \
                        and line[0] or ' '.join(line[0:-1]) 
                   for line in input_list]  
        unique_entries = set(entries)
        #-- ie the defining parameter is never multiply defined ! 
        #   Hence, only a single set of parameters is used here.
        if len(unique_entries) == len(entries):        
            if info_type == 'R_POINTS_MASS_LOSS':
                input_list = ['  '.join(il[1:]) for il in input_list]
            return tuple(input_list)
        else:
            if info_type == 'TRANSITION':
                indices = [i 
                           for i,entry in enumerate(entries) 
                           if entries.count(entry) > 1]
                print 'Identical transition(s) in the transition list: Doubl'+\
                      'es will not be removed even if N_QUAD is the same too!'
                for i in indices:
                    print 'At index %i:  %s' %(i,entries[i])
                raw_input('Abort if identical transitions are not expected. '+\
                          'Press enter otherwise.')
            if info_type == 'R_POINTS_MASS_LOSS':
                #-- This will be a list of R_POINTS_MASS_LOSS sets, where each 
                #   set is defined as a list of radial grid point parameters
                final = []  
                while input_list:
                    if int(input_list[0][0]) != 1:
                        raise IOError('The grid point ID numbers for the ' + \
                                      'R_POINTS_MASS_LOSS keywords do not ' + \
                                      'follow a correct order. Use ID = 1 ' + \
                                      'as a reset for new set of grid points.')
                    final.append([input_list.pop(0)])
                    while input_list and int(input_list[0][0]) != 1:
                        final[-1].append(input_list.pop(0))
                #-- Remove the grid point ID numbers, not relevant anymore
                #   put the rest back into a string
                final = [tuple(['  '.join([num 
                                           for i,num in enumerate(this_point) 
                                           if i != 0]) 
                                for this_point in this_set]) 
                         for this_set in final]      
                return final
            final = [[line 
                      for line,entry in zip(input_list,entries) 
                      if entries.count(entry) == 1]]
            multiples = set([entry 
                             for entry in entries 
                             if entries.count(entry) != 1])
            for entry in multiples:
                this_input = [line 
                              for this_entry,line in zip(entries,input_list) 
                              if this_entry == entry]
                final = [this_list + [extra_line] 
                         for extra_line in this_input 
                         for this_list in final]
            return [tuple(final_list) for final_list in final]


def getChemistryAbundances(filename):
    
    '''
    Reads in the Chemistry abundance output, works for both 
    fractional abundances and number densities.
    
    @param filename: The filename of the abundance output 
                     (csfrac.out or csnum.out)
    @type filename: string
    
    @return: Recursive array containing the abundance per species name.
    @rtype: recarray
    
    '''
    
    ###- Open file, read in all columns
    ##data = DataIO.readCols(filename,start_row=1)
    ###- Join them into one output array
    ##C = []
    ##for c in data[1:]:
        ##C =  np.concatenate((C,c),axis = 0)

    #- Open file, read in all lines
    f = open(filename, 'r')
    lines = f.read().splitlines()
    f.close()
    
    #- Remove first line, empty strings within lists, and empty lists
    data = [filter(None, line.split(' ')) for line in lines[1:]]
    data = filter(None, data)
    
    #- Number of columns is constant throughout the file
    if len(set([len(d) for d in data])) == 1:
        data = DataIO.readCols(filename,start_row=1)
        
        #- Join them into one output array
        C = []
        for c in data[1:]:
            C =  np.concatenate((C,c),axis = 0)
            
        #- Number of calculations per species
        c0 = data[0]
        L = np.where(np.array(c0) == c0[0])[0][1]+1        

    
    #- If the number of columns varies througout the file (e.g. by adding 
    #  species), run a more elaborate method to read in the columns    
    else:
        #- Read in first block of 10 columns and concatenate
        limit = [i for i,d in enumerate(data) if len(d) != 10][0]
        blok = zip(*data[:limit])
        C = []
        for c in blok[1:]:
            C =  np.concatenate((C,c),axis = 0)
        
        # Add the appendix (with less then 10 columns)
        app = zip(*data[limit:])
        for c in app[1:]:
            C =  np.concatenate((C,c),axis = 0)
        
        #- Number of calculations per species
        c0 = blok[0]
        L = np.where(np.array(c0) == c0[0])[0][1]+1        
    
    #- Put the names of the species in an array
    names = []      
    for ii in range(len(C)):
        if ii%L == 0:
            names.append(C[ii])  
    N = len(names)
    
    #- Radii of calculation
    radius = np.array(c0[1:L-1])
    radius = radius.astype(np.float)
    
    #- Output array: [species,[output]]
    species = np.recarray(shape = [L-2,], dtype = zip(names, [float]*N))
    for ii in range(N):
        species[names[ii]] = C[(ii*L)+1:((ii+1)*L)-1].astype(float)
    
    return species 


def getChemistryPhysPar(filename, keyword):
    
    '''
    Reads in the Chemistry physical output
    
    @param filename: The filename of the abundance output 
                     (csphyspar.out)
    @type filename: string
    
    @keyword keyword: The physical parameter in question. Options are:
                      RADIUS, n(H2), TEMP. A_V, RAD. FIELD, 
                      CO K(PHOT), VELOCITY
                      
    @type keyword: string
    
    @return: Recursive array containing the abundance per species name.
    @rtype: recarray
    
    '''
    
    #- Open file, read in all columns
    data = DataIO.readCols(filename,start_row=1)
    
    #- Initialise keyword
    keyword = keyword.upper()
    
    #- Select right columm
    c = [i for i,s in enumerate(data) if keyword in s[0]][0]
    par = [float(p) for p in data[c][1:]]
    
    return par


def getChemistrySpecies(filename,parents=1):
    
    '''
    Reads the species and parent species included in the Chemistry code.
    
    @param filename: The .specs file
    @type filename: string
    @keyword parents: Give parent species as output.
                      If parents=0, all species are given.
                      (default: 1)
    
    @return: (Parent) species included in the Chemistry code.
    @rtype: list
    
    '''
    
    #- Read in file
    data = DataIO.readFile(filename)
    data = [x.split() for x in data]
    
    #- Determine the different sections in the species file
    separators = np.where([len(i)==2 for i in data])[0]
    
    if parents:
        #- Select parent species
        parents = data[separators[2]:]
        parents = [parents[x][0] for x in range(len(parents))]
    else:
        #- Select species included
        species = data[1:separators[0]]
        species = [species[x][1] for x in range(len(species))]
        
    return species
  
  

