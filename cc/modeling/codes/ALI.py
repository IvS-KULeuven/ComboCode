# -*- coding: utf-8 -*-

"""
Running ALI iteratively with the EnergyBalance.

Author: R. Lombaert

"""

import os
import subprocess      

import cc.path
from cc.tools.io import DataIO
from cc.modeling.physics import EnergyBalance as EB



def runEB_ALI(afn,ai=0,ei=0,iTmax=200,iter_conv=0,ALI_args=[],iterT_kwargs={},\
              runALIinit=0,iter_texguess=0.9,eb=None,*args,**kwargs):

    '''
    Run the energy balance module concurrently with the ALI RT code. 
    
    The method iterates between EB and ALI, updating the temperature profile
    and level populations between different iterations.
    
    ALI always starts first, so that the EB always has level populations to work
    with. Make sure the inputfile for ALI and for EB have the same initial T 
    profile.
    
    If more than one molecule is requested, ALI is ran for each molecule, then 
    EB takes into account all molecules in one go, after which the method goes
    back to ALI for one molecule at a time. 
    
    Automatically updates the ALI inputfile with filename for the temperature 
    input. Adds '_iter' to the ALI input filename
    and the level populations as well as temperature output files to help
    differentiate between the initial calculation and the iterations. Note that
    the .pop file is used by the EB, while ALI works with the .popiter file 
    (that doesn't change name at all). They should give the same populations 
    between iterations. The .popiter file is the one given in the ALI inputfile,
    while the .pop file is based on the inputfilename and is created at the end
    of an ALI run. The .pop filenames (one for every molecule) are the input for 
    EB's pop keyword and must be set by the user in either an EB inputfile or as
    keyword in kwargs.
    
    Input includes 
        - ALI input filename
        - Number of iterations to run in ALI and EB (0 if convergence is needed)
        - Less strict ALI convergence criterion for iterations 1 up to n-1
        - Additional arguments for the ALI execution
        - Additional keyword arguments for the temperature iteration (iterT)
        - Additional args/kwargs for the EnergyBalance object initialisation
    
    The EB input template is MCP by default, but a different template can be 
    passed to kwargs. Anything defined in the EB input template can be redefined
    in args and kwargs as well as by passing fn to this function call, which 
    points to the EB inputfile (not the template!).
    
    Note that ALI iterations do not share information, while the EB keeps track
    of all previous iterations. At the end of the method, the EB object is 
    returned in case you want to use it to check the different iterations after
    the calculation is done, including temperature, level populations, heating 
    and cooling terms, etc.
    
    In terms of initial populations, the user can specify whether to use 
    TexGuess = -1 or TexGuess = 0.9, i.e. start from the pops calculated in the 
    previous iteration (keep_pop on), or start from the standard initial 
    conditions.
    
    If you have an EnergyBalance object from a previous iteration or call to 
    this function, you can also pass this, and the method will continue with 
    that object instead. Make sure to adapt your iTmax and such to this object.
    
    @param afn: The inputfile name for ALI. One filename given as a string in 
                case of one molecule, multiple filenames given as strings in a 
                list in case of multiple molecules. 
    @type afn: str/list[str]
    
    @keyword ai: The amount of iterations to do in the ALI calculation. Set to 
                 the default of 0 if you wish to let ALI converge to whatever 
                 convergence criterion set in the ALI input OR by the iter_conv
                 keyword. This is not applicable to the first run of ALI.
                 
                 (default: 0)
    @type ai: int
    @keyword ei: The amount of iterations to do in the temperature calculation
                 during the energy balance. Set to the default of 0 if you wish
                 to let the energy balance converge to whatever convergence 
                 criterion is determined for/by iterT. This is not applicable to
                 the first run of EB (use imax in iterT_kwargs in that case). 
                 Cannot be 1. 
    
                 (default: 0)
    @type ei: int
    @keyword iTmax: The maximum total number of iterations allowed for the EB.
                    This is the sum of all iterations between ALI and EB, and so
                    is different from imax in EB.iterT! It is primarily used to
                    put an upper limit on the while loop of the ALI vs EB 
                    iteration.
                    
                    (default: 200)
    @type iTmax: int
    @keyword iter_conv: The convergence criterion to use in ALI during iteration
                        with the energy balance. Reverts to the requested value
                        in the ALI inputfile at the end of the iteration. If 
                        more strict than the criterion given in the ALI 
                        inputfile this keyword is ignored. Default in case the
                        same convergence criterion should be used as in the ALI
                        inputfile during the iteration.
                        
                        (default: 0)
    @type iter_conv: float 
    @keyword runALIinit: Run the initial iteration of ALI. In some cases, this 
                         model is already calculated, in which case the code 
                         starts from the existing files. 
                         
                         (default: 0)
    @type runALIinit: bool
    @keyword iter_texguess: The TexGuess value for the iterations (ie not the 
                            first calculation, in which case the value is given
                            in the inputfile). Typically this is 0.9 for 
                            standard initial conditions, but alternative could 
                            be -1 to start off from the populations calculated 
                            in the previous iteration.
                       
                            (default: 0.9)
    @type iter_texguess: float
    @keyword ALI_args: Additional arguments that are appended to the execute
                       command separated by spaces. If a single string is given,
                       only that string is added to the command. Default if no 
                       extra arguments are needed.
                       
                       (default: [])
    @type ALI_args: list[str]
    @keyword iterT_kwargs: Additional arguments for the energy balance 
                           temperature iteration. See method iterT in 
                           EnergyBalance for more information.
    
                           (default: {})
    @type iterT_kwargs: dict
    @keyword eb: An EnergyBalance object from a previous call. The method will
                 simply continue with this object rather than making a new one.
                 By default runALIinit will be off. 
                 
                 (default: None)
    @type eb: EnergyBalance()
    
    @return: The EnergyBalance object is returned with all its properties.
    @rtype: EnergyBalance()
    
    '''
    
    #-- 1) Calculate the zeroth ALI iteration for initial level populations
    #-- 2) Set up the EnergyBalance object
    #-- 3) Calculate the first EB iteration and create new input for ALI
    #-- 4) Run ALI with the new input. 
    #-- 5) Rinse and repeat until converged. 
    
    #-- Step 1: Calculate ALI given the pre-defined inputfile for all molecules
    if isinstance(afn,str): afn = [afn]
    if not eb is None: runALIinit = 0
    if runALIinit:
        print('--------------------------------------------')
        print('Running first guess for ALI.')
        print('--------------------------------------------')
        for ifn in afn: execALI(ifn,args=ALI_args)
    
    #-- Step 2: Set up the EnergyBalance object. 
    if eb is None:
        pars = {'template': 'mcp'}
        pars.update(kwargs)
        eb = EB.EnergyBalance(*args,**pars)
            
    #-- Remember the maximum requested iterations if it is present. Otherwise, 
    #   set it to the default in the EnergyBalance. 
    ei = int(ei)
    imax = iterT_kwargs['imax'] if iterT_kwargs.has_key('imax') else 50
    if eb is None:
        iterT_kwargs['imax'] = ei if int(ei) else imax
    else: 
        iterT_kwargs['imax'] = eb.i + ei if int(ei) else eb.i+imax
    
    #-- Step 3: Run the EnergyBalance
    m = 'next' if not eb is None else 'first'
    print('--------------------------------------------')
    print('Running {} guess of EnergyBalance (EB).'.format(m))
    print('--------------------------------------------')
    eb.iterT(**iterT_kwargs)
    
    #-- Prep step 4: Create new filenames for ALI input and Tkin
    fnT = afn[0].replace('.inp','_iter.temp')
    fn_new = [ifn.replace('.inp','_iter.inp') for ifn in afn]
    
    #-- Step 4: Iterate between steps 2 and 3 until the temperature iteration
    #           needs only 1 step to converge. Note that ei cannot be 1 for this
    #           to work.
    ei = 2 if int(ei) == 1 else int(ei)    
    i, iT, iter_texguess = 1, -2, float(iter_texguess)
    while eb.i != iT + 1 and eb.i <= iTmax:
        print('--------------------------------------------')
        print('Running iteration {} of ALI + EB. Current EB'.format(i) + \
              ' iteration index: {}.'.format(eb.i))
        print('--------------------------------------------')
        
        #-- Remember the current temperature iteration
        iT = eb.i
        
        #-- Step 2':
        #-- Update the ALI inputfiles and run. This overwrites the previous temp
        #   and pop files. You have access to both pops and temps of all 
        #   iterations through eb.pop[i] and eb.T_iter[i], respectively, with i  
        #   the iteration. Also sets the (maybe) new convergence criterion, new 
        #   Tkin/pop filename, and TexGuess to -1.
        DataIO.writeCols(filename=fnT,cols=[eb.r,eb.T.eval()])
        for ifn,ifn_new in zip(afn,fn_new):
            updateInputfile(fn=ifn,fnT=fnT,ai=ai,conv=iter_conv,\
                            texguess=iter_texguess,fn_new=ifn_new)
            execALI(ifn_new,args=ALI_args)
        
        #-- Step 3':
        #-- Update the level populations for all molecules. 
        for m,ifn_new in zip(eb.molecules,fn_new):
            eb.updatePop(m=m,fn=ifn_new.replace('.inp','.pop'))
        
        #-- Only if ei is not zero, limit the maximum amount of iterations.
        #   Otherwise set it to imax or the default of imax, in addition to the
        #   current iteration
        iterT_kwargs['imax'] = eb.i+ei if int(ei) else eb.i+imax
        eb.iterT(**iterT_kwargs)

        #-- Increase the EB+ALI iteration index.
        i += 1
    
    #-- If the temperature iteration reaches convergence, run the original ALI
    #   inputfile again, with updated T and pop files (ie with the original 
    #   convergence criterion and number of iterations.
    print('--------------------------------------------')
    print('Running final ALI iteration.')
    print('--------------------------------------------')
    DataIO.writeCols(filename=fnT,cols=[eb.r,eb.T.eval()])
    for ifn,ifn_new in zip(afn,fn_new):
        updateInputfile(fn=ifn,fnT=fnT,texguess=iter_texguess,fn_new=ifn_new)
        execALI(ifn_new,args=ALI_args)

    return eb
    
    
    
def updateInputfile(fn,fnT,ai=0,conv=0,texguess=0.9,fn_new=None):

    '''
    Update the ALI inputfile with new information for the next iteration. 
    
    @param fn: The full original ALI inputfilename
    @type fn: str
    @param fnT: The full filename of the temperature profile
    @type fnT: str
    
    @keyword ai: The amount of iterations to do in the ALI calculation. Set to 
                 the default of 0 if you wish to let ALI converge to whatever 
                 convergence criterion set in the ALI input OR by the iter_conv
                 keyword. 
                 
                 (default: 0)
    @type ai: int
    @keyword conv: The convergence criterion to use in ALI during iteration
                   with the energy balance. If 
                   more strict than the criterion given in the ALI 
                   inputfile this keyword is ignored. Default in case the
                   same convergence criterion should be used as in the ALI
                   inputfile during the iteration.
                   
                   (default: 0)
    @type conv: float 
    @keyword texguess: The TexGuess value for the iterations. Typically this is 
                       0.9 for standard initial conditions, but alternative 
                       could be -1 to start off from the populations calculated 
                       in the previous iteration.
                       
                       (default: 0.9)
    @type texguess: float    
    @keyword fn_new: The new (full) ALI inputfile name. Default if original is
                     to be updated.
                     
                     (default: None)
    @type fn_new: str
    
    '''
    
    ai, conv, tex_guess = int(ai), float(conv), float(texguess)
    
    #-- Change T and pop filenames
    data = changeKey(fn=fn,k='Tkin',v='I {} 1.0'.format(fnT))

    #-- Set the TexGuess to -1 so ALI reads the population filename
    tex_line = getKey(data=data,k='PassBand').split()
    tex_line[1] = str(texguess)
    data = changeKey(data=data,k='PassBand',v=' '.join(tex_line))
    
    #-- Set the number of maximum iterations if full convergence is not needed,
    #   or alternatively change the convergence criterion
    if ai or conv: 
        conv_line = getKey(data=data,k='MaxIter').split()
        if ai: conv_line[0] = str(ai)
        if conv and float(conv_line[5]) < conv: conv_line[5] = str(conv)
        data = changeKey(data=data,k='MaxIter',v=' '.join(conv_line))
    
    #-- Save the inputfile
    DataIO.writeFile(input_lines=data,filename=fn if fn_new is None else fn_new)
    
    

def execALI(fn,args=[]):

    '''
    Call the ALI executable for a given filename. 
    
    Note that the filename should be located in the ALI home folder (given in 
    usr/Path.dat), and that the ouput of the code will given there as well. 
    
    This method changes the working directory to the ALI home folder.
    
    @param fn: The full ALI inputfilename.
    @type fn: str
    
    @keyword args: Additional arguments that are appended to the execute command
                   separated by spaces. If a single string is given, only that
                   string is added to the command. Default if no extra arguments
                   are needed.
                   
                   (default: [])
    @type args: list[str]
    
    '''
    
    if isinstance(args,str): args = [args]
    os.chdir(cc.path.ali)
    fncut = os.path.split(fn.replace('.inp',''))[1]
    subprocess.call(['./ali {}'.format(fncut)+' '.join(args)],shell=True)



def getKey(k,data=None,fn=None):

    '''
    Retrieve data from an ALI inputfile. 
    
    Returns the line following the line that contains the given key. 
    
    @param k: The unique input key word for which the ALI inputfile is 
              searched.
    @type k: str

    @keyword data: The data, ie a file read by readFile with delimiter set to ''
                   A filename must be given if data is None.
                   
                   (default: None)
    @type data: list[str]
    @keyword fn: The ALI input filename. Only used if data is None.
    
                 (default: None)
    @type fn: str
    
    @return: The line following the line that contains given key
    @rtype: str
    
    ''' 
    
    if data is None:
        data = DataIO.readFile(filename=fn,delimiter=None,replace_spaces=0)
    i = DataIO.findKey(0,data,k)
    return data[i+1].replace('\n','')
    


def changeKey(k,v,data=None,fn=None):
    
    '''
    Update an ALI inputfile with new information.
    
    Takes an input key and value. The file is searched for the key, and the 
    value is inserted on the next line. 
    
    The content is then returned. This method does not save the file.

    @param k: The unique input key word for which the ALI inputfile is 
              searched.
    @type k: str
    @param v: The value inserted on the next line after key. Must be a 
              formatted string; no formatting is done here.
    @type v: str

    @keyword data: The data, ie a file read by readFile with delimiter set to ''
                   A filename must be given if data is None.
                   
                   (default: None)
    @type data: list[str]
    @keyword fn: The ALI input filename. Only used if data is None.
    
                 (default: None)
    @type fn: str
    
    @return: The data as read from a file or as adapted from the original 
    @rtype: list[str]
    
    '''
    
    if data is None:
        data = DataIO.readFile(filename=fn,delimiter=None,replace_spaces=0)
        data = [line.replace('\n','') for line in data]
    i = DataIO.findKey(0,data,k)
    data[i+1] = v
    return data
    
    

if __name__ == "__main__":
    
    tests = {1: 'ad', 2: 'drift', 3: 'dg', 4: 'dt', 5: 'a', 6: 'pe', 7: 'h2'}
    
    try:
        plot = int(sys.argv[1])
    except IndexError:
        plot = 0
    
    test = [int(argi) for argi in sys.argv[2:]]
        
    if not test: test = range(1,len(tests)+1)
    
    for t in test: 
        tests[t](int(plot))