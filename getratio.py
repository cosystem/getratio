#!/usr/bin/env python
#
###############################################
#
# python script to calculate the #Observation/#sites ratio
# at specified resolution
#
# ussage: ./getratio.py [hklin] [pdbin] [reso_min] [reso_max] [reso_step]
# output: ratiolist.csv
#
# Author: Guanya PENG
# Aug 13, 2014
#
##############################################

import re
import sys
import csv
import os.path
import cmath, math
import argparse
import numpy as np


class fileopt(object):
    '''class of file options including readlines method so far'''
    def __init__(self, filename):
        self.filename = filename
    def read_all_lines(self, strip=None, lst=None):
        with open(self.filename, 'r') as f_in:
            if strip == None:
                lines = (line for line in f_in)
            elif strip == 'strip':
                lines = (line.strip() for line in f_in)
            else:
                print 'Bad strip option!'
            lines = list(line for line in lines if line)
            if lst == None:
                return lines
            elif lst == 'lst':
                lines = [line.split() for line in lines]
                return lines
            else:
                print 'Bad lst option!'
                exit(1)

class pdbopt(fileopt):
    '''class of pdb options inheritanted from fileopt() reading unit-cell and counting number of user specified sites'''
    def __init__(self, filename):
        self.filename = filename
        self.lines_pdb_list = fileopt(self.filename).read_all_lines()

    def getunitcell(self):
        '''read unit-cell'''
        crystlines = []
        for line in self.lines_pdb_list:
            if re.search('(^CRYST1)', line):
                crystlines.append(line)
            else:
                pass
        if len(crystlines) > 1:
            for line in crystlines:
                if len(line) == 70:
                    uniquecrystline[:0] = [line]
                else:
                    pass
        else:
             uniquecrystline = list(crystlines)
        try:
            return [float(x) for x in uniquecrystline[0][6:54].split()]
        except ValueError:
            print "Unit-Cell reading is wrong!"

    def getAtomNum(self, elementname=None):
        '''count user specified sites. Input should be the element name string'''
        allatoms = []
        for line in self.lines_pdb_list:
            if re.search('(^ATOM|^HETATM)', line):
                allatoms.append(line[76:78].strip())
            else:
                pass
        atomnum = allatoms.count(elementname)
        return atomnum

def abslist(list_in):
    '''generate an absolute list'''
    for index, item in enumerate(list_in):
        if item < 0:
            list_in[index] *= -1
        else:
            pass
    return list_in

def deleteDupList(longlist):
    '''detele the duplicate sub-lists in a list'''
    tuple_set = set(tuple(sublist) for sublist in longlist)
    listunique = [list(subtuple) for subtuple in tuple_set]
    return listunique

def crossList3(x1_list, x2_list):
    '''calculate inner product of vectors'''
    return [
            x1_list[1]*x2_list[2],
            x1_list[0]*x2_list[2],
            x1_list[0]*x2_list[1]
            ]

def complist(x1_list, x2_list):
    '''return element in x1_list which not exist in x2_list'''
    s2 = set(x2_list)
    comp = [x for x in x1_list if x not in s2]
    return comp

def chop(arr_in, tol=1e-10):
    '''chop all values in an array'''
    try:
        arr_in.real[abs(arr_in.real) < tol] = 0.0
        arr_in.imag[abs(arr_in.imag) < tol] = 0.0
    except ValueError:
        arr_in.real[abs(arr_in.real) < tol] = 0.0
    return arr_in

def Tensor(UnitCell_List):
    '''calculate the transformation tensor for h k l -> Sin(theta)^2/lambda^2'''
    abc = UnitCell_List[0:3]
    t = crossList3(abc,abc)
    angles = UnitCell_List[-3:]
    cosAngles=[np.cos(np.radians(x)) for x in angles]
    sinAngles=[np.sin(np.radians(x)) for x in angles]
    # components of cross vectors
    c = crossList3(cosAngles, cosAngles)
    s = crossList3(sinAngles, sinAngles)
    # volumn of unit-cell
    Vcell =(
        reduce(
                lambda x,y: x*y,
               abc
        )
    *
    cmath.sqrt(
        1-reduce(
                lambda x,y: x+y,
                [angle**2 for angle in cosAngles]
                )
        +
        2*reduce(
                lambda x,y: x*y,
                cosAngles
                )
        )
    )
    # conjugate variables in the reciprocal space
    abcStar = [x/Vcell for x in [a*b for a,b in zip(t,sinAngles)]]
    cosAnglesStar=[(a1-a2)/a3 for a1,a2,a3 in zip(c,cosAngles,s)]
    resul = np.diag([x**2 for x in abcStar])
    # calculate the tensor
    for i in xrange(3):
        for j in xrange(i+1,3):
            resul[i,j] = abcStar[i]*abcStar[j]*cosAnglesStar[complist(xrange(3),[i,j])[0]]
            resul[j,i] = resul[i, j]
    return chop(resul, 1e-15)

def Sthol2(tensor_Array, hkl_List):
    '''get sin(theta)^2/lambda^2 from h, k, l indexes'''
    hklarray = [np.array(hkl) for hkl in hkl_List]
    return [np.dot(hkl, np.dot(tensor_Array, hkl)).real/4. for hkl in hklarray]

def newrange(minval, maxval, step):
    a = np.arange(minval, maxval, step)
    if chop(np.array([abs(maxval-a[-1])]), 1e-10) != 0:
        return np.append(a,maxval)
    else:
        return a

def getratiolist(sthol2_List, sitesnum, minreso, maxreso, step):
    '''get the #Observation/#sites list with specified element and resolution range'''
    resolist = newrange(minreso, maxreso, step)
    ratiolist=[]
    for i in resolist:
        try:
            ratiolist.append(len([x for x in sthol2list if np.sqrt(1/(4*x)) >= i])/sitesnum)
        except ZeroDivisionError:
            sys.exit("There is no defined element in the pdb file. Exit.")
    return ratiolist

def is_numstring(s):
    '''test number string'''
    try:
        float(s)
        return True
    except ValueError:
        return False

# read the arguments
try:
    for arg in sys.argv[3:]:
        if not is_numstring(arg):
            sys.exit("Resolution arguments must be numbers. Exit.")
    hklfile = sys.argv[1]
    pdbfile = sys.argv[2]
    resomin = float(sys.argv[3])
    resomax = float(sys.argv[4])
    resostep = float(sys.argv[5])
except Exception:
    print "\nussage: ./getratio.py [hklin] [pdbin] [reso_min] [reso_max] [reso_step]\n"
    sys.exit(1)

# check the existence of input files
if not os.path.isfile(hklfile) and not os.path.isfile(pdbfile):
    print 'No such files: \'%s\' and \'%s\'. Exit' %(hklfile, pdbfile)
    exit(1)
elif not os.path.isfile(hklfile) and os.path.isfile(pdbfile):
    print 'No such file: \'%s\'. Exit.' %hklfile
    exit(1)
elif not os.path.isfile(pdbfile) and os.path.isfile(hklfile):
    print 'No such file: \'%s\'. Exit.' %pdbfile
    exit(1)
else:
    pass

# open hkl file and readlines
hkllines = fileopt(hklfile).read_all_lines()
while True:
    try:
        # get h, k, l list from the lines of hkl file
        hkllist = [[int(x) for x in [line[0:4].strip(), line[4:8].strip(), line[8:12].strip()]] for line in hkllines]
    except ValueError:
        print 'Something wrong with the hkl file: \'%s\'. Please check its content. Exit' %hklfile
        exit(1)
    else:
        # convert all the Friedles to |h|,|k|,|l|
        abshkllist = [abslist(sub) for sub in hkllist]
        # generate unique reflection list by getting rid of duplicate |h|, |k|, |l|
        uniquehkllist = deleteDupList(abshkllist)
        break

while True:
    try:
        # read pdb file and extract unit-cell parameters
        unitcell = pdbopt(pdbfile).getunitcell()
    except IndexError:
        print 'Something wrong with the pdb file: \'%s\'. Please check its content. Exit' %pdbfile
        exit(1)
    else:
        # count the number of sulfur(S)
        numsites = pdbopt(pdbfile).getAtomNum('P')
        # calculate transformation tensor
        tensor = Tensor(unitcell)
        # convert h k l to sin(theta)^2/lambda^2
        sthol2list = Sthol2(tensor, uniquehkllist)
        break

# get the ratio list
while True:
    try:
        resolist = newrange(resomin, resomax, resostep)
    except IndexError:
        print "Unreasonable resolution range: min %s, max %s, step %s\nPlease re-define your resolution range and step. Exit."\
                %(sys.argv[3], sys.argv[4], sys.argv[5])
        sys.exit(1)
    else:
        np.seterr(divide='ignore') # ignore divided by 0 error in np.sqrt(1/(4*x) in the getratiolist function
        ratiolist = getratiolist(sthol2list, numsites, resomin, resomax, resostep)
        break

# format the result and write to output file
outtuples = [(round(tup[0],2), int(round(tup[1]))) for tup in zip(resolist, ratiolist)]
with open('ratiolist.csv', 'wb') as csvfile:
    spamwriter = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
    spamwriter.writerow(['Resolution (A)', 'ratio'])
    for tup in outtuples:
        spamwriter.writerow(tup)

__author__ = "Guanya Peng"
__copyright__ = "Copyright 2014, P-SAD Project"
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Guanya Peng"
__email__ = "guanya.peng@psi.ch"
__status__ = "Production"
