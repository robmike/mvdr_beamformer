#!/usr/bin/python

import os
# import re
from sys import argv


bboxstr = '%%BoundingBox:'
endstr = '%%EndComments'

def getHeader(f):
    endstr = '%%EndComments'
    x = text.lower().find(endstr.lower())
    return text[0:x+len(endstr)]
##     l = f.readline().lower()
##     while(l[0:len(endstr)] != endstr):
##         l = f.readline().lower()
##         print l
##         if l == None:
##             return None


##     return f.read()[0:f.tell()]

## def getBbox(text):
##     bboxstr = '%%BoundingBox:'
##     m = re.search(text,'^' + bboxstr + '.*\n')
##     if m:
##         return m.group() + '\n'
##     else:
##         return None


def getEpsInDir(dirname):
    f =[os.path.join(dirname,f) for f in os.listdir(dirname) if \
        os.path.isfile(os.path.join(dirname,f)) and f.lower().endswith('.eps')]
    return f


def isBadEps(data):
    '''Check if text is badly formatted EPS'''
    head = '%!PS-AdobeFont' # This is not Eps
    # Eps should be something like %!PS-Adobe-3.0 EPSF-3.0
    if data[0:len(head)].lower() == head.lower():
        return True
    else:
        return False


def Error(Exception):
    pass


if __name__ == '__main__':

    # Get directories from arguments

    #fbbx = getEpsInDir(dirbbox)

    #t = 'error-sinr10.eps'
    #fbbx = os.path.join('graphics2',t)
    #fbad = os.path.join('graphics',t)

    fbbx = argv[1]
    fbad = argv[2]

    #fbad = os.path.join(dirbad,os.path.basename(f))
    # TODO: Do not overwrite files. Make new output files
    fbad = open(fbad,'U')
    badeps = fbad.read()
    nline = fbad.newlines #New line type for this file (one of them if many)

    
    
    if isBadEps(badeps):
        # Replace header
        fb = open(fbbx)
        h = fb.readline()
        #print h

        bbox = None

        for l in fb:
            if l[0:min(len(bboxstr),len(l))].lower() == bboxstr.lower():
                bbox = l
                break

        #print bbox



        fb.close()

        
        fbad.seek(0) # Rewind file
        h = fbad.readline() # Discard header

        badeps = fbad.read()

        # FIXME: Figure out how the hell to add the proper new line
        # for this file
        #h = h.strip('\n')
        #bbox = bbox.strip('\n')
        neweps = h + bbox + endstr + nline + '\n' + badeps
        fbad.close()

        #print '%s%s%s%s\n%s' % (h,bbox,endstr,nline,badeps)
        print neweps
        #fbad.write(badeps)
        
    

