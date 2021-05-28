import time
import os
import sys

conffn = "conf_inv_mat_MC"
outfn = "out"
exect = [128]

def getError(correctfn, testedfn):
    i = 0
    acc = 0
    absaccum = 0
    relative = 0
    maxrel = 0
    f1 = open(correctfn, "r")
    f2 = open(testedfn, "r")

    c = f1.readline()
    while c:
        i += 1
        t = f2.readline()
        if not t:
            print("Error: output larger than target")
            return (0,0)
        acc += float(c) - float(t)
        absaccum += abs(float(c) - float(t))
        relative += abs(float(c) - float(t))/float(c)
        maxrel = max(maxrel, abs(float(c) - float(t))/float(c))
        c = f1.readline()
    
    t = f2.readline()
    if t:
        print("Error: output smaller than target")
        return (0,0)
    return (relative / float(i) , absaccum / float(i), maxrel)

if len(sys.argv) == 3:
    print getError(sys.argv[1], sys.argv[2])