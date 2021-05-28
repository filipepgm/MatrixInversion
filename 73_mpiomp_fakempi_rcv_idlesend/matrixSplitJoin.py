import sys

def appendFiletoFile(main, toAppend):
    line = toAppend.readline()
    while line:
        main.write(line)
        line = toAppend.readline()

def mergeFiles(destfn, listFiles):
    f = open(destfn, "w")
    for file in listFiles:
        lilf = open(file, "r")
        appendFiletoFile(f, lilf)
        lilf.close()
    f.close()

def splitFile(fn, nParts):
    extension = "seg"+str(nParts)
    f = open(fn, "r")
    line = f.readline()
    (size, size2) = line.split(" ")
    size = int(size)
    index = line

    line = f.readline()
    for part in range(nParts):
        start = (part * size) / nParts + 1
        end = ((1+part) * size) / nParts + 1

        fnout = fn + "." + str(part)+"."+ extension
        fout = open(fnout, "w")
        fout.write(str(size) + " " + str(size) + "\n")
        
        index += str(start) + " " + str(end) + " " + fnout + "\n"

        while line:
            #print line
            (lin, col, val) = line.split(" ")
            if int(lin) >= end:
                break
            fout.write(line)
            line = f.readline()
        fout.close()

    findex = open(fn + ".index"+str(nParts), "w")
    findex.write(index)
    findex.close()

    f.close()

print sys.argv

if len(sys.argv) == 3:
    splitFile(sys.argv[1], int(sys.argv[2]))

def merge(base, prefix, maxindex):
    listFiles = [base + str(i) + prefix for i in range(maxindex)]
    outfn = base + "all" + prefix

    mergeFiles(outfn, listFiles)


#base = "1048576smallw.dat.index16."
#prefix = [".out", ".256.out", ".512.out", ".1024.out"]
#maxindex = [128, 256, 512, 1024]
#
#for i in range(len(prefix)):
#    merge(base, prefix[i], maxindex[i])
