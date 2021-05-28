import sys

def getWeights(fn):
    f = open(fn, "r")
    
    line = f.readline()
    (size, size2) = line.split(" ")
    size = int(size)
    assert size == int(size2)

    total = 0
    weights = [0.0] * size

    line = f.readline()
    while line:
        (lin, col, val) = line.split(" ")
        col = int(col)
        val = float(val)
        weights[col-1] += val
        total += val
        line = f.readline()
    f.close()

    return (size, weights, total)

def splitFile(fn, nParts):
    size, weights, total = getWeights(fn)
    target = total/nParts

    distrib = [0]
    accum = 0
    for i, w in enumerate(weights):
        accum += w
        if accum >= target:
            distrib.append(i)
            target += total/nParts
    
    distrib = distrib[0:nParts] + [size]
            

    f = open("distrib.conf", "w")
    r = 0
    for v in distrib:
        f.write( str(v) + "\n")
    f.close()


# Accept matrix filename, number of divisions
# Matrix with index starting at 1 !!!!!!!!!!!!!!!
print sys.argv

if len(sys.argv) == 3:
    splitFile(sys.argv[1], int(sys.argv[2]))

