from chompy import *

## TO DEBUG

def pctorus_solid(r=1,R=3,n=8,m=16,p=1):
    dom = cprod([intervals(2*PI)(n),
           intervals(2*PI)(m),
           intervals(1.0)(p)])
    return Map(cart2torus3d(r,R), dom, out=PolytopalComplex)


out = pctorus_solid(r=1,R=3,n=4,m=4,p=2)
out.view()
