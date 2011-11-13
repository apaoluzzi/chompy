from chompy import *

dom1 = intervals(2*pi)
dom2 = intervals(1)
dom = cprod([dom1(12),dom2(1),dom2(1)])
dom.view()

def psurf(point,r=0.5):
    u,v,w = point
    return [(1.0-r*w)*cos(u), (1.0-r*w)*sin(u), v]

Map(psurf, dom, PolytopalComplex).view()

a = Map(psurf, dom, PolytopalComplex)
draw(a.boundary())


#draw(Map(psurf, dom, "PolytopalComplex"))
##
##
##def cylsurf(r=0.5,h=1,usteps=16,vsteps=1):
##    dom1 = intervals(2*PI)
##    dom2 = intervals(1)
##    dom = cprod([dom1(usteps),dom2(vsteps)])
##    return Map(lambda (u,v): [cos(u),sin(u),v], dom,
##               PolytopalComplex).scale([r,r,h])
##
##draw(cylsurf())
##draw(cylsurf(h=3,vsteps=6))
