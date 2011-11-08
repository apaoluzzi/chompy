from chompy import *

dom1 = intervals(2*PI)
dom2 = intervals(1)
dom = cprod([dom1(24),dom2(1)])
dom.view()

def psurf(point):
    u,v = point
    return [cos(u),sin(u),v]

Map(lambda (u,v): [cos(u),sin(u),v],dom).view()
draw(Map(lambda (u,v): [cos(u),sin(u),v],dom,PolytopalComplex))


def cylsurf(r=0.5,h=1,usteps=16,vsteps=1):
    dom1 = intervals(2*PI)
    dom2 = intervals(1)
    dom = cprod([dom1(usteps),dom2(vsteps)])
    return Map(lambda (u,v): [cos(u),sin(u),v], dom,
               PolytopalComplex).scale([r,r,h])

draw(cylsurf())
draw(cylsurf(h=3,vsteps=6))
