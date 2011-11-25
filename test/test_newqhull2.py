from chompy import *

if __name__=="__main__":
    obj = simplexGrid([[1],[1],[1]])
    myprint("decWingedRepr(obj)",decWingedRepr(obj))
    myprint("boundaryWingedRepr(obj)",boundaryWingedRepr(obj))



##def issimplex(obj):
##    """
##    Predicate to check if 'obj' is a simplex (i.e. a singleton complex).
##    Return a truth value.
##    """
##    return (type(obj) == SimplicialComplex) \
##            and (len(obj.cells[-1]) == 1) \
##            and (len(obj.vertices.points) == (obj.dim+1))
##
##
##def join(args):
##    """
##    To compute the convex join of two non-intersecting convex arguments.
##    Return a (decompositive) winged representation.
##    """
##    conv0,conv1 = args
##
##    if conv0.rn == conv1.rn: # same embedding space (same number of coordinates)
##
##        if AND(AA(issimplex)([conv0,conv1])):# conv0,conv1 = simplex,simplex
##            if ((conv0.dim+1) + (conv1.dim+1) <= (conv0.rn+1)):
##                points = AA(eval)(conv0.vertices.dict.keys()+conv1.vertices.dict.keys())
##                cells = [range(len(points))] # TODO: test for dimension independence
##                return SimplicialComplex(points,cells)
##
##        elif OR(AA(issimplex)([conv0,conv1])):
##            if issimplex(conv1):            # conv0,conv1 = simplex,convex
##                conv0,conv1 = conv1,conv0
##
##            """ build  winged representation of conv1"""
##            wingedStorage = [[conv1.cells[-1][k], boundary(conv1)(k)]
##                             for k in range(len(conv1.cells[-1]))]
##            myprint("wingedStorage",wingedStorage)
##                
##            
##
##            #return SimplicialComplex(points,cells)
##
##        else:                               # conv0,conv1 = convex,convex
##            pass
##
##    else: print "error: join no possible; different spaces."
##    return conv1
##
##
##s0 = SimplicialComplex([[0,0,0]],[[0]])
##s1 = SimplicialComplex([[1,0,0],[0,1,0]], [[0,1]])
##s2 = SimplicialComplex([[1,0,0],[0,1,0],[0,0,1]], [[0,1,2]])
##
##draw(s0)
##join([s0,simplexGrid([[1],[1],[1]]).translate([1,1,1])])
##                
##
##assert(issimplex(s0)==True)
##assert(issimplex(s1)==True)
##assert(issimplex(s2)==True)
##assert(issimplex(INSR(join)([s0, translate([0,0,1])(s0), s1]))==True)
##
####draw(s0)
####draw(s1)
####draw(s2)
####draw(join([s0,s1]))
####draw(join([s0,s2]))
####draw(join([s1,rotate(1,2,PI/4)(translate([0.,0.,1.])(s1))]))
####draw(join([s1,translate([0.,0.,1.])(rotate(1,2,PI/4)(s1))]))
####draw(join([s0,translate([0,0,1])(s0)]))
####draw(INSR(join)([s0,translate([0,0,1])(s0),s1]))
