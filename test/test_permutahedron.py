from chompy import *
from scipy.spatial.qhull import *
import numpy as np


def remove(args):
    n,seq = args
    return seq[:n] + seq[n+1:]

remove([5, range(10)])


##    DEF permutahedron (d::IsIntPos) =
##        ( project:1 ∼ rotations ∼ translation ):object
##        WHERE
##            object = (MKPOL ∼ [ID, [INTSTO ∼ LEN], K:<<1>>]):vertices,
##            vertices = permutations:(1 .. (d+1)),
##            center = Meanpoint: vertices,
##            translation = T:(1..(d+1)):(AA:-: center),
##            rotations = COMP:(((CONS ∼ AA:R):(1..d DISTR (d+1))):(PI/4))
##        END;


def permutahedron(d): 
    vertices = np.array(permutations(range(1,d+2)), dtype=double)
    cells = Delaunay(vertices).convex_hull.tolist()
    obj = SimplicialComplex(vertices.tolist(),cells)
    center = centroid(obj,range(len(vertices)))
    translation = t(range(d+1))([-coord for coord in center])
    rotations = COMP(CONS(AA(r)(DISTR([range(1,d+1),d+1])))(pi/4))
    return rotations(translation(obj)).project()


permutahedron(2).view() # OK
permutahedron(3).view() # KO:  approximation error
