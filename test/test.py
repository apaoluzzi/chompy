from chompy import *

def affineCoords(basis):
    """
    To compute the affine coordinates of a point w.r.t. a reference simplex.

    Return the coded string of affine coordinates.
    """
    myprint("basis",basis)
    d = len(basis)
    mapping = mat([vect for vect in basis]).I
    def affineCoords0(point):
        affine = (mat(point+[1.])*mapping).tolist()[0]
        return code(affine)
    return affineCoords0


def simplex(verts,cell):
    """
    To generate a simplex square matrix in homogeneous coordinates.

    Return a list of lists of numbers, taken from verts coordinate list.
    """
    return [array(verts[k] + [1.]) for k in cell]



def greatSimplex(points):
    """
    To extract a reference simplex from a list of 'points'.

    Return a set of d+1 affinely independent points.
    """
    def randomSimplex(verts,d):
        
        def closetozero(number):
            if abs(number)<1.0E-5:
                return True
            else:
                return False
            
        cell,m = [], len(verts)
        while len(cell) < d+1:
            index = int(random.random()*(m))
            if index not in cell:
                cell += [index]
            if len(cell) == d+1 and \
                closetozero(linalg.det(mat(simplex(verts,cell)))):
                cell = []
        return [verts[k] for k in cell]

    myprint("points in greatSimplex",points)
    d = len(points[0])
    reference = randomSimplex(points,d)
    affine = affineCoords([point+[1.] for point in reference])
    points = [(affine(point),point) for point in points]
    out = [max(points,key=lambda x: eval(x[0])[k]) for k in range(d+1)]
    return TRANS(out)[1]


points = [[0.485269, 0.395176, 0.114873], [0.453551, 0.624347, 1.129176], [0.07109, 0.302641, 1.321738], [1.799656, 1.419034, 1.555227], [0.184447, 1.855099, 0.702753]]

greatSimplex(points)
