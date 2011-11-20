from chompy import *


def grading(point):
    def grade(x):
        if x >= 0.0: return'1'
        else: return '0'
    return "0B"+"".join(AA(grade)(eval(point)))


def affineCoords(basis):
    d = len(basis)
    mapping = mat([vect for vect in basis]).I
    def affineCoords0(point):
        affine = (mat(point+[1.])*mapping).tolist()[0]
        out = code(affine)
        return out
    return affineCoords0


def randomPoints(d=2,n=40,scale=2):
    return PointSet([[random.random()*scale for k in range(d)]
                     for point in range(n)])

def simplex(verts,cell):
    """ simplex square matrix in homogeneous coordinates """
    return [array(verts[k] + [1.]) for k in cell]


def greatSimplex(points):
    def randomSimplex(verts,d):
        cell,m = [],2*d
        while len(cell) < d+1:
            index = int(random.random()*m)
            if index not in cell:
                cell += [index]
            if len(cell) == d+1 and \
                linalg.det(mat(simplex(verts,cell))) < 1:
                cell = []
        return [verts[k] for k in cell]   
    d = len(points[0])
    extremes = CAT([(max(points,key=lambda x: x[k]),
                     min(points,key=lambda x: x[k])) for k in range(d)])
    return randomSimplex(extremes,d)

def reframe(points,cell):
    newCoords = affineCoords(simplex(points,cell))
    items = sorted([[CONS([grading,ID])(newCoords(point)),point]
                      for point in points],reverse=True)
    dpoints = [(CONS([eval,ID])(grade),eval(affine),coords)
               for [grade,affine],coords in items]
    return dpoints




if __name__=="__main__":

    rpoints = randomPoints(d=3,n=40,scale=2)
    points = rpoints.points
    myprint("points",points)
    draw(polyline(points),chains=[range(40)])    

    cell = [rpoints.dict[code(p)] for p in greatSimplex(points)]
    draw(SimplicialComplex(points,[cell]),colors=[RED,GREEN,BLUE])

    dpoints = reframe(points,cell)
    myprint("dpoints",dpoints)

