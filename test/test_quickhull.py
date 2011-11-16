from chompy import *

"""
Function to generate a simplicial complex triangulating the
convex hull of a set of d-dimensional points.

Algorithm: QuickHull + Winged representation of polytopes
---------------------------------------------------------

0. Preprocessing: find any d+1 affinely independent points.
1. Construct one initial simplex. Set it as the current one.
2. Compute the barycentric coordinates of points. Sort on the coordinates.
3. Remove the points contained in the current simplex.
4. Enqueue an adjacent simplex, and mark suitably the winged adjaciency.
5. Continue until either all the adjaciencies of the current simplex are
marked or no points exist with the corresponding coordinate < 0.
6. remove the top from the queue and set the next simplex as current.
7. repeat from point 2 until all points are gone.
"""


def quickhull(points):

    def affineCoords(basis):
        d = len(basis)
        mapping = mat([vect for vect in basis]).I
        def affineCoords0(point):
            affine = (mat(point+[1.])*mapping).tolist()[0]
            out = code(affine)
            return out
        return affineCoords0
    
    def basis(cell):
        return [(array(points[k])-cell[0]).tolist() for k in cell[1:]]

    def simplex(cell):
        """ simplex square matrix in homogeneous coordinates """
        return [array(points[k] + [1.]) for k in cell]
    
    def grading(point):
        return sum(AA(SIGN)(eval(point)))

    def randomSimplex(m,d):
        cell = []
        while len(cell) < d+1:
            index = int(random.random()*m)
            if index not in cell:
                cell += [index]
            if len(cell) == d+1 and \
                linalg.det(mat(basis(cell))) == 0.0:
                cell = []
        return cell

    def mostRemote(points,k):
        distance,pivot = 1,-1
        for index,item in enumerate(points):
            affine,coords = item
            if affine[k] < distance:
                distance = affine[k]
                pivot = index           
        return pivot
    
    """
        1. Preprocessing: find any d+1 affinely independent points.
        Build one initial simplex. 
    """
    inputPoints = points
    m,d = len(points),len(points[0])
    cell = randomSimplex(m,d)
    draw(SimplicialComplex(points,[cell]))

    """
        2. Compute the barycentric coordinates of points.
        Sort on such coordinates.
        Remove the points contained in the current simplex.
    """
    newCoords = affineCoords(simplex(cell))
    items = [[CONS([grading,ID])(newCoords(point)),point]
                      for point in points]
    points = [(eval(affine),coords)
              for [grade,affine],coords in items]

    """
        3. Enqueue an adjacent simplex, and update suitably
        the winged adjaciency.
    """
    
    cellStore = []
    lastCell = CCI = 0 ## Current Cell Index
    voids = [-1 for k in range(d+1)]
    wingedCell = [cell,voids]
    cellStore.append(wingedCell)
    
    theCell,theAdjacencies = cellStore[CCI]
    for k,adj in enumerate(theAdjacencies):
        if adj == -1:
            pivotRef = mostRemote(points,k)
            if pivotRef != -1:
                newCell = copy.copy(theCell)
                newCell[k] = pivotRef
                newAdjacencies = voids
                newAdjacencies[k] = CCI
                cellStore.append([newCell,newAdjacencies])
                lastCell += 1
                adj = lastCell
    
            draw(SimplicialComplex(inputPoints,TRANS(cellStore)[0]))
    CCI += 1

    

    return(cellStore)



dataSet = simplexGrid([3*[1],4*[1],2*[1]])
#dataSet = simplexGrid([[1],[1],[1]])
#dataSet = simplexGrid([[1],[1]])

points = dataSet.vertices.rotate(1,3,pi/6).rotate(2,3,pi/4).points
#points = dataSet.vertices.points
cellStore = quickhull(points)
myprint("points",points)
myprint("cellStore",cellStore)
