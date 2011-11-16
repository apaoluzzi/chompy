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
            affine = mat(point+[1.]) * mapping
            return code(affine.tolist()[0])
        return affineCoords0

    def basis(points,cell):
        return [(array(points[k])-cell[0]).tolist() for k in cell[1:]]

    def simplex(points,cell):
        """ simplex square matrix in homogeneous coordinates """
        return [array(points[k] + [1.]) for k in cell]

    def grading(point):
        return sum(AA(SIGN)(eval(point)))

    def mostDistant(points,k):
        distance, pivot = 1, -1
        for index,item in enumerate(points):
            affine,coords = item
            if affine[k] < distance:
                distance = affine[k]
                pivot = index           
        return pivot





    
    """
        1. Preprocessing: find any d+1 affinely independent points.
        Construct one initial simplex. Set it as the current one.
    """

    m,d = len(points),len(points[0])
    points = sorted(points,reverse=True)
    inputPoints = copy.copy(points)    
    myprint("inputPoints",inputPoints)
    cell = range(0, m, m/d)
    myprint("cell",cell)
    if linalg.det(mat(basis(inputPoints,cell))) == 0.0: cell[0] = 1

    myprint("cell",cell)
    if True: draw(SimplicialComplex(inputPoints,[cell]))

    """
        2. Compute the affine coordinates of points.
        Append affine coordinates to Cartesian coordinates.
        Remove the points contained in the current cell,
        i.e. the points with positive coords summing to one.
    """
    
    newCoords = affineCoords(simplex(inputPoints,cell))
    items = sorted([[CONS([grading,ID])(newCoords(point)),point]
                      for point in points],reverse=True)
    #-- remove a subset of points from current set ----------------
    #-- outPoints are external to current triangulation -----------
    outPoints = [(eval(affine),coords)
                 for [grade,affine],coords in items if grade<4]

    """
        3. Enqueue an adjacent simplex, and update suitably
        the winged adjaciency.
    """
    
    cellStore = []
    currentCell = 0
    lastCell = currentCell
    emptyNeighbours = [-1 for k in range(d+1)]
    wingedCell = [cell,emptyNeighbours]
    cellStore.append(wingedCell)
    myprint("cellStore",cellStore)
    theCell,theAdjacencies = cellStore[currentCell]
    for k in range(len(theAdjacencies)):
        if theAdjacencies[k] == -1:
            pivotPoint = mostDistant(outPoints,k)
            newCell = copy.copy(theCell)
            newCell[k] = pivotPoint
            cellStore.append([newCell,emptyNeighbours])
            lastCell += 1
            theAdjacencies[k] = lastCell

            myprint("inputPoints,TRANS(cellStore)[0]",
                    (inputPoints,TRANS(cellStore)[0]))
            if False: draw(SimplicialComplex(inputPoints,TRANS(cellStore)[0]))

    

    return(outPoints)



dataSet = simplexGrid([3*[1],4*[1],2*[1]])
dataSet = simplexGrid([[1],[1],[1]])

points = dataSet.vertices.rotate(1,3,pi/6).rotate(2,3,pi/4).points
myprint("points",points)
points = quickhull(points)
