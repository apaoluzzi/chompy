from chompy import *

"""
Function to generate a simplicial complex triangulating the
convex hull of a set of n-dimensional points.

Algorithm: QuickHull + Winged representation of polytopes
---------------------------------------------------------

0. Preprocessing: find any n+1 affinely independent points.
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

    inputPoints = copy.copy(points)

    def localCoords(basis):
        d = len(basis)
        mapping = mat([vect for vect in basis]).I
        def localCoords0(point):
            affine = (mat(point+[1.])*mapping).tolist()[0]
            standard = (array(affine[:-1])/(d+1)).tolist() + [1-sum(affine[:-1])] 
            return code(affine)
        return localCoords0
    

    def basis(cell):
        return [(array(points[k])-cell[0]).tolist() for k in cell[1:]]

    def simplex(cell):
        """ simplex square matrix in homogeneous coordinates """
        return [array(points[k] + [1.]) for k in cell]
    

    def grading(point):
        return sum(AA(SIGN)(eval(point)))
    """
        1. Preprocessing: find any n+1 affinely independent points.
        Construct one initial simplex. Set it as the current one.
    """
    m,n = len(points),len(points[0])
    points = sorted(points,reverse=True)
    cell = [0, m/3,2*m/3, m-1]
    if linalg.det(mat(basis(cell))) == 0.0: cell[0] = 1

    myprint("cell",cell)
    draw(SimplicialComplex(points,[cell]))

    """
        2. Compute the barycentric coordinates of points.
        Sort on such coordinates.
        Remove the points contained in the current simplex.
    """
    newCoords = localCoords(simplex(cell))
    items = sorted([[CONS([grading,ID])(newCoords(point)),point]
                      for point in points],reverse=True)
    points = [(eval(affine),coords)
              for [grade,affine],coords in items if grade<4]
    myprint("points",points)

    
    pol = polyline(TRANS(points)[1])
    draw(pol,[range(len(pol.cells[0]))])

    """
        3. Enqueue an adjacent simplex, and update suitably
        the winged adjaciency.
    """
    def mostDistant(points,k):
        distance = 1
        pivot = -1
        for index,item in enumerate(points):
            affine,coords = item
            if affine[k] < distance:
                distance = affine[k]
                pivot = index           
        return pivot
        
    cellStore = []
    currentCell = 0
    lastCell = currentCell
    emptyNeighbours = [-1 for k in range(n+1)]
    wingedCell = [cell,emptyNeighbours]
    cellStore.append(wingedCell)
    myprint("cellStore",cellStore)
    theCell,theAdjacencies = cellStore[currentCell]
    for k in range(len(theAdjacencies)):
        if theAdjacencies[k] == -1:
            pivotPoint = mostDistant(points,k)
            newCell = copy.copy(theCell)
            newCell[k] = pivotPoint
            cellStore.append([newCell,emptyNeighbours])
            lastCell += 1
            theAdjacencies[k] = lastCell

            myprint("inputPoints,TRANS(cellStore)[0]",
                    (inputPoints,TRANS(cellStore)[0]))
            draw(SimplicialComplex(inputPoints,TRANS(cellStore)[0]))

    

    return(points)



dataSet = simplexGrid([3*[1],4*[1],2*[1]])
dataSet = simplexGrid([[1],[1],[1]])

points = dataSet.vertices.rotate(1,3,pi/6).rotate(2,3,pi/4).points
points = quickhull(points)
myprint("points",points)
