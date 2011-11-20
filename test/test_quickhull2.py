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
                linalg.det(mat(basis(cell))) > 1.0:
                cell = []
        return cell

    def mostRemote(points,k):
        distance,pivot = 0.0,-1
        for index,item in enumerate(points):
            affine,coords = item
            myprint("k,affine,coords",(k,affine,coords))
            if affine[k] < distance:
                distance = affine[k]
                pivot = index
                myprint("pivot",pivot)
        return pivot
    
    """
        1. Preprocessing: find any d+1 affinely independent points.
        Build one initial simplex. 
    """
    inputPoints = points
    m,d = len(points),len(points[0])
    cell = [10, 0, 21] #randomSimplex(m,d)
    draw(SimplicialComplex(inputPoints,[cell]))

    """
        2. Compute the barycentric coordinates of points.
        Sort on such coordinates.
        Remove the points contained in the current simplex.
    """
    def reframe(points,cell):
        newCoords = affineCoords(simplex(cell))
        items = sorted([[CONS([grading,ID])(newCoords(point)),point]
                          for point in points],reverse=True)
        myprint("items",items)
        points = [(eval(affine),coords)
                  for [grade,affine],coords in items if grade<=d+1]
        return points,cell

    points,cell = reframe(points,cell)

    """
        3. Enqueue an adjacent simplex, and suitably update
        the winged adjaciency.
    """
    
    cellStore = []
    lastCell = CCI = 0 ## Current Cell Index
    voids = [-1 for k in range(d+1)]
    wingedCell = [cell,voids]
    cellStore.append(wingedCell)
    myprint("cellStore",cellStore)

    def exploreSimplex(cellStore,lastCell,CCI):
        theCell,theAdjacencies = cellStore[CCI]
        voids = [-1 for k in range(d+1)]
        for k,adj in enumerate(theAdjacencies):
            if adj == -1:
                myprint("k",k)
                pivotRef = mostRemote(points,k)
                
                if pivotRef != -1:
                    newCell = copy.copy(theCell)
                    newCell[k] = pivotRef
                    newAdjacencies = copy.copy(voids)
                    newAdjacencies[k] = CCI
                    lastCell += 1
                    theAdjacencies[k] = lastCell
                    cellStore.append([newCell,newAdjacencies])
                    myprint("cellStore",cellStore)
        return cellStore,lastCell,CCI

    while True:
        cellStore,lastCell,CCI = exploreSimplex(cellStore,lastCell,CCI)
        draw(SimplicialComplex(inputPoints,TRANS(cellStore)[0]))
        CCI += 1
        if CCI == lastCell+1: break

    return(cellStore)



#dataSet = simplexGrid([3*[1],4*[1],2*[1]])
#dataSet = simplexGrid([[1],[1],[1]])
#dataSet = simplexGrid([[1],[1]])
dataSet = simplexGrid([5*[1],5*[1]])

#points = dataSet.vertices.rotate(1,3,pi/6).rotate(2,3,pi/4).points
points = dataSet.vertices.points
myprint("points",points)
cellStore = quickhull(points)
myprint("points",points)
myprint("cellStore",cellStore)
