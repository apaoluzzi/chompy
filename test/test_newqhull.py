from chompy import *
"""
A novel algorithm for computing the convex hull of d-dimensional points.
"""


def grading(point):
    """
    To compute the "grade" of a point in affine coordinates.

    Return a binary number (d+1 bits) that gives a classification of a
    d-point in affine coordinates in the proper "tiling" of d-space.
    """
    def grade(x):
        if x >= 0.0: return'1'
        else: return '0'
    return "0B"+"".join(AA(grade)(eval(point)))


def affineCoords(basis):
    """
    To compute the affine coordinates of a point w.r.t. a reference simplex.

    Return the coded string of affine coordinates.
    """
    d = len(basis)
    mapping = mat([vect for vect in basis]).I
    def affineCoords0(point):
        affine = (mat(point+[1.])*mapping).tolist()[0]
        return code(affine)
    return affineCoords0



def randomPoints(d=2,n=40,scale=2):
    """
    To produce a random PointSet in [0,scale]^d space.

    Return a PointSet instance with n points, each with d coordinates.
    """
    return PointSet([[random.random()*scale for k in range(d)]
                     for point in range(n)])

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
            print "eccomi"
            index = int(random.random()*(m+1))
            if index not in cell: cell += [index]
            if len(cell) == d+1 and \
                closetozero(linalg.det(mat(simplex(verts,cell)))):
                cell = []
        return [verts[k] for k in cell]
    
    d = len(points[0])
    extremes = CAT([(max(points,key=lambda x: x[k]),
                     min(points,key=lambda x: x[k])) for k in range(d)])
    return randomSimplex(extremes,d)



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

    d = len(points[0])
    reference = randomSimplex(points,d)
    affine = affineCoords([point+[1.] for point in reference])
    points = [(affine(point),point) for point in points]
    out = [max(points,key=lambda x: eval(x[0])[k]) for k in range(d+1)]
    out = TRANS(out)[1]
    if len(out) == len(list(set(AA(code)(out)))): return out
    else: return TRANS(points)[1][:d+1]
    
 

def reframe(points,cell):
    """
    Sorted classification of d-points w.r.t. a reference 'cell' (d+1 simplex). 

    Return a (hierarchical) dictionary of disjoint regions.
    """
    newCoords = affineCoords(simplex(points,cell))
    items = sorted([[CONS([grading,ID])(newCoords(point)),point]
                      for point in points],reverse=True)
    
    """ dictionary of tiled points construction"""
    def regionKey(grade): return tuple(CONS((eval,ID))(grade))
    def regionKey(grade): return grade
    def doubleCoords(affine,coords): return eval(affine),coords
    dpoints = [(regionKey(grade),doubleCoords(affine,coords))
               for [grade,affine],coords in items]
    regionDict = {}
    for point in dpoints:
        thecode,value = point
        if regionDict.has_key(thecode): regionDict[thecode].append(value)
        else: regionDict[thecode] = [value]
    """ recursive contruction of dictionaries in crowded subregions"""
    return regionDict


def selectBasis(points):
    """
    Select the basis elements from the "points" data structure
    (list of pairs (Cartesian, affine) of point coordinates).

    Return the Cartesian representation of the basis elements (in affine
    coordinates).
    """
    def inBasis(affinePoint):
        return OR([True if coord==1.0 else False for coord in affinePoint])
    
    return [point for point in points if inBasis(point[0])]


def pointClassify(pointset):
    """
    To compute the classification of a PointSet w.r.t. to a maximal
    subset of affinely independent elements (greatSimplex).

    Return the dictionary of non-empty regions of the tiling.
    Each tile has a binary key, and as value the subset of points
    contained within.
    """
    points = pointset.points
    simplexPoints = greatSimplex(points)
    cell = [pointset.dict[code(p)] for p in simplexPoints]
    theRegionDict = reframe(points,cell)
    return theRegionDict


def drawRegions(regionDict):
    """
    To draw stepwise a region dictionary with different colors.

    No return value. Just display side effects.
    """
    
    thecolors=[CYAN,MAGENTA,WHITE,YELLOW,RED,
               GREEN,BLUE,GRAY,BROWN,BLACK,ORANGE,PURPLE]
    col=0
    for region in regionDict:
        if type(regionDict[region]) == list:
            
            pointSubset = [point[1] for point in regionDict[region]]
            d = len(pointSubset[0])
            
            if eval(region) == (2**(d+1) - 1):
                simplex = selectBasis(regionDict[region])
                out = SimplicialComplex(simplex,[range(len(simplex))])
            else:
                if 1 <= len(pointSubset) <= d+1:
                    out = SimplicialComplex(pointSubset,[range(len(pointSubset))])
                else:
                    out = SimplicialComplex(pointSubset,AA(LIST)(range(len(pointSubset))))
                    
            draw(out,color=col%len(thecolors),colors=thecolors)
        col += 1


def traverseDict(dictionary):
    """
    Scheme of DFS traversal of a (recursive) dictionary of dictionaries.
    """
    def elaborate (value):
        pass
    
    for key in dictionary.keys():
        if key == special:
            # do something
            pass
        elif type(dictionary[key]) == dict:
            traverseDict(dictionary[key])
            # do elaborate value
        elaborate(dictionary[key])
    return None



def makeRegionDict(pointSet,d):
    regionDict = pointClassify(pointSet)
    for key in regionDict.keys():
        if len(regionDict[key]) > d+1:
            if eval(key) == (2**(d+1) - 1):
                regionDict[key] = selectBasis(regionDict[key])
            else:
                pointSubset = [point[1] for point in regionDict[key]]
                regionDict[key] = makeRegionDict(PointSet(pointSubset),d)
        else: pass
    return regionDict



def drawRegionDict(dictionary):
    """
    Scheme of DFS traversal of a (recursive) dictionary of dictionaries.
    """
    def elaborate (value):
        points = TRANS(value)[1]
        cells = [range(len(points))]
        draw(SimplicialComplex(points,cells))
    
    for key in dictionary.keys():
        if type(dictionary[key]) == dict:
            drawRegionDict(dictionary[key])
        else: elaborate(dictionary[key])
    return None


if __name__=="__main__":

    rpoints = randomPoints(d=3,n=40,scale=2)
    rpoints = randomPoints(d=2,n=40,scale=3)
    d = rpoints.rn
    points = rpoints.points
    draw(polyline(points),chains=[range(len(polyline(points).cells[0]))])

    regionDict = makeRegionDict(rpoints,d)
    myprint("regionDict",regionDict)
    drawRegionDict(regionDict)
            
    
