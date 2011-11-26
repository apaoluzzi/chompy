from chompy import *


def Boundary(obj,d=-1):
    """
    Compute the boundary function 'c_d \mapsto \{ c_{d-1} \} for all the cells
    in the d-skeleton.

    Return an array of lists, indexed by the d-cells. The default is the last
    value of d.
    """
    assert(type(obj) == SimplicialComplex or type(obj) == PolytopalComplex)
    out = [list() for k in range(len(obj.cells[d]))]
    for k,val in obj.homology[d]:
        out[k].append(val)
    return out


def Coboundary(obj,d=-2):
    """
    Compute the coboundary function 'c_d \mapsto \{ c_{d+1} \} for all the cells
    in the d-skeleton.

    Return an array of lists, indexed by the d-cells. The default is the
    penultimate value of d, so that the function returns subsets of the
    last d-skeleton.
    """
    assert(type(obj) == SimplicialComplex or type(obj) == PolytopalComplex)
    out = [list() for k in range(len(obj.cells[d]))]
    for k,val in obj.homology[d+1]:
        out[val].append(k)
    return out


def boundaryEdgeFace(obj):
    """
    Compute the subset of incidence relation edge-face restricted to the boundary.

    Return a list of lists, indexed by the edges.
    """
    assert(type(obj) == SimplicialComplex or type(obj) == PolytopalComplex)
    coboundary_2 = Coboundary(obj)
    coboundary_1 = Coboundary(obj,-3)
    edgeFace = [[face for face in coboundary_1[edge]
                 if len(coboundary_2[face])==1]
                for edge in range(len(obj.cells[-3]))]
    return edgeFace



def decWingedRepr(obj):
    """
    Compute the **decompositive** *winged representation* of the
    SimplicialComplex 'obj'.
    (see Paoluzzi et al., Dimension-independent modeling with simplicial
    complexes, ACM Transactions on Graphics, 12 (1993), page 63).

    Return a list of triples:
    '[vertices, facets, adjCells]' indexed by d-cell indices.
    where:
    'vertices'  = list of (indices of) vertices of d-cell;
    'facets'    = list of (indices of) incident (d-1)-cells of d-cell;
    'adjCells'  = list of (indices of) adjacent d-cells of d-cell.
    
    A value '-1' in 'adjCells' denotes "no adjacency" in that position, i.e.
    the position of a *boundary facet*.
    """
    assert(type(obj) == SimplicialComplex or type(obj) == PolytopalComplex)
    boundary_d = Boundary(obj,-1)
    coboundary_dminus1 = Coboundary(obj,-2)
    def cellSelect(cellList,cell):
        if len(cellList)==1: out = -1
        elif cellList[0]==cell: out = cellList[1]
        else: out = cellList[0]
        return out
    tuples = []
    for cell,facets in enumerate(boundary_d):
        adjCells = [cellSelect(coboundary_dminus1[facet],cell) for facet in facets]
        vertices = obj.cells[-1][cell]
        tuples += [[vertices, facets, adjCells]]
    return tuples


def boundaryWingedRepr(obj):
    """
    Compute the **boundary** *winged representation* of the
    SimplicialComplex 'obj'.
    (see Paoluzzi et al., Dimension-independent modeling with simplicial
    complexes, ACM Transactions on Graphics, 12 (1993), page 63).

    Return a list of triples:
    '[vertices, facets, adjCells]' indexed by d-cell indices.
    where:
    'vertices'  = list of (indices of) vertices of (d-1)-cell;
    'facets'    = list of (indices of) incident (d-2)-cells of (d-1)-cell;
    'adjCells'  = list of (indices of) adjacent (d-1)-cells of (d-1)-cell.
    
    A value '-1' in 'adjCells' denotes "no adjacency" in that position, i.e.
    the position of a *boundary facet*. Such values are not allowed if the
    boundary is closed (i.e. without boundary, like the external "shells"
    of a solid)
    """
    assert(type(obj) == SimplicialComplex or type(obj) == PolytopalComplex)
    boundary_dminus1 = Boundary(obj,-2)
    coboundary_dminus1 = Coboundary(obj,-2)
    edgeFace = boundaryEdgeFace(obj)
    
    tuples = []
    for cell,cocells in enumerate(coboundary_dminus1):        
        vertices = obj.cells[-2][cell]
        facets = boundary_dminus1[cell]
        adjCells = []
        if len(cocells) == 1:
            for facet in facets:
                adjCells += [f for f in edgeFace[facet] if f !=cell]
            tuples += [[vertices, facets, adjCells]]
    return tuples


def orientedBoundaryWingedRepr(obj):
    """
    Compute the **oriented** **boundary winged** representation of the
    SimplicialComplex 'obj'.
    (see Paoluzzi et al., Dimension-independent modeling with simplicial
    complexes, ACM Transactions on Graphics, 12 (1993), page 63).

    Return a list of triples:
    '[vertices, facets, adjCells]' indexed by d-cell indices.
    where:
    'vertices'  = list of (indices of) vertices of (d-1)-cell;
    'facets'    = list of (indices of) incident (d-2)-cells of (d-1)-cell;
    'adjCells'  = list of (indices of) adjacent (d-1)-cells of (d-1)-cell.
    
    A value '-1' in 'adjCells' denotes "no adjacency" in that position, i.e.
    the position of a *boundary facet*. Such values are not allowed if the
    boundary is closed (i.e. without boundary, like the external "shells"
    of a solid)
    """

    def boundaryEdgeFace(obj):
        coboundary_1 = Coboundary(obj,-2)
        edgeFace = [[face for face in coboundary_1[edge]]
                    for edge in range(len(obj.cells[-2]))]
        return edgeFace

    assert(type(obj) == SimplicialComplex or type(obj) == PolytopalComplex)
    boundary_dminus1 = Boundary(obj)
    coboundary_dminus1 = Coboundary(obj,-2)
    edgeFace = boundaryEdgeFace(obj)
    
    tuples = []
    for cell,vertices in enumerate(obj.cells[-1]):
        edges = boundary_dminus1[cell]
        facets = range(len(obj.cells[-1]))
        adjCells = [[f for f in edgeFace[edge] if f != cell] for edge in edges] 
        tuples += [[vertices, edges, CAT(adjCells)]]
    return tuples


def setUpChains(listOfLists,d):
    """
    To set up a single tuple of the winged representation of a complex
    to be displayed by the 'draw' function.

    Usage:

    >>> obj = simplexGrid([2*[1],2*[1],[1]])
    >>> objRepr = decWingedRepr(obj)
    >>> for k in range(len(objRepr)):
            draw(obj,chains=setUpChains(objRepr[k],3))

    Return a list of k-chains, for k in range(d+1).
    """
    chains = (d+1)*[[]]
    chains[0] = listOfLists[0]
    for k in range(1,d-2): chains[k]=[]
    chains[d-1] = listOfLists[1]
    chains[d] = [el for el in listOfLists[2] if el != -1]
    myprint("chains",chains)
    return chains


if __name__=="__main__":
    obj = simplexGrid([2*[1],2*[1],[1]])
    obj = cuboid([1,1,1])
    obj = cprod([intervals(2.0)(2), intervals(2.0)(2), intervals(2.0)(2)])
    myprint("decWingedRepr(obj)",decWingedRepr(obj))
    myprint("boundaryWingedRepr(obj)",boundaryWingedRepr(obj))


if False:
    obj = cuboid([1,1,1])
    objRepr = orientedBoundaryWingedRepr(obj.boundary())
    obj = simplexGrid([2*[1],2*[1],[1]]).boundary()
    objRepr = orientedBoundaryWingedRepr(obj)
    myprint("objRepr",objRepr)

if False:
    def setUpChains(listOfLists,d):
        chains = (d+1)*[[]]
        chains[0] = listOfLists[0]
        for k in range(1,d-2): chains[k]=[]
        chains[d-1] = listOfLists[1]
        chains[d] = [el for el in listOfLists[2] if el != -1]
        myprint("chains",chains)
        return chains
    
    for k in range(len(objRepr)):
        draw(obj,chains=setUpChains(objRepr[k],2))
        

