from chompy import *


def Boundary(obj,d=-1):
    assert(type(obj) == SimplicialComplex)
    out = [list() for k in range(len(obj.cells[d]))]
    for k,val in obj.homology[d]:
        out[k].append(val)
    return out


def Coboundary(obj,d=-2):
    assert(type(obj) == SimplicialComplex)
    out = [list() for k in range(len(obj.cells[d]))]
    for k,val in obj.homology[d+1]:
        out[val].append(k)
    return out


def boundaryEdgeFace(obj):
    assert(type(obj) == SimplicialComplex)
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
    assert(type(obj) == SimplicialComplex)
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
    assert(type(obj) == SimplicialComplex)
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



if __name__=="__main__":
    obj = simplexGrid([[1],[1],[1]])
    myprint("decWingedRepr(obj)",decWingedRepr(obj))
    myprint("boundaryWingedRepr(obj)",boundaryWingedRepr(obj))

