from chompy import *



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



if __name__=="__main__":
    obj = simplexGrid([[1],[1],[1]])
    obj = cuboid([1,1,1])
    objRepr = orientedBoundaryWingedRepr(obj.boundary())
    myprint("objRepr",objRepr)
    for k in range(len(objRepr)):
        draw(obj.boundary(),chains=objRepr[k])
        
