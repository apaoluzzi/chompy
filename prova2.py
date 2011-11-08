from chompy import *


def triangle_array(m,n,points):
    """
    Primitive constructor of a 2-indexed array of points.
    
    """
    out = simplexGrid([m*[1.],n*[1.]])
    return SimplicialComplex(CAT(points),out.cells[2])

#/////////////////////////////////////////////////////////

if __name__ == "__main__":

    points = [ 
        [[0,0,0.1],[1,0,-0.1],[2,0,0.0],[3,0,0.2]], 
        [[0,1,-0.4],[1,1,0.1],[2,1,-0.1],[3,1,0.1]], 
        [[0,2,0.1],[1,2,0.0],[2,2,0.1],[3,2,0.1]], 
        [[0,3,-0.2],[1,3,0.1],[2,3,-0.3],[3,3,0.1]]
    ]

    ta = triangle_array(3,3,points)
    ta.view()
    tb = ta.boundary()
    tb.view()
    draw(tb)


draw(tb,
     chains=[(list(set(CAT(tb.cells[1])))),
             range(len(tb.cells[1]))])
