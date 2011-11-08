from chompy import *

c = cuboid([1,2,3])
c = CuboidalComplex([[1],[2],[3]]) 
c.boundary().view(2)

draw(c.boundary())


##cells = [[0,2,4,6],[0,1,2,4],[2,3,4,6],[0,4,5,6],[0,2,6,7]]
##verts = [[0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0],
##             [1.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0],
##             [0.0, 1.0, 0.0], [0.0, 1.0, 1.0]]
##SimplicialComplex(verts,cells).view(2)
##SimplicialComplex(verts,cells).view(3)
##draw(SimplicialComplex(verts,cells),expl=[3,3,3])



##
##cells = [[0, 1, 2, 3], [1, 2, 3, 4], [2, 3, 4, 5],
##         [3, 4, 5, 6], [4, 5, 6, 7], [5, 6, 7, 0]]
##myprint("cells",cells)
####cells = [[cell[1]]+[cell[0]]+cell[2:] for cell in cells]
####myprint("cells",cells)
##
##verts = c.vertices.points
##cubo = SimplicialComplex(verts,cells)
##draw(cubo.boundary())
###cubo = cubo.transform([[1,0,0],[1./2,1./2,0],[1./3,0,1./3]])
###draw(cubo.boundary())
##cubo.boundary().boundary()
##myprint("cubo.III(0,0,0)",cubo.III(0,0,0))
##myprint("cubo.II(0,0,0)",cubo.II(0,0,0))
##
##simplices = [SimplicialComplex(verts,[cell]) for cell in cells]
##volume = sum([simplex.III(0,0,0) for simplex in simplices])
##myprint("volume cube",volume)
##volume = sum([simplex.transform([[1,0,0],[1./2,1./2,0],[1./3,0,1./3]]).III(0,0,0) for simplex in simplices])
##myprint("volume transformed cube",volume)
##volume = [simplex.boundary().transform([[1,0,0],[1./2,1./2,0],[1./3,0,1./3]]).III(0,0,0) for simplex in simplices]
##myprint("volume boundary transformed cube",volume)
##
##
