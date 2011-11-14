from chompy import *

from scipy.spatial.qhull import *

help(scipy.spatial.qhull)

points = [[0.0, 1.414214], [1.207107, -0.707107], [1.207107, 0.707107],
       [-1.207107, -0.707107], [0.0, -1.414214], [-1.207107, 0.707107]]

Delaunay(points).convex_hull

mycycle = [4, 3, 5, 0, 2, 1, 4]
pol = polyline([points[k] for k in mycycle])

myprint("pol.vertices",pol.vertices)
myprint("pol.cells",pol.cells)

pol.view()
