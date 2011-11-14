from chompy import *
from scipy.spatial.qhull import *

points = [CONS([cos,sin])(alpha) for alpha in scipy.linspace(0.0,2*pi,7)]
myprint("points",points)

trianglefan(points).view(2)
pts = PointSet(points).translate([2,1])
points = pts.points + [pts.points[0]]
trianglefan(points).view(2)
