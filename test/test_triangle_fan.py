from chompy import *
from scipy.spatial.qhull import *

points = [[cos(alpha),sin(alpha)] for alpha in scipy.linspace(0.0,2*pi,7)]
myprint("points",points)
points = (array(points) + [2.,1.]).tolist()
myprint("translated points",points)

trianglefan(points).view()
trianglefan(points).view(2)
trianglefan(points).boundary().view()
