from chompy import *

circle = circle2d()
              
circle.view()
circle.view(2)
circle.boundary().translate([1.,2.]).view()

t([0,1])([0.,1.])(s([0,1])([.5,.5])(circle.boundary())).view()
