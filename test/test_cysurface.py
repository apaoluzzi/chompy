from chompy import *


c = cylsurface(0.5,3,m=4)

c.view()
c.view(2)
c.boundary().view()

draw(c)
draw(c.boundary())

