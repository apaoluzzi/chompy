from chompy import *


c = torus_surface()

c.view()
c.view(2)
c.boundary().view()

draw(c)
draw(c.boundary())


