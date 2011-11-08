from chompy import *
    

c = torus_solid(#n=16,m=24
    )

c.view()
c.view(2)
d = c.boundary()
d.view()
d.view(2)
draw(d)
draw(d.boundary())

