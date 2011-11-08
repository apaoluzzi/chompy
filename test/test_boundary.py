from chompy import *
    


c = cylsurface(0.5,3,m=4)

c.view()
c.view(2)
c.boundary().view()

draw(c)
draw(c.boundary())

c = cylsolid(0.5,3)

c.view()
c.view(2)
c.view(3)
c.boundary().view()

draw(c)
draw(c.boundary())
    
