from chompy import *


c = cuboid([1,2,3])

c.view()
c.boundary().view()
c.boundary().view(2)
draw(c.boundary(),expl=[2,2,2])



c = cprod([intervals(1)(3),
           intervals(2)(3),
           intervals(3)(3)])

draw(c,expl=[2,2,2])
draw(c.boundary(),expl=[2,2,2])
