from chompy import *


a = simplexGrid([4*[1.0]])
draw(a)
a = intervals(4.0)(4)
#draw(a)

b = simplexGrid([4*[1.0],4*[1.0]])
draw(b)
b = cprod([a,a])
draw(b)

c = simplexGrid([4*[1.0],4*[1.0],2*[1.0]])
c.view()
c = cprod([a,a,a])
draw(c)

#draw(c.project())
draw(c.project().project())
draw(c.project().project().project())
