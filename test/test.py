from chompy import *

a = simplexGrid([3*[1], 3*[1], 1*[1]])

a = a.rotate(0,1,PI/6)
print a.vertices.dict
a.view()
a.view(2)
a.view(3)

draw(a)
draw(a,expl=[1.5,1.5,1.5])
