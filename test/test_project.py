from chompy import *
""" Projection of cellular complexes and point sets """

points = cuboid([1,1,1]).vertices
myprint("points",points)
points.project(1)
myprint("points.project(1)",points.project(1))

    

cuboid([1,1,1]).view()

a = cuboid([1,1,1]).project()
draw(a)

simplex(3).view()
a = simplex(3).project()
draw(a)



a = simplexGrid([4*[1],3*[1],2*[1]])
a.view()
a = a.project()
a.view()
draw(a)
