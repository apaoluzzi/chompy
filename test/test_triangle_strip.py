from chompy import *


p = triangle_strip([[0,0],[1,0.1],[0.1,1],[1,1]]).translate([0.1])
p.view()
p.view(2)
p.boundary().view()

points = [[0,3],[1,2],[3,3],[2,2],[3,0],[2,1],[0,0],[1,1],[0,3],[1,2]]
triangle_strip(points).view()
triangle_strip(points).view(2)
triangle_strip(points).boundary().view()

draw(triangle_strip(points).boundary().extrude([1,1,1]))
draw(triangle_strip(points).extrude([1,1,1]).boundary())
