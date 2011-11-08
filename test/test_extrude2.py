from chompy import *

points = [[0,3],[1,2],[3,3],[2,2],[3,0],[2,1],[0,0],[1,1],[0,3],[1,2]]
pol = triangle_strip(points).extrude([1,-1,1]).boundary()
draw(pol, expl=[1,1,3])#,
     #chains=[[],[],range(len(pol.cells[2]))])
