from chompy import *


c = simplexGrid([2*[1],2*[1],[1]])
c.view()

print c.II(0,0,0)
print c.III(0,0,0)
print c.III(1,0,0)/c.III(0,0,0)


c = simplexGrid([2*[1],2*[1]]).embed(1)
c.view()

print c.II(0,0,0)
print c.II(1,0,0)/c.II(0,0,0)
    
