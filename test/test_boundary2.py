from chompy import *

c = simplexGrid([2*[1.0],4*[1.0],1*[1.0]])
          

c.view(2)
c.boundary().view(2)
draw(c.boundary())
        
