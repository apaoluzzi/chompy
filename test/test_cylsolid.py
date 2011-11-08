from chompy import *
    

c = cylsolid(1/SQRT(2),1,n=4)

c.view()
c.view(2)
c.view(3)
c.boundary().view()
c.boundary().view(2)

draw(c)
draw(c.boundary(),expl=[2,2,2])
draw(c,expl=[2,2,2],chains=[[],[],[],range(len(c.cells[3]))] )
draw(c.boundary(),expl=[2,2,2],chains=[[],[],range(len(c.boundary().cells[2]))] )

print c.boundary().II(0,0,0)
print c.boundary().III(0,0,0)

        
    
