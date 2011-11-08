from chompy import *

c0 = SimplicialComplex([[0.,0.],[1.,0.],[0.,1.],[1.,1.]],[[0,1,2],[1,3,2]])
c1 = SimplicialComplex([[0.0],[1.0],[2.0]],[[0,1],[1,2]])
c2 = SimplicialComplex([[0.0],[1.0],[1.5],[2.5],[4.0]],[[0,1],[1,2],[3,4]])


cc1 = complexProd([c1,c2])     
cc1.view()
cc1.view(2)
#draw(cc1)
cc1.boundary().view()
#draw(cc1.boundary())  DEBUG !!!  impedisce visualizzazioni successive ... !!

cc2 = complexProd([cc1,c2])
cc2.view()
cc2.view(2)
cc2.boundary().view(2)
cc2.boundary().view()
draw(cc2.boundary())

cc3 = complexProd([c0,c2])
cc3.view()
cc3.view(2)
cc3.view(3)
cc3.boundary().view(2)

draw(cc3.boundary(),expl=[1.5,1.5,3])
