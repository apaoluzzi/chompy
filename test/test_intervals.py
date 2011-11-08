from chompy import *


domain = intervals(1)(10)
cprod([domain,domain]).view()


graph(intervals(2*pi)(20))([ID,COS]).view()
draw(graph(intervals(2*PI)(20))([ID,COS]))


graph(intervals(2*pi)(20))([SIN,COS]).view()
draw(graph(intervals(2*PI)(20))([SIN,COS]))
draw(circumpherence(0.5))


graph(intervals(4*PI)(40))([SIN,COS,ID]).view()
draw(graph(intervals(4*PI)(40))([SIN,COS,RAISE(DIV)([ID,K(4*PI)])]))
draw(helix())
draw(helix(radius=0.5,pitch=1./4,n=12,turns=4))

