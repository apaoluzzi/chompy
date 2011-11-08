from chompy import *

cuboid([1,1]).view()
draw(cuboid([1,1]))

cuboid([1,1,1]).view()
draw(cuboid([1,1,1]))



out = schlegel(cuboid([1,1,1]).translate([-0.5,-0.5,0.5]))
out.view()
draw(out)



out = schlegel(cuboid([1,1,1,1]).translate([-0.5,-0.5,-0.5,0.5]))
out.view()
draw(out,chains=[range(len(out.cells[0])),
                 range(len(out.cells[1]))])

