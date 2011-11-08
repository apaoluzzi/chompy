from config import *


from config import *


# /////////////////////////////////////////////
# celle dw di una cella
# /////////////////////////////////////////////
def DOWNCELLS(g):
    def DOWNCELLS0(cell):
        it=g.goDw(cell);ret=[]
        while not it.end(): ret+=[it.getNode()];it.goForward()
        return ret
    return DOWNCELLS0

# /////////////////////////////////////////////
# celle up di una cella
# /////////////////////////////////////////////
def UPCELLS(g):
    def UPCELLS0(cell):
        it=g.goUp(cell);ret=[]
        while not it.end(): ret+=[it.getNode()];it.goForward()
        return ret
    return UPCELLS0

# /////////////////////////////////////////////
#  celle ad un certo livello del grafo
# /////////////////////////////////////////////
def CELLSPERLEVEL(g):
    def CELLSPERLEVEL0(level):
        it=g.each(level) 
        ret=[]
        while not it.end():
            ret+=[it.getNode()];it.goForward()
        return REVERSE(ret)  
    return CELLSPERLEVEL0



# /////////////////////////////////////////////
# Generation of 1D polyhedron
# /////////////////////////////////////////////

def Quote(g,numList):
    """ To create the graph of a 1D cell complex """
    sizes = [abs(num) for num in numList]
    points = [Vecf([1.0, x]) for x in AL([0,PROGRESSIVESUM(sizes)])]
    for point in points: node = g.addNode(0); g.setVecf(node,point)
    nodes = CELLSPERLEVEL(g)(0)
    edges = [[nodes[k],nodes[k+1]] for k in range(len(nodes)-1)]
    for edge in edges: g.addNode(1)
    nodes = CELLSPERLEVEL(g)(1)
    for k in range(len(edges)):
        if numList[k] > 0:
            g.addArch(edges[k][0],nodes[k]);g.addArch(edges[k][1],nodes[k])
            #aggiungi la doppia connettivita' dei nodi a livello top
            g.addArch(nodes[k],edges[k][0]);g.addArch(nodes[k],edges[k][1])
        else: g.remNode(nodes[k])
    return g
    

# /////////////////////////////////////////////
# Generation of nD grids
# /////////////////////////////////////////////

def Grid(listOfListsOfNum):
    numList = listOfListsOfNum[0]
    g = Quote(Graph(1),numList)
    for numList in listOfListsOfNum[1:]:
        g1 = Quote(Graph(1),numList)
        g = Graph.power(Matf(1),Matf(1),  g,None,None,  g1,None,None)
    return g



###################################################################

# /////////////////////////////////////////////
# Drawing Graph instances with spheres and cylinders
# /////////////////////////////////////////////


def DRAW(g,expl=[1,1,1]):

    n = g.getMaxDimCells()
    m = g.getPointDim()

    def offset(point,expl=[1,1,1]):
        scaledpoint = [point[k]*expl[k] for k in range(3)]
        vect = VECTDIFF([scaledpoint,point])
        return vect
        

    def spheres(points,expl=[1,1,1]):
        batches = []
        sx = 0.05
        myprint("Batch.openObj('sphere18x27.obj')",Batch.openObj("sphere18x27.obj"))
        unitSphere = Batch.openObj("sphere18x27.obj")[0]
        for point in points:
            batchSphere = Batch(unitSphere)
            vect = offset(point,expl)
            if len(point) == 2:
                point = point + [0.0]
                vect = vect + [0.0]
            batchSphere.matrix =  Mat4f.translate(*vect) * \
                Mat4f.translate(*point)*Mat4f.scale(sx,sx,sx)
            batchSphere.diffuse=CYAN
            batches += [batchSphere]
        return batches


    def transfCylr(batchCylinder,pointpair,expl=[1,1,1]):
        vect,point = VECTDIFF(REVERSE(pointpair)),pointpair[0]
        sx = 0.025
        
        def vectTransform(vect):
            qz = UNITVECT(vect)
            qx = UNITVECT(VECTPROD([ vect,[0,0,1] ]))
            qy = VECTPROD([ qz,qx ])
            Rot = TRANS([qx,qy,qz]) 
            Rot = CAT([ Rot[0]+[0.], Rot[1]+[0.], Rot[2]+[0.], [0.,0.,0.,1.] ])
            h = VECTNORM(vect)
            
            def isclose (a,b,filter_threshold=0.5):
                if abs(a-b)<filter_threshold: return True
                else: return False

            if isclose (Mat4f.determinant(Mat4f(*Rot)),
                        0.0, 1E-5):
                return h,Mat4f.scale(1,SIGN(vect[1]),SIGN(vect[2]))
            else: return h,Mat4f(*Rot)
            
        h,rot = vectTransform(vect)
        center = [c/2.0 for c in VECTSUM(pointpair)]
        vect = offset(center,expl)
        batchCylinder.matrix = Mat4f.translate(*vect) *\
            Mat4f.translate(*point) * rot * Mat4f.scale(sx,sx,h)
        batchCylinder.diffuse = MAGENTA
        return batchCylinder

    def cylinders(batches,edgepoints,expl=[1,1,1]):
        unitCylinder = Batch.openObj("cylinder4x27.obj")[0]        
        vects = [VECTDIFF(edge) for edge in edgepoints]
        for pointpair in edgepoints:
            batchCyl = Batch(unitCylinder)
            batchCyl = transfCylr(batchCyl,pointpair,expl)
            batches += [batchCyl]
        return batches

    def planecells(batches,facepoints,expl=[1,1,1]):
        for points in facepoints:
            n = len(points)
            center = [coord/float(n) for coord in VECTSUM(points)]
            vect = offset(center,expl)
            points = [[point[k]+vect[k] for k in range(3)] for point in points]
            def sign(points):
                return SIGN(VECTPROD(AA(C(VECTDIFF)(center))(points[2:0:-1])))
            face = MKPOL([points,[range(1,n+1)],None])
            faceBatch = Plasm.getBatches(face)
            faceBatch[0].diffuse = WHITE
            batches += faceBatch
        return batches

    def cells(batches,cellpoints,expl=[1,1,1]):
        for points in cellpoints:
            n = len(points)
            center = [coord/float(n) for coord in VECTSUM(points)]
            vect = offset(center,expl)
            points = [[point[k]+vect[k] for k in range(3)] for point in points]
            cell = MKPOL([points,[range(1,n+1)],None])
            cellBatch = Plasm.getBatches(cell)
            cellBatch[0].diffuse = YELLOW
            batches += cellBatch
            # view rotation
            rot = ROTN([ ACOS(INNERPROD([ [1,1,1],[0,0,1] ])), VECTPROD([ [1,1,1],[0,0,1] ]) ])
            batches += Plasm.getBatches(STRUCT([rot, MK([1,1,1])]))
        return batches

    def DRAW0(chain=range(1,g.getNumNode()+1)):

        m = g.getPointDim()
        d = g.getMaxDimCells()
        mapping = [dict(zip(CELLSPERLEVEL(g)(k),range(len(CELLSPERLEVEL(g)(k)))))
                   for k in range(d+1)]

        chains = [[],[],[],[]]
        [chains[g.Level(node)].append(node) for node in chain[::-1]]
        #myprint("chains",chains)
        nodepoints = [[g.getVecf(node)[i] for i in range(1,m+1)]
                      if m>2 else
                      [g.getVecf(node)[i] for i in range(1,m+1)]+[0.0]
                      for node in CELLSPERLEVEL(g)(0)]
        myprint("nodepoints",nodepoints)

        def translate(pointIds):
            #myprint("pointIds",pointIds)
            return [nodepoints[mapping[0][vert[1]]] for vert in 
                          enumerate(pointIds)]
            
        if m==2: m+=1
        if chains[0] != []: vertpoints = translate(chains[0])
        if chains[1] != []:
            edges = [DOWNCELLS(g)(edge) for edge in chains[1]]
            #myprint("edges",edges)
            edgepoints = AA(translate)(edges)
            #myprint("edgepoints",edgepoints)
        if chains[2] != []:
            facesAsEdges = [DOWNCELLS(g)(face) for face in chains[2]]
            facesAsVerts = [list(set(CAT(AA(DOWNCELLS(g))(face))))
                            for face in facesAsEdges]
            facepoints = AA(translate)(facesAsVerts)
        if d == 3:
            if chains[3] != []:
                solidsAsVerts = [UPCELLS(g)(cell) for cell in chains[3]]
                cellpoints = AA(translate)(solidsAsVerts)

        #this is the list of batches you want to display
        batches = []
        if chains[0] != []:
            myprint("vertpoints",vertpoints)
            batches = spheres(vertpoints,expl)
        if chains[1] != []:
            batches = cylinders(batches,edgepoints,expl)
        if chains[2] != []:
            batches = planecells(batches,facepoints,expl)
        if n == 3:
            if chains[3] != []: batches = cells(batches,cellpoints,expl)

        # organize the batch in a loose octree
        octree=Octree(batches)

        # create the viewer and run it
        viewer=Viewer(octree)
        viewer.Run()

    return DRAW0

###################################################################

if __name__ == "__main__":

    g = Grid([[1],[1],[1]])
    DRAW(g,[2.0,2.0,2.0])()
