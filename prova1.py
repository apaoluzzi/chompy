from chompy.chompyScomplex import *


def offset(point,expl=[1,1,1]):
    scaledpoint = [point[k]*expl[k] for k in range(3)]
    vect = VECTDIFF([scaledpoint,point])
    return vect
    

def spheres(batches,points,expl=[1,1,1]):
    sx = 0.05
    unitSphere = Batch.openObj("sphere18x27.obj")[0]
    points = CAT(points)
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

#/////////////////////////////////////////////////////////

def draw (c,chains=4*[[]],expl=[1,1,1]):

    embeddingDim = 3 - c.rn     # possible curation of embedding dimension
    if embeddingDim != 0: c.vertices.embed(embeddingDim)

    if len(chains)<(c.dim+1):   # possible curation of input chains
        chains = chains + ((c.dim+1)-len(chains))*[[]]

    if chains == [[],[],[],[]]: chains = c.cells # draw the whole complex
    else:                                        # draw only the input chains
        chains = [[c.cells[h][k] for k in chains[h]] for h in range(len(chains))]
    primitives = [spheres, cylinders, planecells, cells]
    batches = []
    
    def cellverts(c,cells):
        return [[eval(c.vertices.ind[v]) for v in cell] for cell in cells]

    def addBatches(objType,batches,items,expl):
        return objType(batches,items,expl)

    for k in range(c.dim + 1):
        if chains[k] != []:
            items = cellverts(c,chains[k])
            batches = addBatches(primitives[k],batches,items,expl)

    octree=Octree(batches)
    viewer=Viewer(octree)
    viewer.Run()

#/////////////////////////////////////////////////////////

c = simplexGrid([2*[1.], 2*[1.], 1*[1.]])
c.view()
draw(c)
draw(c,expl=[1,1,2],chains=[[],range(len(c.cells[1]))])


