from chompy import *

##    DEF remove (n::IsIntPos; seq::IsSeq) =
##        (CONS ∼ AA:SEL ∼ CAT):< 1 .. n - 1, n + 1 .. LEN:seq >:seq;
##
##    DEF permutations (seq::IsSeq) =
##        COMP:(AA:CAT AL N:n:(CAT ∼ AA:permute)):<< <>, seq>>
##        WHERE
##            n = LEN:seq,
##            extract = AA: remove ∼ DISTR ∼ [INTSTO ∼ LEN, ID],
##            permute = TRANS ∼ [ AA:AR ∼ DISTL, extract ∼ S2 ]
##        END;
##
##    permutations:<1,2,3> 
##


def remove(args):
    n,seq = args
    return seq[:n] + seq[n+1:]

remove([5, range(10)])


def permutations (SEQ):
    if len(SEQ) <= 1: return [SEQ]
    ret = []
    for i in range(len(SEQ)):
        element = SEQ[i]
        rest    = permutations(SEQ[0:i] + SEQ[i+1:])
        for r in rest: ret+=[[element] + r]
    return ret

print permutations([1,2,3,4])

def triangle_strip(points):
	"""
	Primitive constructor of a strip of triangles.
	
	"""
	ptk,verts = remap(points)
	cells = [AA(ptk)([k,k+1,k+2]) for k in range(len(points)-2)]	
	return SimplicialComplex(points,cells)
    

##    DEF permutahedron (d::IsIntPos) =
##        ( project:1 ∼ rotations ∼ translation ):object
##        WHERE
##            object = (MKPOL ∼ [ID, [INTSTO ∼ LEN], K:<<1>>]):vertices,
##            vertices = permutations:(1 .. (d+1)),
##            center = Meanpoint: vertices,
##            translation = T:(1..(d+1)):(AA:-: center),
##            rotations = COMP:(((CONS ∼ AA:R):(1..d DISTR (d+1))):(PI/4))
##        END;


def permutahedron(d): 
    vertices = permutations(range(1,d+2))
    myprint("vertices",vertices)
    obj = PolytopalComplex(vertices,[range(len(vertices))])
    center = centroid(obj,range(len(vertices)))
    myprint("center",center)
    translation = t(range(d+1))([-x for x in center])
    rotations = COMP(CONS(AA(r)(DISTR([range(1,d),d+1])))(pi/4))
    return rotations(translation(obj))#.project()


permutahedron(2).view()
