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
##    DEF permutahedron (d::IsIntPos) =
##        ( project:1 ∼ rotations ∼ translation ):object
##        WHERE
##            object = (MKPOL ∼ [ID, [INTSTO ∼ LEN], K:<<1>>]):vertices,
##            vertices = permutations:(1 .. (d+1)),
##            center = Meanpoint: vertices,
##            translation = T:(1..(d+1)):(AA:-: center),
##            rotations = COMP:(((CONS ∼ AA:R):(1..d DISTR (d+1))):(PI/4))
##        END;


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


##a = cuboid([1,1,1]).vertices
##print a
##a.project(1)
##print a.project(1)

    


project(cuboid([1,1,1])).view()

a = project(cuboid([1,1,1]))
draw(a)

simplex(3).view()
a = project(simplex(3))
draw(a)



a = project(simplexGrid([4*[1],3*[1],2*[1]]))
a.view()
draw(a)





    



def permutahedron(d): pass
##    obj =
##    vertices =
##    center =
##    translation =
##    rotations = 
##
##    return COMP([project(1),rotations,translation])(obj)
