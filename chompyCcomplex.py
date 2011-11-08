## This file provides a library for polytopal and simplicial complexes.
## Includes input, output, skeleton and boundary evaluation, linear
## extrusion, boundary and coboundary operators, and linear combination
## of chains.
## Author: Alberto Paoluzzi (paoluzzi@dia.uniroma3.it)
## Copyright (C) 2009,2010 Dipartimento Informatica e Automazione,
## Università Roma Tre, Rome, Italy.

## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.

## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.

## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

from chompyScomplex import *


def incr(listOfLists): return (AA(AA(C(SUM)(1))))(listOfLists)

def findsubsets(S,m):
    return set(itertools.combinations(S, m))

def offset(shape):
    """ Translation from multidimensional array index to linear offset.
        Formula: offset(D) = D0xS1x...xSn + D1xS2x...xSn + D(n-1)xSn + Dn
    """
    n = len(shape)
    def offset0(index):
        return SUM([index[d] * PROD(shape[d+1:n])
                    for d in range(n-1)]) + index[n-1]
    return offset0

def index(shape):
    """ index(shape)(offset) is the inverse of offset(shape)(index).
        To convert an integer offset to a multidimensional array index.
    """
    n = len(shape)
    def index0(offset):
        ind = []
        for k in range(n):
            ind += [offset % shape[n-k-1]]
            offset = offset / shape[n-k-1]
        return REVERSE(ind)
    return index0


if __name__ == "__main__" and __debug__:
    assert (offset([3,3,4])([0,0,0]) == 0) 
    assert (index([3,3,4])(0) == [0,0,0])
    assert (offset([3,3,4])([0,2,2]) == 10)
    assert (index([3,3,4])(10) == [0,2,2])
    assert (offset([3,3,4])([1,1,2]) == 18)
    assert (index([3,3,4])(18) == [1,1,2])
    assert (offset([3,3,4])([2,2,2]) == 34)
    assert (index([3,3,4])(34) == [2,2,2])

def pgrid(sides):
    """ To generate a d-dimensional (hyper-)cuboid simplexGrid.
        shape = shape of the node array.
    """

    def checkBounds (cell):
        return AND([AND(array(point) < array(shape)) for point in cell])

    shape = [s+1 for s in AA(len)(sides)]
    nverts = PROD(shape)
    verts = CART(AA(progressive_sum)(sides))
    dim = len(sides)
    cells = AA(LIST)(range(dim + 1))
    cells[0] = [[k] for k in range(nverts)]

    for d in range(1,dim+1):
        axes = AA(list)(list(findsubsets(range(dim),d)))
        for k in range(nverts):
            candidates = CART([ [index(shape)(k)], axes ])
            for pair in candidates:
                baseIndex = pair[0]
                newIndices = [baseIndex]
                for i in pair[1]:
                    displacement = dim * [0]
                    displacement[i] = 1
                    newIndices += [displacement]
                cell0 = [newIndices[0]]
                cell1 = []
                for j in range(1,len(newIndices)):
                    cell0 += [VECTSUM([c,newIndices[j]]) for c in cell0]
                    if checkBounds(cell0):
                        cell2 = AA(offset(shape))(cell0)
                        if len(cell2) == 2**d: cells[d] += [cell2]
        cells[d] = AA(sorted)(TAIL(cells[d]))
    return verts,cells

def cfacets (ccomplex):
    def cfacets0 (cell):
        """ To compute the (non oriented) boundary of a d-cuboid.

        Uses a common method of facet extraction:
        -  convert from offset to multiindex representation;
        -  min-max computation for each one of the d axes.
        Both the d-cuboid and its 2*d (d-1)-facets are represented as
        an ordered list of vertices.

        Return a list of facets, i.e. of (d-1)-dimensional cuboids.
        """
        def mask(cell): return AA(EQ)(TRANS(cell))

        def ISNOT(value): return not value

        def multiFilter(cell):
                boolmask = mask(cell)
                return [[p[k] for k in range(len(p))
                         if ISNOT(boolmask[k])] for p in cell]

        def inverseMFilter(oldCell):
            def inverseMFilter0(newCell):
                boolmask = mask(oldCell)
                h,out = 0,[]
                for k in range(len(boolmask)):
                    if boolmask[k] == TRUE:
                        out += [oldCell[0][k]]
                    else:
                        out += [newCell[h]]
                        h += 1
                return out
            return inverseMFilter0

        dim = int(LN(len(cell))/LN(2))
        if dim == 1: return AA(LIST)(cell)

        cell0 = [index(ccomplex.shape)(k) for k in cell]
        cell1 = multiFilter(cell0)
        
        faces = []
        
        for d in range(dim):
            min_index = array(TRANS(cell1))[d].min()
            max_index = array(TRANS(cell1))[d].max()
            faces += [[v for v in cell1 if v[d] == min_index]]
            faces += [[v for v in cell1 if v[d] == max_index]]
            
        faces0 = AA(AA(inverseMFilter(cell0)))(faces)
        faces1 = AA(AA(offset(ccomplex.shape)))(faces0)
        return faces1
    return cfacets0


#########################################################
class CuboidalComplex(PolytopalComplex):
#########################################################
    """ The data type to represent a grid of (hyper)cuboids.
        Efficient generation of every subset of d-dimensional faces.
    """    

    def mktables (cells):
        """ To build the cell database of a cell complex.

        Return a list of d+1 dictionaries, where d is the dimension of the complex.
        Key is the 'repr' of the cells; values are integers in range(0,k[d]) where
        k[d] is the cardinality of the d-skeleton.
        """
        # resetting dictionaries
        # (so that each call to mktables results in a new set of dictionaries)
        dictos = []
        
        
        def revert(cell):
            """ Change the sign of a permutation by exchanging the first and the last elements
            """
            if len(cell) > 1:
                cell = [cell[-1]] + cell[1:-1] + [cell[0]]
            return cell

        d = 0
        for skeleton in cells:
            dictos.append({})
            [operator.setitem(dictos[d], repr(cell), skeleton.index(cell)) \
                for cell in skeleton \
                    if not (dictos[d].has_key(repr(cell)) or \
                             dictos[d].has_key(repr(revert(cell))))]
            d += 1
        return dictos




    ## -- __init__ Method -------------------------------
    def __init__(self, sides):


        def homology_maps (dictos):
            """ Compute the homology map of a (simplicial) cell complex.
            The map is given as d+1 lists (4 lists in this prototype implementation)
            [b_0, b_1, b_2, b_3] corresponding to the boundary operators

            b_d: K_d -> K_(d-1).  Notice that b_0 == [] by definition.
            Each map is given as an (ordered) list of pairs.
             
            Return a list of lists of pairs (tuples) of positive integer indices.
            The first element is the index of a d-cell; the second element is the
            index of an incident (d-1)-facet.
            """

            def revert(cell):
                """ Change the sign of a permutation by exchanging the first and the last elements
                """
                if len(cell) > 1:
                    cell = [cell[-1]] + cell[1:-1] + [cell[0]]
                return cell

            homology = [[],[]]
            skeleton = dictos[1]
            for cell in skeleton:
                homology[1].extend([(skeleton[cell], facet[0], ) \
                   for facet in cfacets(self)(eval(cell))])
            d = 1
            for skeleton in dictos[2:]:
                homology.append([])
                d += 1
                for cell in skeleton:
                    for facet in cfacets(self)(eval(cell)):
                        try:
                            key = dictos[d-1][repr(facet)]
                            homology[d].extend([(skeleton[cell], key)])
                        except KeyError:
                            key = dictos[d-1][repr(revert(facet))]
                            homology[d].extend([(skeleton[cell], key)])
            return homology
        

        def convert(points, cells):
            cellsByVertices = \
                [[ points[v-1] for v in cell] for cell in cells]
            return [[self.vertices.dict[code(p)] for p in cell]
                    for cell in cellsByVertices]

        def pack(cells):
            return sorted(AA(eval)(list(set(AA(repr)(cells)))))



        verts, cells = pgrid(sides)

        self.shape = [s+1 for s in AA(len)(sides)]
        dim = len(sides)
        self.dim = int(LN(len(cells[-1][0]))/LN(2))
        assert dim == self.dim 
        self.rn = len(verts[0])

        self.vertices = PointSet(verts)
        self.cells = AA(LIST)(range(self.dim + 1))
        self.n = len(self.vertices.dict)

        for d in range(self.dim+1):
            self.cells[d] = cells[d]
            self.cells[d] = pack(AA(sorted)(self.cells[d]))

        self.dictos = mktables(self.cells)
        self.inv_dict = [
            dict([[v,k] for k,v in self.dictos[i].items()])
                for i in range(self.dim + 1) ]

        self.homology = homology_maps(self.dictos)
        #myprint("self.homology",self.homology)
                        

    ## -- Skeleton extraction Method ---------------------
    def skeleton(self,d):
        """ returns the d-skeleton of a cuboidal complex """
        return self.cells[d]

    def chain_complex (self):
        """ To compute d adjacency matrices that implement the boundary maps.

        Return a list of sparse matrices in CSR (Compressed Sparse Row) format.
        """
        def incidence_coeff (i,j):
            #myprint("(i,j)",(i,j))
            face = self.inv_dict[d-1][i]
            #myprint("face",face)
            faces = map( repr, cfacets(self)(eval(self.inv_dict[d][j])) )
            #myprint("faces",faces)
            index = next((i for i in range(len(faces))
                          if faces[i] == face), None)
            if index is not None:
                if ISEVEN(index): value = 1
                else: value = -1
            #myprint("value",value)
            return value

        matrices, d = [], 0
        for pair_set in self.homology[1:]:
            d += 1
            a = array(pair_set)
            # implicit transposition of pairs
            J = a[:,0]
            I = a[:,1]
    ##        V = array( len(a)*[1] )
            V = array([ incidence_coeff(I[h], J[h]) for h in range(len(a)) ])

            A = sparse.coo_matrix((V,(I,J)), shape=(max(I)+1, max(J)+1))
            matrices.append(A.tocsr())
        return matrices

    


#     ## -- boundary Method -------------------------------
#     def boundary (self):
#         """ To compute the boundary of a cuboidal d-complex.
# 
#         Return the closed (polytopal) boundary (d-1)-complex .
#         """
#         obj = self.copy()
#         cells = obj.cells
#         vertices = obj.vertices.points
#         d = obj.dim
#         dictos = obj.dictos
#         h = obj.homology
#         a = array(h[d])
# 
#         V = array( len(a)*[1] )
#         J = a[:,0]
#         I = a[:,1]
# 
#         A = sparse.coo_matrix((V,(I,J)), shape=(max(I)+1, max(J)+1)).tocsr()
# 
#         boundary_indices = [i for i in range(A.shape[0]) if A[i].sum() == 1]
#         
#         
#         inv_dict = dict([[v,k] for k,v in dictos[d-1].items()])
#         faces = [eval(inv_dict[k]) for k in boundary_indices]
#         
#         return PolytopalComplex(vertices, faces)
# 
# 
	## -- boundary Method -------------------------------
	def boundary (self):
		""" To compute the boundary of a simplicial d-complex.

		Return the closed (d-1)-complex triangulating the boundary of the input complex.
		"""
		
		myprint(">>>>>>>>>>"," eccomi in Ccomplex!!")
		
		obj = copy.deepcopy(self)
		cells = obj.cells
		if obj.dim == -1:  return SimplicialComplex([], [])
		vertices = obj.vertices.points
		d = obj.dim
		dictos = obj.dictos
		h = obj.homology
		a = array(h[d])

		V = array( len(a)*[1] )
		J = a[:,0]
		I = a[:,1]

		A = sparse.coo_matrix((V,(I,J)), shape=(max(I)+1, max(J)+1)).tocsr()

		# make boundary orientation coherent  --------------------------

		def simplex(cell):
			point = obj.vertices.ind
			return [eval(point[k])+[1.0] for k in cell]
		
		def volume(cell):
			return linalg.det(mat(simplex(cell)))

		def orientation():
			if d == obj.rn:   # solid complex
				out = [volume(cell) for cell in cells[-1]]
			else:				# embedded complex
				out = [linalg.det(linalg.qr(mat(simplex(cell)))[1][:,:-1])
					   for cell in cells[-1]]  # DEBUG (choose minor with det(minor != 0	))
			return out		

		
		boundary_indices = [i for i in range(A.shape[0]) if A[i].sum() == 1 ] 
			
		boundary_signs = orientation()   

		boundary_pairs = [(i,j) for (i,j) in a 
			if (j in boundary_indices)] 


		facetsdict = dict([[v,k] for k,v in dictos[d-1].items()])
		cellsdict = dict([[v,k] for k,v in dictos[d].items()])
						
		def invertOrientation(facet):
			facet[0],facet[-1] = facet[-1],facet[0]
			return facet


		facets = [eval(facetsdict[k]) for k in boundary_indices]
 		facets = [invertOrientation(eval(facetsdict[facet]))  
 					if boundary_signs[face]<0 else eval(facetsdict[facet])
 					for face,facet in boundary_pairs]
				
		# remapping section -----------------------------------------------
		
		if facets != []:
			oldinds = list(set(CAT(facets)))
			newverts = PointSet([eval(obj.vertices.ind[k]) for k in oldinds])
			newfacets = [[ newverts.dict[obj.vertices.ind[k]] 
							for k in facet] for facet in facets]
			return SimplicialComplex(newverts.points, newfacets)
		else: return SimplicialComplex([], [])


