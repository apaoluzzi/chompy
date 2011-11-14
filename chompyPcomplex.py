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
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
## Lesser General Public License for more details.

## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA	 02111-1307	 USA

from chompyPointSet import *

#########################################################
class PolytopalComplex(object):
#########################################################
	""" The data type to represent a bounded convex set,
		supporting a cell complex of polytopal boundary faces. """
		
	# TODO (alberto):  substitute pyplasm with scipy.spatial.qhull
	###############################################################

	## -- __init__ Method -------------------------------
	def __init__(self, points, cells=[]):
	
		def homology_maps (dictos):
			""" Compute the homology map of a (polytopal) cell complex.
			The map is given as d+1 lists
			[b_0, b_1, ..., b_d] corresponding to the boundary operators

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
			if self.dim > 0:
				skeleton = dictos[1]
				for cell in skeleton:
					homology[1].extend([(skeleton[cell], facet[0], ) \
					   for facet in pfacets(self)(eval(cell))])
				d = 1
				for skeleton in dictos[2:]:
					homology.append([])
					d += 1
					for cell in skeleton:
						for facet in pfacets(self)(eval(cell)):
							try:
								key = dictos[d-1][repr(facet)]
								homology[d].extend([(skeleton[cell], key)])
							except KeyError:
								key = dictos[d-1][repr(revert(facet))]
								homology[d].extend([(skeleton[cell], key)])
			return homology

		def convert(points, cells):
			cellsByVertices = \
				[[ points[v-1] for v in cell] for cell in cells]	 ### CHECK !!!
			#myprint("cellsByVertices",cellsByVertices)
			return [[self.vertices.dict[code(p)] for p in cell]
					for cell in cellsByVertices]

		def pack(cells):
			packedcells = sorted(AA(eval)(list(set(AA(repr)(cells)))))
			return packedcells

		#-- __init__ body -------------------------------------
		
		# TODO (alberto):  remove pyplasm and use qhull
		###############################################
		
		if cells != []:
			hpc = MKPOL([ points, [[v+1 for v in cell] for cell in cells], None])
		else:
			hpc = JOIN(AA(MK)(points))
		self.dim = DIM(hpc)
		self.rn = RN(hpc)
		
		self.cells = (self.dim + 1)*[[]]
		d = self.dim		
		
		self.vertices = PointSet(UKPOL(SKELETON(d)(hpc))[0])
		self.n = len(self.vertices.dict)
		
		for d in range(self.dim+1):
			skel = UKPOL(SKELETON(d)(hpc))
			self.cells[d] = convert(skel[0], skel[1])
			self.cells[d] = pack(AA(sorted)(self.cells[d]))						

		self.dictos = mktables(self.cells)
		self.inv_dict = [
			dict([[v,k] for k,v in self.dictos[i].items()])
				for i in range(self.dim + 1) ]
		
		self.homology = homology_maps(self.dictos)
		self.properties = {}
		self.properties["boundary"] = []

	def chain_complex (self):
		""" To compute d adjacency matrices that implement the boundary maps.
		Return a list of sparse matrices in CSR (Compressed Sparse Row) format.
		"""
		def incidence_coeff (i,j):
			face = self.inv_dict[d-1][i]
			faces = map( repr, pfacets(self)(eval(self.inv_dict[d][j])) )
			index = next((i for i in range(len(faces))
						  if faces[i] == face), None)
			if index is not None:
				if ISEVEN(index): value = 1
				else: value = -1
			return value

		matrices, d = [], 0
		for pair_set in self.homology[1:]:
			d += 1
			a = array(pair_set)
			# implicit transposition of pairs
			J = a[:,0]
			I = a[:,1]
	##		  V = array( len(a)*[1] )
			V = array([ incidence_coeff(I[h], J[h]) for h in range(len(a)) ])

			A = sparse.coo_matrix((V,(I,J)), shape=(max(I)+1, max(J)+1))
			matrices.append(A.tocsr())
		return matrices
			

	## -- __repr__ Method -------------------------------
	def __repr__(self):
		'''Return a string representation of this CellComplex in the form of a list of cells.'''
		return ('\nPolytope: \n\t.n: {0}'+ 
						'\n\t.dim: {1}'+ 
		'').format(self.n, self.dim)


	## -- __str__ Method -------------------------------

	def __str__(self):
		'''Return a string representation of this CellComplex in the form of a list of cells.'''

		if self.n < 100:
			return ('\nPolytope: \n\t.vertices: {0}'+ 
					'\n\t.cells: {1}'+ 
					'').format(self.vertices, self.cells)
		else:
			verts = self.vertices
			cells = self.cells
			return ('\nPolytope: \n\t.vertices: {0}'+ 
					'\n\t... ...'+ 
					'{1}'+ 
					'\n\t.cells: {2}'+ 
					'\n\t... ...'+ 
					'{3}'+ 
					'').format(
						verts[:20], verts[-20:],
						cells[:20], cells[-20:]
						)

	## -- copy Method -------------------------------
	def copy(self):
		""" To make a copy of a cellComplex object. """
		return copy.deepcopy(self)

	## -- view Method -------------------------------
	def view(self,d=1):
		""" To give Plasm a PolytopalComplex to visualize.

		Return 'None'.
		"""
		def plasm(d, verts, cells):
			""" To print the PLASM definition of a 'd'-skeleton.

			Returns a string, to be written on the output .psm file. The skeleton name is generated
			starting from the input file 'datafile'(.dat)
			"""
			Plasm.View(Plasm.mkpol(d, verts, cells))
			
		cells = self.cells
		if not (self.vertices == [] and cells == []):
			vertices = CAT(AA(eval)(self.vertices.ind.values()))
			dim = self.dim
			rn = self.rn
			plasm(rn, vertices, cells[d])
		else: print "cannot view an empty complex!"
		
	## -- HPC constructor -------------------------------
	def hpolc(self):
		""" To transform a PolytopalComplex into a Hierarchical
		Polyhedral Complex (HPC).

		Return an object of Hpc class.
		"""			   
		cells = self.cells
		vertices = CAT(self.vertices.points)	
		dim = self.dim
		m = self.n
		return Plasm.mkpol(dim, vertices, cells[dim])
			

	## -- export Method -------------------------------
	def ukpol (self):
		""" To the "UKPOL" of a simplicial d-complex.

		Return a tuple verts, cells[d].
		"""
		return self.vertices.points, self.cells[self.dim]

	## -- write Method -------------------------------
	def write(self, file_name):
		""" To write a '.psm' output file, with plasm syntax.

		Return 'None'.
		"""
		def plasm_def(name, d, verts, cells):
			""" To print the PLASM definition of a 'd'-skeleton.

			Returns a string, to be written on the output .psm file. The skeleton name is generated
			starting from the input file 'datafile'(.dat)
			"""
			if cells != []:
				cells = [ map(int,c) for c in cells ]
				definition = "DEF " + file_name + "_K" + str(d) + " = MKPOL:" \
						 + printer([ verts, cells, "< 1 .. " + str(len(cells)) + ">" ]) \
						 + ";\n\n"
			else: definition = ""
			return definition

		cells = self.cells
		verts = [eval(self.vertices.ind[k])
				 for k in range(len(self.vertices.dict))]

		output_file = open(file_name + ".psm", "w")
		output_file.write("% "+ str(datetime.date.today()) + \
						  " PLASM converted file " + " %\n\n" )
		output_file.write("DEF points = " + printer(verts) + ";\n\n" )
		output_file.write("DEF " + file_name + "_K0 = " +
			"MKPOL:<points, (aa:list~intsto~len):points, (aa:list~intsto~len):points>;\n\n")

		dim = self.dim
		if dim > 0:
			for d in range(1,dim+1):
				output_file.write(plasm_def( file_name, d, "points",
							[AA(C(SUM)(1))(cell) for cell in cells[d]] ))		
			for d in range(dim+1):
				output_file.write("\nVIEW: " + file_name + "_K" + str(d) + ";\n")
		output_file.close()



	## -- translate Method -------------------------------
	def translate(self, vect):
		""" To translate the object by the vect.

		Return a new object.
		"""
		obj = copy.deepcopy(self)
		obj.vertices.translate(vect)
		return obj


	## -- scale Method -------------------------------
	def scale(self, vect):
		""" To scale the object by the vect.

		Return a new object.
		"""
		obj = copy.deepcopy(self)
		obj.vertices.scale(vect)
		return obj


	## -- rotate Method -------------------------------
	def rotate (self, axis1, axis2, angle):
		""" To rotate the object by the 'angle'.
		The changed coords are 'axis1', 'axis2'

		Return a new object.
		"""
		obj = copy.deepcopy(self)
		obj.vertices.rotate(axis1, axis2, angle)
		return obj


	## -- transform Method -------------------------------
	def transform (self, matrix):
		""" To rotate the object by the 'angle'.
		The changed coords are 'axis1', 'axis2'

		Return a new object.
		"""
		obj = copy.deepcopy(self)
		obj.vertices.transform(matrix)
		return obj


	def boundary(self):
		""" To compute the boundary of a polytopal d-complex.
	
		Return the closed (polytopal) boundary (d-1)-complex .
		"""
	
		""" local copy of the input complex ------------------------------ """
		cells = self.cells
		points = self.vertices.points
		d = self.dim
		dictos = self.dictos
		boundary3 = array(self.homology[d])
	
		""" Computation of boundary faces -------------------------------- """
		V = array( len(boundary3)*[1] )
		J = boundary3[:,0]
		I = boundary3[:,1]
		A = sparse.coo_matrix((V,(I,J)), shape=(max(I)+1, max(J)+1)).tocsr()
		boundary_indices = [i for i in range(A.shape[0]) if A[i].sum() == 1]
		facetdict = dict([[index,eval(key)] for key,index in dictos[d-1].items()
			if index in boundary_indices])
	
		""" Computation of facet orientation ----------------------------- """
		orderedFacets = []
		
		for facetIndex in boundary_indices:
			facet = self.cells[d-1][facetIndex]
			facetPoints = [eval(self.vertices.ind[k]) for k in facet]
			facetCenter = centroid(self,facet)
			facetVectors = (array(facetPoints) - facetCenter).tolist()
			
			cellIndex = [c for (c,f) in boundary3 if (f==facetIndex)][0]
			cell = self.cells[d][cellIndex]
			cellCenter = centroid(self,cell)
			cellVector = VECTDIFF([cellCenter,facetCenter])
	
			transform_3d_2d = mat(facetVectors[0:2] + [cellVector])
			if linalg.det(transform_3d_2d) == 0.0:
				transform_3d_2d = mat(facetVectors[1:3] + [cellVector])
			transform_3d_2d = transform_3d_2d.I
	
			facetVectors3d = (facetVectors * transform_3d_2d).tolist()
			facetVectors2d = AA(RTAIL)(facetVectors3d)
			ordering = sorted([[atan2(vect),k]
							   for k,vect in enumerate(facetVectors2d)])
			order = [ind for angle,ind in ordering]
			orderedFacet = [facet[k] for k in order]
			
			vol = linalg.det(transform_3d_2d)
			if vol < 0: orderedFacet = REVERSE(orderedFacet)
			orderedFacets += [orderedFacet]
	
		orderedFacetsByPoints = [[self.vertices.ind[k] for k in facet]
								 for facet in orderedFacets]
			
	
		""" Remove unused cells from input complex (w remapping) --------- """
		
		out = PolytopalComplex(points, facetdict.values())
		orderedFacets = [[out.vertices.dict[pointKey] for pointKey in face]
						 for face in orderedFacetsByPoints]
		out.cells[-1] = orderedFacets
		return out
		
		
	def project(self):
		"""
		Projection of the d-complex boundary on the subspace of the first n-1 coordinates.

		Return a polytopal (d-1)-complex.
		"""
		obj = self.boundary()
		
		# remove the back-faces from the boundary ----------------------------
		outvert = (array(obj.vertices.max()) + ((obj.rn-1)*[0]+[1000])).tolist()  # to make better

		def frontFace(verts):
			verts = [array(v) - verts[0] for v in (verts[1:3] + [outvert])]
			transformation = mat(verts)
			transformation = transformation.I
			if linalg.det(transformation) > 0.0: return True
			else: return False
					
		frontFaces = [face for face in obj.cells[-1] 
						if frontFace([obj.vertices.points[k] for k in face])]
		
		points = obj.vertices.project(1).points
		ptk,verts = remap(points)
		
		def pred(face): 
			"""
			Predicate to test if there are no duplicates in a sequence.
			Return `True` if there are no duplicate elements in the input parameter.
			"""
			return COMP([C(EQ)(len(face)),len,set])(face)
			
		cells = [AA(ptk)(face) for face in frontFaces if pred(AA(ptk)(face))]
		return PolytopalComplex(verts.points,cells)
		



def complexProd(args):
	pol1,pol2 = args
	p1, p2 = pol1.vertices.points, pol2.vertices.points
	m1, m2 = len(p1), len(p2)
	points = [ p1[i]+p2[j] for i in range(m1)
			  for j in range(m2) ]

	def ind(args):
		i,j = args
		return i*m2 + j

	c1, c2 = pol1.cells[-1], pol2.cells[-1]
	n1, n2 = len(c1), len(c2)
	cells = [ AA(ind)(CART([c1[i],c2[j]])) for i in range(n1)
			  for j in range(n2) ]

	return PolytopalComplex(points,cells)

