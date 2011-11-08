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
				[[ points[v-1] for v in cell] for cell in cells]   ### CHECK !!!
			
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
		""" To translate the SimplicialComplex by the vect as side effect.

		Return self.
		"""
		self.vertices.translate(vect)
		return self


	## -- scale Method -------------------------------
	def scale(self, vect):
		""" To scale the SimplicialComplex by the vect as side effect.

		Return self.
		"""
		self.vertices.scale(vect)
		return self


	## -- rotate Method -------------------------------
	def rotate (self, axis1, axis2, angle):
		""" To rotate the SimplicialComplex by the 'angle' as side effect.
		The changed coords are 'axis1', 'axis2'

		Return self.
		"""
		self.vertices.rotate(axis1, axis2, angle)
		return self


	## -- transform Method -------------------------------
	def transform (self, matrix):
		""" To rotate the SimplicialComplex by the 'angle' as side effect.
		The changed coords are 'axis1', 'axis2'

		Return self.
		"""
		self.vertices.transform(matrix)
		return self


	## -- boundary Method -------------------------------
	def boundary (self):
		""" To compute the boundary of a polytopal d-complex.

		Return the closed (polytopal) boundary (d-1)-complex .
		"""
		myprint(">>>>>>>>>>"," eccomi in Pcomplex!!")

		obj = self.copy()
		cells = obj.cells
		vertices = obj.vertices.points
		d = obj.dim
		dictos = obj.dictos
		h = obj.homology
		a = array(h[d])

		V = array( len(a)*[1] )
		J = a[:,0]
		I = a[:,1]

		A = sparse.coo_matrix((V,(I,J)), shape=(max(I)+1, max(J)+1)).tocsr()

		boundary_indices = [i for i in range(A.shape[0]) if A[i].sum() == 1]

		# make boundary orientation coherent  --------------------------
		boundary_pairs = [(i,j) for (i,j) in a if (j in boundary_indices)] 
		facetsdict = dict([[v,k] for k,v in dictos[d-1].items()])
		cellsdict = dict([[v,k] for k,v in dictos[d].items()])
		pairs = [(cellsdict[c], facetsdict[f]) for c,f in boundary_pairs]
		
		def getPointUnderFace(pair):
			cell,face = AA(eval)(pair)
			vert = list(set(cell).difference(face))[0]
			belowPoint = eval(obj.vertices.ind[vert])
			faceBarycenter = centroid(obj,face)
			facePoints = [eval(obj.vertices.ind[v]) for v in face]
			return [belowPoint,faceBarycenter,facePoints]

		triples = [getPointUnderFace(pair) for pair in pairs]
		belowPoints,faceBarycenters,facePoints = TRANS(triples)
		faceVectors = [[list(point - triple[1]) for point in AA(array)(triple[2])] 
					for triple in triples]
		transforms_3d_2d = AA(mat)(
			[vect3d[0:2] + [list(array(belowPoints[k])-faceBarycenters[k])] 
					for k,vect3d in enumerate(faceVectors)])
		transforms_3d_2d = [transform.I for transform in transforms_3d_2d]
		faceVectors3d = [(faceVectors[k] * transform).tolist()
						for k,transform in enumerate(transforms_3d_2d)]
		faceVectors2d = [AA(RTAIL)(vect3d) for vect3d in faceVectors3d]
		
		def simplex(cell):
			out = [point+[1.0] for point in cell]
			return out
		
		def volume(cell):
			return linalg.det(mat(simplex(cell)))
			
		faces = []
		for k,face in enumerate(facePoints):
			vol = volume(face[:-1]+[belowPoints[k]])
			verts = eval(pairs[k][1])
			if vol>0:  faces += [[ verts[0], verts[1], verts[3], verts[2] ]]
			else:  faces += [[ verts[2], verts[3], verts[1], verts[0] ]]

#		inv_dict = dict([[v,k] for k,v in dictos[d-1].items()])
#		faces = [eval(inv_dict[k]) for k in boundary_indices]
		return PolytopalComplex(vertices, faces)


#	## -- boundary Method -------------------------------
#	def boundary (self):
#		""" To compute the boundary of a simplicial d-complex.
# 
#		Return the closed (d-1)-complex triangulating the boundary of the input complex.
#		"""
#		obj = self.copy()
#		cells = obj.cells
#		vertices = obj.vertices.points
#		d = obj.dim
#		dictos = obj.dictos
#		h = obj.homology
#		a = array(h[d])
# 
#		V = array( len(a)*[1] )
#		J = a[:,0]
#		I = a[:,1]
# 
#		A = sparse.coo_matrix((V,(I,J)), shape=(max(I)+1, max(J)+1)).tocsr()
# 
#		# make boundary orientation coherent  --------------------------
# 
#		def simplex(cell):
#			point = self.vertices.ind
#			return [eval(point[k])+[1.0] for k in cell]
#		
#		def volume(cell):
#			return linalg.det(mat(simplex(cell)))
# 
#		def orientation():
#			if d == self.rn:   # solid complex
#				out = [volume(cell) for cell in cells[-1]]
#			else:				# embedded complex
#				out = [linalg.det(linalg.qr(mat(simplex(cell)))[1][:,:-1])
#					   for cell in cells[-1]]  # DEBUG (choose minor with det(minor != 0	))
#			return out		
# 
#		
#		boundary_indices = [i for i in range(A.shape[0]) if A[i].sum() == 1 ] 
#			
#		boundary_signs = orientation()	 
# 
#		boundary_pairs = [(i,j) for (i,j) in a 
#			if (j in boundary_indices)] 
# 
# 
#		facetsdict = dict([[v,k] for k,v in dictos[d-1].items()])
#		cellsdict = dict([[v,k] for k,v in dictos[d].items()])
#						
#		def invertOrientation(facet):
#			facet[0],facet[-1] = facet[-1],facet[0]
#			return facet
# 
# 
#		facets = [eval(facetsdict[k]) for k in boundary_indices]
#		facets = [invertOrientation(eval(facetsdict[facet]))  
#					if boundary_signs[face]<0 else eval(facetsdict[facet])
#					for face,facet in boundary_pairs]
#				
#		# remapping section -----------------------------------------------
#		
#		if facets != []:
#			oldinds = list(set(CAT(facets)))
#			newverts = PointSet([eval(obj.vertices.ind[k]) for k in oldinds])
#			newfacets = [[ newverts.dict[obj.vertices.ind[k]] 
#							for k in facet] for facet in facets]
#			return SimplicialComplex(newverts.points, newfacets)
#		else: return SimplicialComplex([], [])


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

