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

from chompyPcomplex import *
from copy import deepcopy



def remap(points):
	verts = PointSet(points)
	def remap0(k): return verts.dict[code(points[k])]
	return remap0,verts



def _init_CellComplex (verts, d_simplices):
	
		vertices = PointSet(verts)

		def identify(d_simplices):
			new_simplices = []
			for simplex in d_simplices:
				new_simplex = []
				for v in simplex:
					if type(v) == list: v = v[0]
					point = map(lambda x: round_or_zero(x),
								vertices.points[v])
					index = vertices.dict[code(point)]
					new_simplex = new_simplex + [index]
				new_simplices += [new_simplex]
			return new_simplices
		
		d_simplices = identify(d_simplices)
		return (map(eval, vertices.ind.values()), d_simplices)

#########################################################
class SimplicialComplex(PolytopalComplex):
#########################################################
	""" The data type to represent a simplicial complex.

	The SimplicialComplex is a dimension-independent ADT (Abstract Data Type),
	characterized by an `intrinsic` and an `embedding dimension`, a set of `vertices`,
	and sets of `cells` (of various dimensions), given as lists of vertex indices.

	The constructor returns an **object** of the class.
	"""
	## -- __init__ Method -------------------------------
	def __init__(self, vertices, d_simplices):
		""" A new simplicial complex with given 0-cells and d-cells. """

		if vertices != []:
			vertices, d_simplices = _init_CellComplex (vertices, d_simplices)
			self.vertices = PointSet(vertices)
			self.rn = len(vertices[0])
			
			self.dim = len(d_simplices[0]) - 1
			self.cells = (self.dim + 1)*[[]]
			self.cells[0] = map(lambda x: [x], self.vertices.ind.keys())
			self.cells[self.dim] = d_simplices
			
			self.dictos = mktables(cell_complex(self.cells))
			self.inv_dict = [
				dict([[v,k] for k,v in self.dictos[i].items()])
					for i in range(self.dim + 1) ]
			self.homology = homology_maps(self.dictos)
			self.cells[-1] = map(eval, self.inv_dict[-1].values())
			
		else:
			self.vertices = []
			self.rn = -1
			self.dim = -1
			self.cells = (self.dim + 1)*[[]]
			self.inv_dict = []	
	
	## -- __repr__ Method -------------------------------
	def __repr__(self):
		'''Return a string representation of this SimplicialComplex in the form of a list of cells.'''
		return 'SimplicialComplex: %s' % self.inv_dict

	## -- __str__ Method -------------------------------
	def __str__(self):
		'''Return a string representation of this SimplicialComplex in the form of a list of cells.'''
		return ('SimplicialComplex: \n\t.n: {0}'+ 
						'\n\t.dim: {1}'+ 
						'\n\t.vertices: {2}'+ 
						'\n\t.cells: {3}'+ 
						'\n\t.dictos: {4}'+ 
		'').format(self.rn, self.dim, self.vertices, self.cells, self.dictos)



	## -- boundary Method -------------------------------
	def boundary (self):
		""" To compute the boundary of a simplicial d-complex.

		Return the closed (d-1)-complex triangulating the boundary of the input complex.
		"""
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

	## -- extrude Method -------------------------------
	def extrude (self, hlist):
		""" To complute the multiple linear extrusion of a d-complex.
		Map R^d -> R^(d+1), according to: Ferrucci & Paoluzzi, CAD 1991.
		'cells' is a simplicial complex, given as a list of lists;
		'hlist' is a list of heights in the added dimension.

		Return a (d+1)-dimensional simplicial complex.
		Only the 0- and (d+1)- skeletons are computed.
		In order to get the remaining skeletons, the 'SimplicialComplex' function
		is used.
		"""

		cells = self.cells
		dim = self.dim
		verts = self.vertices.points
		lastcoords = progressive_sum(AA(abs)(hlist))
		
		if dim == 0:
			cells = [[],[]]
			vertices = AA(LIST)(lastcoords)
			cells[1] = [[i,i+1] for i in range(1,len(hlist)+1)]
		else:
			simplexes = cells[dim]
			nverts = len(verts)
			nsteps = len(lastcoords)
			sliced_vertices = nsteps*[verts]
			def coords_distribute(x):
				return COMP([ CAT, AA(COMP([ AA(AR), DISTR ])) ])(x)
			vertices = coords_distribute(TRANS([nsteps*[verts],lastcoords]))

			extruded_simplices = []
			for cell in simplexes:
				vertPtrs = cell + map(lambda x: x+nverts, cell)  
				extruded_simplices += subcomplex(dim+2,vertPtrs)
		
			final_simplices = []
			for i in range(nsteps-1):
				if hlist[i] > 0:
					simplex_layer = shift(nverts*i,extruded_simplices)
					final_simplices += simplex_layer
			cells.append(final_simplices)
			
		return SimplicialComplex(vertices, cells[-1])


	## -- Embed Method -------------------------------
	def embed(self,dimensionIncrement):
		pol = deepcopy(self)
		points = [point + dimensionIncrement*[0] for point in pol.vertices.points]
		return SimplicialComplex(points,pol.cells[-1])

	## -- II Method -------------------------------
	def II(self,alpha,beta,gamma):
		""" To compute surface integrals of a three-variate monomial on the complex boundary.

		Return a number.
		"""
		if self.rn != 3:
			print "Integration error: 'n' != 3"
		elif self.dim == 3:
			self = self.boundary()
		elif self.dim == 2:
			pass
		else:
			print "warning: obj.dim != 3"

		surface = complex2surface(self)
		w = 0.0
		for triangle in surface:
			w += T3(triangle,alpha,beta,gamma)
		return w


	## -- III Method -------------------------------
	def III(self,alpha,beta,gamma):
		""" To compute volume integrals of a three-variate monomial on the complex boundary.

		Return a number.
		"""
		if self.rn != 3:
			print "Integration error: 'n' != 3"
		elif self.dim == 3:
			self = self.boundary()
		elif self.dim == 2:
			pass
		else:
			print "warning: obj.dim != 3"

		surface = complex2surface(self)
		
		def magnitude(vect):
			return math.sqrt(sum(x*x for x in vect))

		w = 0.0
		for triangle in surface:
			t = mat(triangle)
			a = (t[1] - t[0]).tolist()[0]
			b = (t[2] - t[0]).tolist()[0]
			c = [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]
			w += (c[0] / magnitude(c)) * T3(triangle,alpha+1,beta,gamma)
		return w/(alpha + 1)


	## -- multi-resolution Method -------------------------------
	def multires (self, C2):
		""" To compute a multiresolution cellcomplex generated by embedding
		the C2 cellcomplex within self.
		C2 must be a SUBCOMPLEX of self.

		Return a new SimplicialComplex.
		"""
		
		def pointset (simplices):
			""" To compute the list of (indices of) vertices of a simplex.

			Return a list of lists of integers.
			"""
			points = []
			for simplex in simplices:
				points += simplex
			return set(points)

		def swap(args):
			""" To swap the two (first) elements of a list.

			Return a list of two elements.
			"""
			return [args[1], args[0]]

		C1 = self.copy()
		
		# splitta il complesso interno
		S2 = C2.split()
		C2 = SimplicialComplex(S2.vertices.points, C2.cells[2])

		# calcola il bordo interno iniziale e l'insieme dei suoi punti
		C23 = C2.boundary()
		C23_points = pointset(C23.cells[-1])
		
		# calcola il bordo interno splittato e l'insieme dei suoi punti
		S23 = C23.split()
		C2_points = pointset(C2.cells[-1])
	  
		# calcola le celle (da splittare) adiacenti al bordo
		C3_cells = [s for s in C1.cells[-1]
					if len(set(s).difference(C2_points)) == 1 ]

		if C3_cells != []:

			C3 = SimplicialComplex(S23.vertices.points, C3_cells)
			C3_sets = map(set, C3_cells)

			# splitta il complesso adiacente al bordo interno
			k = -1
			S3 = []
			for cell in C23.cells[-1]:
				v = [c.difference(cell) for c in C3_sets
					 if len(c.difference(cell)) == 1 ]
				k += 1
				if v != []:
					v = v[0]
					S3 += [ S23.cells[-1][2*k] + list(v), S23.cells[-1][2*k+1] + list(v) ]

			C2 = C2.add(C3)
			C1 = C1.subtract(C2)

			return SimplicialComplex (S23.vertices.points, C1.cells[-1] + S2.cells[-1] + S3)		
		else:
			return SimplicialComplex (S23.vertices.points, S2.cells[-1])





	## -- trivial split Method -------------------------------
	def trivial_split (self):

		""" To compute the baricentric decomposition of a d-dimensional
		obj SimplicialComplex.

		Return a simplicial SimplicialComplex object with k_d * (d+1)! d-cells, where
		k_d is the cardinality of the d-skeleton of obj. 
		"""
		# initialize
		obj = self.copy()
		verts = obj.vertices.points
		dim = obj.dim
		n = obj.rn
		cells = obj.cells
		dictos = obj.dictos
		inv_dict = obj.inv_dict


		def facenode (face):
			""" To compute the index of the node associated to a face.

			Return an integer number.
			"""
			point = centroid(self,face)
			nodekey = inv_dict[0][ str(point) ]
			return nodekey

		def split1_(face):
			""" To barycentrically split a face and its boundary.

			Return a list of nested pairs.
			"""
			node = facenode(face)
			return [ [node] + f for f in facets(face)]

		def sign (face):
			""" To compute the sign of a d-simplex in a d-complex.

			Return either 1 or -1 according with the permutation class of the face.
			"""
			vertices = [verts[v] + [1.0] for v in face]
			return cmp(linalg.det(vertices),0)
			

		# initialize
		obj = self.copy()
		verts = obj.vertices.points
		dim = obj.dim
		n = obj.rn
		cells = obj.cells
		dictos = obj.dictos
		inv_dict = obj.inv_dict
		
		# add new vertices
		base=len(verts)-1
		k = 0
		for cell in cells[-1]:
			p = centroid(self,cell)
			verts.append(p)
			k += 1
			operator.setitem(inv_dict[0], str(p), base+k)
			operator.setitem(dictos[0], str([base+k]), base+k )

		# build the new d-cells
		ncells = dim+1 #number of splitted cells per d-cell
		oldcells = list(cells[-1])
		cells[-1] = []
		for cell in oldcells:
			cells[-1] += split1_(cell)
						
		 # coherent orientation of the d-cells
		for i in range(len(cells[-1])):
			cell = cells[-1][i]
			if dim == n == len(cell) - 1:
				if sign(cell) < 0:
					cells[-1][i] = [cell[-1]] + cell[1:-1] + [cell[0]]			  

		return SimplicialComplex(verts, cells[-1])


   ## -- split Method -------------------------------
	def split (self):
		""" To compute the baricentric decomposition of a d-dimensional
		obj SimplicialComplex.

		Return a simplicial SimplicialComplex object with k_d * (d+1)! d-cells, where
		k_d is the cardinality of the d-skeleton of obj. 
		"""

		def facenode (face):
			""" To compute the index of the node associated to a face.

			Return an integer number.
			"""
			point = centroid(self,face)
			nodekey = inv_dict[0][ str(point) ]
			return nodekey

		def split_(face):	
			if len(face) > 1:
				return ([ facenode(face), AA(split_)(facets(face)) ])
			return face

		def sign (face):
			""" To compute the sign of a d-simplex in a d-complex.

			Return either 1 or -1 according with the permutation class of the face.
			"""
			vertices = [verts[v] + [1.0] for v in face]
			return cmp(linalg.det(vertices),0)

		# initialize
		obj = self.copy()
		verts = obj.vertices.points
		dim = obj.dim
		n = obj.rn
		cells = obj.cells
		dictos = obj.dictos
		inv_dict = obj.inv_dict
		
		# add new vertices
		base=len(verts)-1
		k = 0
		for skeleton in cells[1:]:
			for face in skeleton:
				p = centroid(self,face)
				verts.append(p)
				k += 1
				operator.setitem(inv_dict[0], str(p), base+k)
				operator.setitem(dictos[0], str([base+k]), base+k )

		# build the new d-cells
		ncells = PROD(INTSTO(dim+1)) #number of splitted cells per d-cell
		oldcells = list(cells[-1])
		cells[-1] = []
		for cell in oldcells:
			store = []
			traversal(split_(cell), [], store)
			cells[-1] += store

		# coherent orientation of the d-cells
		for i in range(len(cells[-1])):
			cell = cells[-1][i]
			if dim == n == len(face) - 1:
				if sign(cell) < 0:
					cells[-1][i] = [cell[-1]] + cell[1:-1] + [cell[0]]			  
			

		return SimplicialComplex(verts, cells[-1])



	## -- offspring Method -------------------------------
	def offspring(self, germs):
		""" To compute a subcomplex of self, with d-cells incident on some 0-cell
		in the list of face lists 'germs'.
		
		Return a new cellComplex with the 0-cells of self and a subset of d-cells.
		"""
		subcomplex = self.copy()
		points = subcomplex.vertices.points
		cells = subcomplex.cells
		
		germ_verts = []
		for germ in germs:
			germ_verts += germ
			
		germ_verts = set(germ_verts)
		offsimplices = [s for s in cells[-1]
						if germ_verts.intersection(s) != set([]) ]	  
		return SimplicialComplex(points, offsimplices)


	## -- subtract Method -------------------------------
	def subtract(self, scomplex):
		""" To compute the subtraction of d-cells of 'scomplex' from self. 'scomplex'
		must be a subcomplex of self.

		Return a new subcomplex of self.
		"""
		newsub = self.copy()
		points = newsub.vertices.points
		cells = newsub.cells

		arg1 = cells[-1]
		arg2 = scomplex.cells[-1]
		new_d_cells = [d_cell for d_cell in arg1 if d_cell not in arg2]
		return SimplicialComplex(points, new_d_cells)



	## -- add Method -------------------------------
	def add(self, scomplex):
		""" To compute the addition of d-cells of 'scomplex' to self. 'scomplex'
		must be a subcomplex of self.

		Return a new subcomplex of self.
		"""
		newsub = self.copy()
		points = newsub.vertices.points
		cells = newsub.cells

		arg1 = cells[-1]
		arg2 = scomplex.cells[-1]
		arg2 = [d_cell for d_cell in arg2 if d_cell not in arg1]
		new_d_cells = arg1 + arg2
		return SimplicialComplex(points, new_d_cells)


	## -- normalize Method -------------------------------
	def normalize(self):
		""" To map the SimplicialComplex so that its containment box is mapped to the
		standard unit cube.

		Return a new subcomplex of self.
		"""
		newsub = self.copy()
		newsub.vertices = newsub.vertices.normalize()
		return newsub


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

##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
		

def simplexGrid (listofquotes):
	""" To generate a simplexGrid of d-simplices.

	Return the complex as a 'cells' list of lists  of simplices.
	"""
	zero_simplex = [[]]
	cells = zero_simplex
	for quote in listofquotes:
		cells = extrude(cells,quote)
	return SimplicialComplex(cells[0],cells[-1])


def chain_complex (obj):
	""" To compute d adjacency matrices that implement the boundary maps.

	Return a list of sparse matrices in CSR (Compressed Sparse Row) format.
	"""
	def incidence_coeff (i,j):
		face = obj.inv_dict[d-1][i]
		faces = map( repr, facets(eval(obj.inv_dict[d][j])) )
		if d == 1:
			if face == faces[0]: value = 1
			elif face == faces[1]: value = -1
		elif face in faces:
			value = 1
		else:
			value = -1
		return value

	cells = obj.cells
	homology = obj.homology
	matrices = []
 
	d = 0
	for pair_set in homology[1:]:
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

## --------------------------------------------------
## --file reading------------------------------------
## --------------------------------------------------


def complex_read(file_name):
	""" To read a '.dat' input file, with CSV syntax.
	Each cell is read from a single line.

	Return the complex as a list of lists 'cells' of simplices.
	"""
	file_name = str(file_name)
	input_file = open(file_name + ".dat", "r")
	reader = csv.reader(input_file)
	cells = [[]]
	
	## ---sequential reading---------
	for row in reader:
		if row[0][0] != "#":
			d = len(cells)
			cell_dim = int(row[0][0])
			n = cell_dim - d + 1
			if row[0] == "0D" and row[-1] == "0D" :
				pass
			elif row[0][0] == "0":
				for i in range(n): cells.append([])
				cells[cell_dim].append(map(float,row[1:]))
			elif row[0][1] == "D":
				for i in range(n): cells.append([])
				cells[cell_dim].append(map(int,row[1:]))

	## ---skeleton extraction--------
	input_file.close()
	return SimplicialComplex(cells[0],cells[-1])


	
def Map (mapping, pol, out=SimplicialComplex):
	""" To Apply the "mappping" to the vertices of the SimplicialComplex pol.
	
	Be careful to write "Map", with the first letter in Upper case.
	Return a new SimplicialComplex.
	"""
	newpoints = map(mapping, pol.vertices.points)
	ptk,verts = remap(newpoints)
	newcells = [AA(ptk)(cell) for cell in pol.cells[-1]]
	""" remove degenerate cells ------------------- """
	
	def simplex(cell):
		point = verts.ind
		return [eval(point[k])+[1.0] for k in cell]
	
	def volume(cell):
		return linalg.det(mat(simplex(cell[:verts.dim+1])))

	newcells = [cell for cell in newcells if (len(cell) == len(set(cell)))]
# 	if verts.dim == pol.rn:  
# 		newcells = [cell  for cell in newcells if volume(cell) != 0.0 ]
				
	return out(newpoints, newcells)



def merge_complex(c1,c2):
	c = SimplicialComplex(c1.vertices.points + c2.vertices.points, c1.cells[-1])
	ind = [c.vertices.dict[repr(v)] for v in c2.vertices.points]
	def coding (k): return coding[k]
	verts = c.vertices.points
	newcells = remove_duplicates(
		c1.cells[-1] + [map(code, cell) for cell in c2.cells[-1]] )
	return SimplicialComplex(verts,newcells)



#//////////////////////////////////////////////////////////

from chompy.chompyScomplex import *


def offset(point,expl=[1,1,1]):
	scaledpoint = [point[k]*expl[k] for k in range(3)]
	vect = VECTDIFF([scaledpoint,point])
	return vect
	

def spheres(batches,points,expl=[1,1,1]):
	sx = 0.05
	points = CAT(points)
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

#/////////////////////////////////////////////////////////

def draw (c,chains=4*[[]],expl=[1,1,1]):

	if c.dim != -1:
		embeddingDim = 3 - c.rn		# possible curation of embedding dimension
		if embeddingDim != 0: c.vertices.embed(embeddingDim)
	
		if len(chains)<(c.dim+1):	# possible curation of input chains
			chains = chains + ((c.dim+1)-len(chains))*[[]]
	
		if chains == [[],[],[],[]]: chains = c.cells # draw the whole complex
		else:										 # draw only the input chains
			chains = [[c.cells[h][k] for k in chains[h]] for h in range(len(chains))]
		
		def cellverts(c,cells):
			return [[eval(c.vertices.ind[v]) for v in cell] for cell in cells]
	
		def addBatches(objType,batches,items,expl):
			return objType(batches,items,expl)
	
		primitives = [spheres, cylinders, planecells, cells]
		batches = []
		for k in range(c.dim + 1):
			if chains[k] != []:
				items = cellverts(c,chains[k])
				batches = addBatches(primitives[k],batches,items,expl)
	
		octree=Octree(batches)
		viewer=Viewer(octree)
		viewer.Run()
	else: print "cannot draw an empty complex!"

#/////////////////////////////////////////////////////////

if __name__ == "__main__":

	c = simplexGrid([2*[1.], 2*[1.], 1*[1.]])
	c.view()
	draw(c,expl=[1,1,2],chains=[[],range(len(c.cells[1]))])

