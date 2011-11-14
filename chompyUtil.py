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

__debug = False

ROUND_ZERO = 1E-06

from numpy import array,average
from scipy import sparse,mat,eye,linalg
from pyplasm import *
import operator,copy,datetime,chompyIntegr,scipy,csv
from chompyIntegr import T3
import itertools,string

def myprint(arg1,arg2): print "\n"+arg1+" =", arg2

def round_or_zero (x,prec=6):
	"""
	Decision procedure to approximate a small number to zero.

	Return either the input number or zero.

	"""
	def myround(x):
		return eval(('%.'+str(prec)+'f') % round(x,prec))
	xx = myround(x)
	if abs(xx) < ROUND_ZERO: return 0.0
	else: return xx

def evalprint__(string):
	"""
	Simple tool to evaluate a language expression passed as a string.

	Side effect: print the string and its value.
	Return ``None``. 

	"""
	print "\n" + string + " => ", eval(string)

def prepKey (args): return "["+", ".join(args)+"]"

def fixedPrec(value):
	if abs(value - int(value))<1E-5: value = int(value)
	out = ('%0.6f'% value).rstrip('0')
	if out == '-0.': out = '0.'
	return out
	
	
def code (vect): 
	"""
	To generate a string representation of a number array.
	
	Used to generate the vertex keys in PointSet dictionary, and other similar operations.
	"""
	return prepKey(AA(fixedPrec)(vect))

		
def cat(listOfLists): 
	"""
	To catenate a list of lists (flat:  no recursion).
	
	Return a list.
	"""
	return [element for list in listOfLists for element in list]
		
	
def atan2(point):
	"""
	To compute the ATAN function starting from COS and SIN.
	
	Return a number between between the interval [-pi,pi].
	"""
	x,y = point
	return math.atan2(x,y)


## --------------------------------------------------
## --Utility functions-------------------------------
## --------------------------------------------------	 

def __evalprint__(string):
	"""
	Simple tool to evaluate a language expression passed as a string.

	Side effect: print the string and its value.
	Return ``None``. 

	"""
	print "\n" + string + " => ", eval(string)


def printer (args):
	""" To make an 'internal printing' of a (possibly nested) list in plasm format.

	Return a string, with squaare brackets substituted by angle brackets
	"""
	def one (list_of_strings):
		"""To catenate a flat list of string.

		Return a single string in plasm format.
		"""
		return '<' + ', '.join(list_of_strings) + '>'
	
	if isinstance(args, list):
		return one(map(printer, args))
	else:
		return str(args)


def progressive_sum(List):
	""" To compute the coordinates associated to a list of (Num) differences.
	Starting from zero.

	Return a list of Num
	"""
	return [sum(List[0:i]) for i in range(len(List)+1)]


def centroid (obj,face):
	""" To compute the centroid od a d-face.

	`face` is the canonical representation of face (list of vertex indices)
	Return a n-point, convex combination of d+1 n-points.
	"""
	simplex = [ obj.vertices.points[v]  for v in face ]
	d = len(simplex)
	A = mat(simplex, dtype=scipy.float32)
	C = mat(d*[1.0/d], dtype=scipy.float32)
	point = (C * A).tolist()[0]
	return point
	
def Centroid (args): return centroid (*args)


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

## --------------------------------------------------
## --Skeletons of a complex--------------------------
## --------------------------------------------------

def dimension(cells):
	""" To compute the (intrinsic) dimension of a cell complex.

	Return the last non zero cardinality of k-skeletons, with
	k <= dim.
	"""
	k = map(len, cells)
	dim = len(cells) - 1
	while k[dim] == 0 and dim > 0:
		dim = dim - 1
	return dim

def ispair (arg):
	return len(arg) == 2

def issingleton (arg):
	return len(arg) == 1

def traversal(tree, path, store):
	
	path = path + [tree[0]]
	
	if ispair(tree):
		for item in tree[1]:
			traversal(item, path, store)
			
	elif issingleton(tree):
		store.append(path)
		path = []


def remove_duplicates (hasDupes):
	""" To remove the duplicates from a list of cells.
	Consider as duplicate the cell with opposite orientation (using the
	canonical representation)
	The mapping on a dictionary is used, with key the repr of the elements.

	Returns a new list without duplicates.
	"""
	def revert(cell):
		""" Change the sign of a permutation by exchanging the first
		and the last elements
		"""
		if len(cell) > 1:
			cell = [cell[-1]] + cell[1:-1] + [cell[0]]
		return cell

	noDupes = {}
	[operator.setitem(noDupes,repr(cell), cell) for cell in hasDupes \
		 if not (noDupes.has_key(repr(cell))  or \
				 noDupes.has_key(repr(revert(cell))))]
	return noDupes.values()



def skeleton (h,cells):
	""" To compute the h-skeleton of a complex (given as a list of lists of 'cells').

	Return the list of h-cells, without duplicates.
	"""
	if h < len(cells):
		A = []
		for cell in cells[h]:
			A.extend(facets(cell))
	return remove_duplicates(A)



def cell_complex (cells):
	""" To compute and store all the skeletons of a cell complex.
	The input data structure is only required to contain the 0- and the
	highest dimensional skeletons (0- and n-, with dim == n)

	Return the (possibly upgraded) list of lists of k-dimensional cells,
	with 0 <= k <= dim.
	"""	   
	dim = dimension(cells)
	k = map(len, cells)
	for h in range(dim,1,-1):
		cells[h-1] = skeleton(h,cells)
		k[h-1] = len(cells[h-1])
	return cells


def mktables (Complex):
	""" To build the cell database of a cell complex.

	Return a list of d+1 dictionaries, where d is the dimension of the complex.
	Key is the 'repr' of the cells; values are integers in range(0,k[d]) where
	k[d] is the cardinality of the d-skeleton.
	"""
	# resetting dictionaries
	# (so that each call to mktables results in a new set of dictionaries)
	dictos = []
	
	def rotate(cell):
		""" Rotate a permutation into another in the same permutation class
		"""
		m = min(cell)
		local = list(cell)
		while local[0] != m:
			local = local[1:] + [local[0]]
		cell = local
		return cell
	
	def revert(cell):
		""" Change the sign of a permutation by exchanging the first and the last elements
		"""
		if len(cell) > 1:
			cell = [cell[-1]] + cell[1:-1] + [cell[0]]
		return cell

	d = 0
	for skeleton in Complex[:]:
		dictos.append({})
		[operator.setitem(dictos[d], repr(cell), skeleton.index(cell)) \
			for cell in skeleton \
				if not (dictos[d].has_key(repr(cell)) or \
						 dictos[d].has_key(repr(revert(cell))))]
		d += 1
	return dictos

def facets (cell):
	""" To compute the (non oriented) boundary of a d-simplex.
	Uses the standard method of algebraic topology:
	elimination of a vertex index to generate every face.
	Both the simplex and its facets are represented as an ordered list of vertices.

	Return a list of facets, i.e. of (d-1)-dimensional simplexes.
	"""
	def rotate(cell):
		""" Rotate a permutation into another in the same permutation class
		"""
		m = min(cell)
		local = list(cell)
		while local[0] != m:
			local = local[1:] + [local[0]]
		cell = local
		return cell

	def revert(cell):
		""" Change the sign of a permutation by exchanging the first and the last elements
		"""
		if len(cell) > 1:
			cell = [cell[-1]] + cell[1:-1] + [cell[0]]
		return cell

	faces = []
	# extracts facets of odd index ------------
	for i in range(1, len(cell), 2):
		v = cell[i]
		face = list(cell)
		face.remove(v)
		faces.append(rotate(face))


	# extracts facets of even index ------------
	for i in range(0, len(cell), 2):
		v = cell[i]
		face = list(cell)
		face.remove(v)
		faces.append(revert(rotate(face)))
	return faces

if __name__ == "__main__" and __debug__:
	print "\n --boundary operator on a single cell-----------"
	__evalprint__(""" facets([1,2]) """)
	__evalprint__(""" facets([1,2,3]) """)
	__evalprint__(""" facets([1,2,3,4]) """)
	__evalprint__(""" facets([1,2,3,4,5]) """)


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

	global homology
	homology = [[],[]]
	skeleton = dictos[1]
	for cell in skeleton:
		homology[1].extend([(skeleton[cell], facet[0], ) \
		   for facet in facets(eval(cell))])
	d = 1
	for skeleton in dictos[2:]:
		homology.append([])
		d += 1
		for cell in skeleton:
			for facet in facets(eval(cell)):
				try:
					key = dictos[d-1][repr(facet)]
					homology[d].extend([(skeleton[cell], key)])
				except KeyError:
					key = dictos[d-1][repr(revert(facet))]
					homology[d].extend([(skeleton[cell], key)])
	return homology


