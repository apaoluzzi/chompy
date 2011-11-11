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

from chompyUtil import *

#########################################################
class PointSet(list):  ## NEW
#########################################################
    """ The data type to represent a set of points,
        supporting a cell complex. """


    ## -- __init__ Method --------NEW-----------------------
    def __init__(self, points):
                            
        self.dict = {}
        k = -1
        points = [ AA(float)(p) for p in points ]
        points = [eval(code(p)) for p in points]
        for i in range(len(points)):
            if code(points[i]) not in self.dict:
                k += 1
                operator.setitem(self.dict,code(points[i]), k)
                
        self.ind = dict([[v,k] for k,v in self.dict.items()])
#        self.points = points
        self.points = AA(eval)(self.ind.values())
#        self.points = [eval(self.ind[k]) for k in range(len(self.dict))]
        self.dim = len(self.points[0])
        self.m = len(self.dict)


    ## -- __repr__ Method -------------------------------
    def __repr__(self):
        '''Return a string representation of this SimplicialComplex
            in the form of a list of cells.'''
        return 'Pointset: %s' % self.ind

    ## -- __str__ Method -------------------------------
    def __str__(self):
        '''Return a string representation of this SimplicialComplex
            in the form of a list of cells.'''
        indices = self.ind.keys()
        pts = [eval(self.ind[k]) for k in indices]
        return 'Pointset: {0} x {1}: {2}'.format(self.m, self.dim, pts)

    ## -- copy Methods -------------------------------
    def copy(self):
        """ To make a copy of a PointSet object. """
        return copy.deepcopy(self)

    ## -- insert Method -------------------------------
    def insert(self,point):
        """ To add a point (if not already present) to a PointSet object. """
        key = code(point)
        try: key in self.dict[key]
        except KeyError: 
            self.m += 1
            self.dict[key] = self.m
            self.ind[self.m] = key
            self.points = self.points + [eval(code(point))]
        return self

    def update(self,modify):
        """ To update a PointSet object without change the element ordering. """
        index = dict([[k,modify(v)] for k,v in self.ind.items()])
        self.ind.update(index)
        self.dict = dict([[v,k] for k,v in self.ind.items()])
        self.points = [eval(self.ind[v]) for v in range(self.m)]
        return self

    ## -- embed Method -------------------------------
    def embed (self, n):
        def modify(v): return code(eval(v) + n*[0.0])
        self.update(modify)
        self.dim += n
        return self
    
    def transform (self, matrix):
        def modify(v): return code(mat(matrix) * mat(AA(LIST)(eval(v))))
        self.update(modify)
        return self
    
    def translate (self, vect):
        def modify(v): return code(array(eval(v)) + vect)
        self.update(modify)
        return self
        
    def scale (self, vect):
        def modify(v): return code(array(eval(v)) * vect)
        self.update(modify)
        return self

    def rotate (self, axis1, axis2, angle):
        axis1 += -1
        axis2 += -1
        dim = self.dim
        rotation = eye(dim,dim)
        c = math.cos(angle)
        s = math.sin(angle)
        rotation[axis1,axis1] = c
        rotation[axis2,axis2] = c
        rotation[axis1,axis2] = -s
        rotation[axis2,axis1] = s
        def modify(v):
            return code((mat(eval(v)) * mat(rotation).T).tolist()[0])
        self.update(modify)
        return self

    def max (self):
        return map(max, array(self.points).T)

    def min (self):
        return map(min, array(self.points).T)

    def Med (self):
        Max = array(map(max, array(self.points).T))
        Min = array(map(min, array(self.points).T))
        return (Max + Min)/2

#--TODO---------------------------------------
    def normalize (self):
        Min = array(map(min, array(self.points).T))
        Max = array(map(max, array(self.points).T))
        translation = -Min
        size = Max-Min
        scaling = array([1/s if (s != 0) else 1.0 for s in size])
        self = self.translate(translation.tolist())
        self = self.scale(scaling.tolist())
        return self
    
##<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


def remap(points):
	verts = PointSet(points)
	def remap0(k): return verts.dict[code(points[k])]
	return remap0,verts


##if __name__ == "__main__" and __debug__:
##    covectors = S1(UKPOLF(PROD([QUOTE([1,1,1]),QUOTE([1,1,1])])))
##    covectors = PointSet(covectors)
##    print "\ncovectors =", covectors


def embed (pointset, n):
    added_coords = n * [0.0]
    points = map(lambda p: p + added_coords, pointset.points)
    return PointSet(points)


def subcomplex (d,List):
    """ To compute the list of adjacent d-tuples obtained by shifting the
    list List by one element.
    
    Return the list of n d-lists, with n = len(List)-d+1
    """
    return [List[i:i+d] for i in range(len(List)-d+1)]


def shift (n, listoflists):
    """ To increment by (int) 'n' all the (int) elements of a list of lists.

    Returns a list of lists of int.
    """
    return [map(lambda x: x+n, list) for list in listoflists]


def extrude (cells, hlist):
    """ To complute the multiple linear extrusion of a d-complex.
    Map R^d -> R^(d+1), according to: Ferrucci & Paoluzzi, CAD 1991.
    'cells' is a simplicial complex, given as a list of lists;
    'hlist' is a list of heights in the added dimension.

    Return the (d+1)-dimensional simplicial complex.
    Only the 0- and (d+1)- skeletons are computed.
    In order to get the remaining skeletons, the 'cell-complex' function
    must be invoked.
    """
    dim = dimension(cells)
    if dim == 0:
        cells = [[],[]]
        lastcoords = progressive_sum(hlist)
        cells[0] = AA(LIST)(lastcoords)
        cells[1] = [[i,i+1] for i in range(0,len(hlist))]
    else:
        verts = cells[0]
        simplexes = cells[dim]
        nverts = len(verts)
        lastcoords = progressive_sum(hlist)
        nsteps = len(lastcoords)
        sliced_vertices = nsteps*[verts]
        coords_distribute = lambda x: \
            COMP([ CAT, AA(COMP([ AA(AR), DISTR ])) ])(x)
        cells[0] = coords_distribute(TRANS([nsteps*[verts],lastcoords]))
        extruded_simplexes = []
        for cell in simplexes:
            I = cell + map(lambda x: x+nverts, cell)
            LL = subcomplex(dim+2,I)
            extruded_simplexes += LL
        final_simplices = []
        for i in range(nsteps-1):
            simplex_layer = shift(nverts*i,extruded_simplexes)
            final_simplices += simplex_layer
        cells.append(final_simplices)
    return cell_complex(cells)

# def extrude (cells, hlist):
# 	""" To complute the multiple linear extrusion of a d-complex.
# 	Map R^d -> R^(d+1), according to: Ferrucci & Paoluzzi, CAD 1991.
# 	'cells' is a simplicial complex, given as a list of lists;
# 	'hlist' is a list of heights in the added dimension.
# 
# 	Return a (d+1)-dimensional simplicial complex.
# 	Only the 0- and (d+1)- skeletons are computed.
# 	In order to get the remaining skeletons, the 'SimplicialComplex' function
# 	is used.
# 	"""
# 
# 	dim = dimension(cells)
# 	verts = cells[0]
# 	myprint("hlist",hlist)
# 	lastcoords = progressive_sum(AA(ABS)(hlist))
# 	
# 	if dim == 0:
# 		cells = [[],[]]
# 		vertices = AA(LIST)(lastcoords)
# 		cells[1] = [[i,i+1] for i in range(1,len(hlist)+1)]
# 	else:
# 		simplexes = cells[dim]
# 		nverts = len(verts)
# 		nsteps = len(lastcoords)
# 		sliced_vertices = nsteps*[verts]
# 		def coords_distribute(x):
# 			return COMP([ CAT, AA(COMP([ AA(AR), DISTR ])) ])(x)
# 		vertices = coords_distribute(TRANS([nsteps*[verts],lastcoords]))
# 
# 		extruded_simplices = []
# 		for cell in simplexes:
# 			vertPtrs = cell + map(lambda x: x+nverts, cell)  
# 			extruded_simplices += subcomplex(dim+2,vertPtrs)
# 	
# 		final_simplices = []
# 		for i in range(nsteps-1):
# 			if hlist[i] > 0:
# 				simplex_layer = shift(nverts*i,extruded_simplices)
# 				final_simplices += simplex_layer
# 		cells.append(final_simplices)
# 		
# 	return cell_complex(cells)

## --------------------------------------------------
## --Format conversion-------------------------------
## --------------------------------------------------

def complex2surface (obj):
    """ Transforms the 2-skeleton of a simplicial complex in
    a triangulated surface.

    Return a list of triples of surface vertices (nD points).
    """
    cells = obj.cells
    verts = obj.vertices.points
    surface = [[verts[i],verts[j],verts[k]] for [i,j,k] in cells[2]]
    return surface


def pfacets (pcomplex):
    def pfacets0 (cell):
        """ To compute the (non oriented) boundary of a d-polytope.
        Both the d-polytope and its (d-1)-facets are represented as
        an ordered list of vertices.
        Return a list of facets, i.e. of (d-1)-polytopes.
        """
        def pack(cells):
            return sorted(AA(eval)(list(set(AA(repr)(cells)))))

        verts = [eval(pcomplex.vertices.ind[v]) for v in cell]

        pol = Plasm.mkpol(len(verts[0]),CAT(verts),[range(len(cell))])
        pol = Plasm.skeleton(pol, Plasm.getPointDim(pol)-1)

        v = StdVectorFloat()
        cells = StdVectorStdVectorInt()
        dim = Plasm.ukpol(pol, v, cells)
        points = [[v[h] for h in range(i,i+dim)]
                  for i in xrange(0, len(v), dim)]

        facets = [pack([pcomplex.vertices.dict[code(points[k])] for k in cell])
                  for cell in cells]
        return facets
    return pfacets0

