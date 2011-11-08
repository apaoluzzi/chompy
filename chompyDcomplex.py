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

from chompyCcomplex import *


#########################################################
class CovertexSet(PointSet):  
#########################################################
    """ The data type to represent a set of points,
        supporting a cell complex. """


    ## -- __init__ Method -------------------------------
    def __init__(self, coverts):
                            
        self.dict = {}
        k = -1
        coverts = [ AA(float)(p) for p in coverts ]
        coverts = [ AA(round_or_zero)(p) for p in coverts ]
        for i in range(len(coverts)):
            if coverts[i][0] < 0.0:
                coverts[i] = SCALARVECTPROD([-1,coverts[i]])
        for i in range(len(coverts)):
            if code(coverts[i]) not in self.dict:
                k += 1
                operator.setitem(self.dict,code(coverts[i]), k)
                
        self.ind = dict([[v,k] for k,v in self.dict.items()])
        self.points = [eval(self.ind[k]) for k in range(len(self.dict))]
        self.dim = len(self.points[0])
        self.m = len(self.dict)


#########################################################
class CofaceComplex(PolytopalComplex):  
#########################################################

            
    ## -- __init__ Method -------------------------------
    def __init__(self, coverts,cocells=[],dim=3):

        def convert(coverts, cocells):
            cocells = [[ coverts[f-1] for f in cocell] for cocell in cocells]
            out = []
            for k in range(len(cocells)):
                for j in range(len(cocells[k])):
                    if cocells[k][j][0] < 0.0:
                        cocells[k][j] = SCALARVECTPROD([-1,cocells[k][j]])
                out += [[self.vertices.dict[code(p)] for p in cocells[k]]]
            return out

        def pack(cofaces):
            return AA(eval)(list(set(AA(repr)(cofaces))))

        self.vertices = CovertexSet(coverts)
        self.dim = dim
        self.rn = len(coverts[1])

        self.cells = (self.dim + 1)*[[]]
        self.n = len(self.vertices.dict)

        self.cells[-1] = convert(coverts, cocells)
        self.cells[-2] = [[f] for f in self.vertices.ind]
                        


#########################################################
class DualComplex(PolytopalComplex):  
#########################################################
    """ Input: a hierarchical polyhedral complex (HPC data structure) """

    def pointMembership (p, coface):
        return INNERPROD([[1.0]+p,coface])==0.0
    
            
    ## -- __init__ Method -------------------------------
    def __init__(self, hpc):


        def hpcExtract(obj,Plasm_fun,incr):
            v = StdVectorFloat()
            u = StdVectorStdVectorInt()
            vdim = Plasm_fun(obj, v, u) + incr
            feature = []
            for i in xrange(0, len(v), vdim):
                feature += [[v[i] for i in range(i,i+vdim)]]
            return feature,[k for k in u]

        verts,faces = hpcExtract(hpc,Plasm.ukpol,0)
        
        self.primal = PolytopalComplex(verts,faces)
        self.primal.view()
        dim = self.primal.dim

        hfaces,u = hpcExtract(hpc,Plasm.ukpolf,1)
        cofaces = CovertexSet(hfaces)

        verts,faces = hpcExtract(Plasm.skeleton(hpc,dim-1),Plasm.ukpol,0)
        facesByVerts = [[verts[w] for w in face] for face in faces]
        newfaces = [sorted([self.primal.vertices.dict[code(w)] for w in face])
                     for face in facesByVerts]

        def normalize (face):
            if round_or_zero(face[0])<0.0: face = SCALARVECTPROD([-1,face])
            return face
        def primalFacetIndex(face): return self.primal.dictos[dim-1][repr(face)]
        def dualFacet(i):
            return cofaces.dict[code(normalize(hfaces[i]))]
        
        pairs = sorted([( primalFacetIndex(face),dualFacet(i) )
                    for i,face in enumerate(newfaces)
                    if repr(face) in self.primal.dictos[dim-1]])

        # remove duplicates
        self.dualCocells = [pairs[k] for k in range(len(pairs)-1)
                       if pairs[k] != pairs[k+1]] + [pairs[-1]]
        dualCocells = [[pair[1]] for pair in self.dualCocells]

        self.dual = PolytopalComplex(cofaces.points, dualCocells)

        

    ## -- __repr__ Method -------------------------------TODO:
    def __repr__(self):
        '''Return a string representation of this CellComplex in the form of a list of cells.'''
        return ('\nDual polytope: \n\t.primal.n: {0}'+ 
                        '\n\t.primal.dim: {1}'+ 
                        '\n\t.dual.n: {2}'+ 
                        '\n\t.dual.dim: {3}'+ 
        '').format(self.primal.n, self.primal.dim, self.dual.n, self.dual.dim)


    ## -- __str__ Method -------------------------------

    def __str__(self):
        '''Return a string representation of this CellComplex in the form of a list of cells.'''
        return ('\nDual polytope: \n\t.primal.vertices: {0}'+ 
                        '\n\t.primal.cells: {1}'+
                        '\n\t.dual.vertices: {2}'+ 
                        '\n\t.dual.cells: {3}'+
                        '\n\t.primal_dual.cocells: {4}'+
        '').format(self.primal.vertices, self.primal.cells[-2],
                   self.dual.vertices, self.dual.cells[-1],
                   self.dualCocells)

        


