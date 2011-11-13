from chompy.chompyDcomplex import *

cprod = INSR(complexProd)

def idnt(d):
	"""
	Constructor of the identity matrix with d rows.
	
	Return a list of lists of integers.
	"""
	return [k*[0]+[1]+(d-k-1)*[0] for k in range(d)]

def cuboid(sides):
	"""
	Primitive constructor of a d-dimensional cuboid.
	
	"""
	unit = SimplicialComplex([[0.],[1.]],[[0,1]])
	d = len(sides)
	return cprod(d * [unit]).scale(sides)
	
def simplex(d):
	"""
	Primitive constructor of a d-dimensional simplex.
	
	"""
	vertices = [d*[0]] + idnt(d)
	cells = [range(d+1)]
	return SimplicialComplex(vertices, cells)


def polyline(points):
	"""
	Primitive constructor of a polygonal line.
	
	"""
	cells = [[k,k+1] for k in range(len(points)-1)]
	return SimplicialComplex(points,cells)

def triangle_strip(points):
	"""
	Primitive constructor of a strip of triangles.
	
	"""
	ptk,verts = remap(points)
	cells = [AA(ptk)([k,k+1,k+2]) for k in range(len(points)-2)]	
	return SimplicialComplex(points,cells)

def triangle_array(m,n,points):
	"""
	Primitive constructor of a 2-indexed array of points.
	
	"""
	out = simplexGrid([m*[1.],n*[1.]])
	return SimplicialComplex(CAT(points),out.cells[2])


def cart2cyl2d(point):
	u,v = point
	return [COS(u),SIN(u),v]



def cylsurface(r=1,h=1,n=16,m=2):
	return Map(cart2cyl2d, simplexGrid([n*[2*PI/n], m*[1.0/m]])).scale([r,r,h])
	


def cart2cyl3d(point):
	u,v,w = point
	return [v*COS(u), v*SIN(u), w]


def cylsolid(r=1,h=1,n=16,m=1,p=1):
	domain = simplexGrid([n*[2*PI/n],m*[1.0/m],p*[1.0/p]])
	return Map(cart2cyl3d, domain).scale([r,r,h])


def cart2torus2d(r,R):
	def cart2torus2d0(point):
		u,v = point
		return [(R+r*COS(v))*COS(u),(R+r*COS(v))*SIN(u),r*SIN(v)]
	return cart2torus2d0


def torus_surface(r=1,R=3,n=12,m=8):
	return Map(cart2torus2d(r,R), simplexGrid([n*[2*PI/n],m*[2*PI/m]]))


def cart2torus3d(r,R):
	def cart2torus3d0(point):
		u,v,w = point
		return [(R+r*w*COS(u))*COS(v),
				(R+r*w*COS(u))*SIN(v),
				r*w*SIN(u)]
	return cart2torus3d0

def torus_solid(r=1,R=3,n=8,m=16,p=1):
	return Map(cart2torus3d(r,R), simplexGrid([n*[2*PI/n],m*[2*PI/m],p*[1.0/p]]))


def schlegel(pol):
    def project(point):
        return [coord/point[-1] for coord in point[:-1]]
    verts = [project(v) for v in pol.vertices.points]
    cells = pol.cells[-2]
    return PolytopalComplex(verts,cells)


def intervals(tip):
    def intervals0(n):
        points = [[x] for x in scipy.linspace(0.0, tip, n+1)]
        return polyline(points)
    return intervals0


def graph(domain):
    def graph0(funs):
        def vfun(point): return CONS(funs)(point[0])
        return Map(vfun,domain)
    return graph0


def circumpherence(r,nsides=24):
    return graph(intervals(2*math.pi)(nsides))([math.cos,math.sin]).scale([r,r])


def helix(radius=1,pitch=1,n=24,turns=1):
    return graph(intervals(2*PI*turns)(n*turns))(
        [SIN,COS,RAISE(DIV)([ID,K(2*PI/pitch)])]).scale([radius,radius,1])




if __name__=="__main__":
	
	points = [[0,3],[1,2],[3,3],[2,2],[3,0],[2,1],[0,0],[1,1],[0,3],[1,2]]
	pol = triangle_strip(points).extrude([1,-1,1]).boundary()
	draw(pol, expl=[1,1,3])#,
		 #chains=[[],[],range(len(pol.cells[2]))])
