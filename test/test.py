from chompy import *

def randomPoints(d=2,n=40,scale=2):
    return PointSet([[random.random()*scale for k in range(d)]
                     for point in range(n)])

if __name__=="__main__":

    rpoints = randomPoints(d=3,n=40,scale=2)
    points = rpoints.points
    myprint("points",points)
    draw(SimplicialComplex(points,AA(LIST)(range(len(points)))))
