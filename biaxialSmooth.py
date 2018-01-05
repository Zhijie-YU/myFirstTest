""" Author: Ning Guo <ceguo@connect.ust.hk>
    run `mv biaxialSmooth.yade.gz 0.yade.gz`
    to generate initial RVE packing
"""
from esys.escript import *
from esys.finley import Rectangle
from esys.weipa import saveVTK
from esys.escript.pdetools import Projector
from msFEM2D import MultiScale
from saveGauss import saveGauss2D
import time


vel = -0.0001; confining=-1.e5;
lx = 0.05; ly = 0.1; # sample size, 50mm by 100mm
nx = 2; ny = 4; # sample discretization, 2 by 4 quadrilateral elements

mydomain = Rectangle(l0=lx,l1=ly,n0=nx,n1=ny,order=2,integrationOrder=2)#rectangle is defined in escript

dim = mydomain.getDim()#escript,returns the number of spatial dimensions of the Domain.Here returns 2
k = kronecker(mydomain)#Here returns ([1,0],[0,1])
numg = 4*nx*ny; # number of Gauss points, 4 GP each element (reduced integration)
nump = 4; # number of processes for multiprocessing

prob = MultiScale(domain=mydomain,ng=numg,np=nump,random=False,rtol=1e-2,usePert=False,pert=-2.e-5,verbose=True)

disp = Vector(0.,Solution(mydomain)) #Solution is from escript, returns solutions of a PDE;??
t=0

#fout2=open('./result/displacement.dat','w')
#fout2.write(str(t)+' '+str(disp)+'\n')

stress = prob.getCurrentStress() # initial stress. getCurrentStress is imported from msFEM2D.py
proj = Projector(mydomain) #from escript. 
sig = proj(stress) # project Gauss point value to nodal value
sig_bounda = interpolate(sig,FunctionOnBoundary(mydomain)) 
# interpolate(a, where): interpolates argument a into the FunctionSpace where.
#from escript. FunctionOnBoundary(domain): returns the boundary FunctionSpace on the Domain domain. Data objects in this type of general FunctionSpace are defined on the boundary of the geometric region defined by domain.
traction = matrix_mult(sig_bounda,mydomain.getNormal()) # boundary traction
#getNormal(): If the domain of functions in the FunctionSpace is a hyper-manifold (e.g. the boundary of a domain),the method returns the outer normal at each of the data sample points. Otherwise an exception is raised.
#matrix_mult(a0, a1): returns the matrix product of a0 and a1.
x = mydomain.getX() # nodal coordinate;returns the locations in the Domain. x[0] means x and x[1] means y for 2D
bx = FunctionOnBoundary(mydomain).getX()  #??why is the value of bx not consistent with x??

topSurf = whereZero(bx[1]-sup(bx[1]))
#sup(a): returns the maximum value over all components and all data sample points of a.
tractTop = traction*topSurf # traction at top surface,'stress'
forceTop = integrate(tractTop,where=FunctionOnBoundary(mydomain)) # resultant force at top
lengthTop = integrate(topSurf,where=FunctionOnBoundary(mydomain)) # length of top surface, here returns 0.05
#integrate(a[ , where=None ]) returns the integral of a, where the domain of integration is defined by the FunctionSpace of a.
fout=open('./result/biaxial_surf.dat','w')  #where is the file??./we should create folder 'result' by ourseles.
fout.write('0 '+str(forceTop[1])+' '+str(lengthTop)+'\n')
'''
print('topSurf=',topSurf)
print('traction=',traction)
print('tractTop=',tractTop)
print('forceTop=',forceTop)
print('lengthTop=',lengthTop)
exit()
'''
# Dirichlet BC positions, smooth at bottom and top, fixed at the center of bottom
#define where the displacement boundary is,x[0] means x,x[1] denotes y, for the 1st [0,1], 0 means no constraint in x direction,1 means there is constraint in y direction
Dbc = whereZero(x[1])*[0,1]+whereZero(x[1]-ly)*[0,1]+whereZero(x[1])*whereZero(x[0]-.5*lx)*[1,1] 
# Dirichlet BC values
#define how much the dis boundary is,[0,vel] means no displacement in x direction and the disp in y direction is vel every step
Vbc = whereZero(x[1])*[0,0]+whereZero(x[1]-ly)*[0,vel]+whereZero(x[1])*whereZero(x[0]-.5*lx)*[0,0]
# Neumann BC, constant confining pressure
#define the stress boundary, confine denotes stress, and right is positive for x direction
#bx denotes the boundary
Nbc = whereZero(bx[0])*[-confining,0]+whereZero(bx[0]-lx)*[confining,0] 
#for boundary displace and stress, positive/negative is consistent with coordinate.
time_start = time.time()#from python, return current time e.g. 1507688377.83. the value means how many seconds has passed since 1970.1.1
while t < 0: # apply 100 load steps
   
   timebegin=time.time()
   prob.initialize(f=Nbc, specified_u_mask=Dbc, specified_u_val=Vbc) # initialize BC
   t += 1
   du=prob.solve(iter_max=100) # get solution: nodal displacement:!!!Key Calculation Step!!!

   disp += du
   stress=prob.getCurrentStress()# stress on GP
   
   #fout2.write(str(t)+' '+str(disp)+'\n')
 
   dom = prob.getDomain() # domain is updated Lagrangian formulation
   proj = Projector(dom)
   sig = proj(stress)# project Gauss point value to nodal value

   sig_bounda = interpolate(sig,FunctionOnBoundary(dom))
   traction = matrix_mult(sig_bounda,dom.getNormal())
   tractTop = traction*topSurf
   forceTop = integrate(tractTop,where=FunctionOnBoundary(dom))
   lengthTop = integrate(topSurf,where=FunctionOnBoundary(dom))
   time1iter=time.time()-timebegin
   fout=open('./result/biaxial_surf.dat','a')
   fout.write(str(t*vel/ly)+' '+str(forceTop[1])+' '+str(lengthTop)+' '+str(time1iter)+'s'+'\n')
      
   vR=prob.getLocalVoidRatio()
   fabric=prob.getLocalFabric()
   strain = prob.getCurrentStrain()
   saveGauss2D(name='./result/gauss/time_'+str(t)+'.dat',strain=strain,stress=stress,fabric=fabric)
   volume_strain = trace(strain)
   dev_strain = symmetric(strain) - volume_strain*k/dim
   shear = sqrt(2*inner(dev_strain,dev_strain))
   saveVTK("./result/vtk/biaxialSmooth_%d.vtu"%t,disp=disp,shear=shear,e=vR)

   

prob.getCurrentPacking(pos=(),time=t,prefix='./result/packing/')
time_elapse = time.time() - time_start
fout.write("#Elapsed time in hours: "+str(time_elapse/3600.)+'\n')   
fout.close()
prob.exitSimulation()
