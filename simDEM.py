__author__="Ning Guo, ceguo@connect.ust.hk"

"""
DEM part for multiscale simulation which sets up a packing representing
a material point (RVE) at the Gauss point of the FEM domain and returns
constitutive responses at this point.
"""
# import yade modules
import sys
sys.path.append('/home/xiaoyu/yadeNew/install/bin') # path where you have yadeimport.py!!! this should be changed under different conditions
from yadeimport import *
import numpy

def initLoad(ID=0): # where ID identifies the Gauss point location
   if 1:
      # All Gauss points import 0.yade.gz resulting in a uniform sample
      Omega().load('0.yade.gz')
   else:
      # Otherwise load different packings to generate random field
      # resulting in an inherently heterogeneous sample
      Omega().load(str(ID)+'.yade.gz')
   Omega().tags['id']=str(ID)
   return Omega().sceneToString()

def outputPack(param):
   if len(param) != 3:
      raise RuntimeError,"No. of param. should be exactly 3. 0: RVE scene; 1: step count; 2: name prefix"
   Omega().stringToScene(param[0])
   pos = Omega().tags['id'] #this tags is defined in 'initLoad' above
   Omega().save(param[2]+'packing_'+pos+'_'+str(param[1])+'.yade.gz')

def numOfParticles(scene):
   Omega().stringToScene(scene)
   return len(Omega().bodies) # !!! for spherical particle only

# Apply deformation on DEM packing
def shear2D(param):
   if len(param) != 2:  #len(a),for matrix a, this returns the num of rows
      raise RuntimeError,"No. of param. should be exactly 2. 0: RVE scene; 1: strain."
   Omega().stringToScene(param[0])
   ns=int(max(1e5*numpy.max(numpy.abs(param[1])),2))
   dstrain = utils.Matrix3(param[1][0],param[1][1],0, param[1][2],param[1][3],0, 0,0,0)
   Omega().cell.velGrad=dstrain/(ns*Omega().dt)
   Omega().run(ns,True)
   Omega().cell.velGrad=utils.Matrix3.Zero
   return Omega().sceneToString()

def shear3D(param):
   if len(param) != 2:
      raise RuntimeError,"No. of param. should be exactly 2. 0: RVE scene; 1: strain."
   Omega().stringToScene(param[0])
   ns=int(max(1e5*numpy.max(numpy.abs(param[1])),2))
   dstrain = utils.Matrix3(param[1][0],param[1][1],param[1][2], param[1][3],param[1][4],param[1][5], param[1][6],param[1][7],param[1][8])
   Omega().cell.velGrad=dstrain/(ns*Omega().dt)
   Omega().run(ns,True)
   Omega().cell.velGrad=utils.Matrix3.Zero
   return Omega().sceneToString()

""" Used for perturbation method only (2D)
def getStressTensor(scene):
   Omega().stringToScene(scene)
   stress = utils.getStress()
   stress = .5*(stress+stress.transpose())
   return [[stress[0],stress[1]],[stress[3],stress[4]]]
"""

def getFabric2D(scene):
   Omega().stringToScene(scene)
   f = utils.fabricTensor(splitTensor=False)[0]
   return [[f[0,0],f[0,1]],[f[1,0],f[1,1]]]

def getFabric3D(scene):
   Omega().stringToScene(scene)
   f = utils.fabricTensor(splitTensor=False)[0]
   return [[f[0],f[1],f[2]],[f[3],f[4],f[5]],[f[6],f[7],f[8]]]

""" Used for clumped particle model only
def getParOriFabric(scene):
   Omega().stringToScene(scene)
   fab = utils.Matrix3.Zero
   numPar = 0
   for b in Omega().bodies:
      if b.isClump:
         numPar += 1
         keys = b.shape.members.keys()
         pos1 = Omega().bodies[keys[0]].state.pos
         pos2 = Omega().bodies[keys[1]].state.pos
         ori = (pos1-pos2).normalized()
         fab += ori.outer(ori)
   fab /= numPar
   return [[fab[0],fab[1]],[fab[3],fab[4]]]
"""

def getStressAndTangent2D(scene):
   Omega().stringToScene(scene)
   st = utils.getStressAndTangent(symmetry=True)
   s = st[0]
   s = .5*(s+s.transpose())
   t=[[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]]
   t[0][0][0][0]=st[1][0,0]
   t[1][1][0][0]=t[0][0][1][1]=st[1][0,1]
   t[0][1][0][0]=t[0][0][0][1]=t[1][0][0][0]=t[0][0][1][0]=st[1][0,5]
   t[1][1][1][1]=st[1][1,1]
   t[1][1][0][1]=t[0][1][1][1]=t[1][1][1][0]=t[1][0][1][1]=st[1][1,5]
   t[0][1][0][1]=t[0][1][1][0]=t[1][0][0][1]=t[1][0][1][0]=st[1][5,5]
   return [[s[0,0],s[0,1]],[s[1,0],s[1,1]]],t

def getStressAndTangent3D(scene):
   Omega().stringToScene(scene)
   st = utils.getStressAndTangent(symmetry=True)
   s = st[0]
   s = .5*(s+s.transpose())
   t = numpy.zeros((3,3,3,3))
   t[0][0][0][0]=st[1][0]
   t[0][0][1][1]=t[1][1][0][0]=st[1][1]
   t[0][0][2][2]=t[2][2][0][0]=st[1][2]
   t[0][0][1][2]=t[0][0][2][1]=t[1][2][0][0]=t[2][1][0][0]=st[1][3]
   t[0][0][0][2]=t[0][0][2][0]=t[0][2][0][0]=t[2][0][0][0]=st[1][4]
   t[0][0][0][1]=t[0][0][1][0]=t[0][1][0][0]=t[1][0][0][0]=st[1][5]
   t[1][1][1][1]=st[1][7]
   t[1][1][2][2]=t[2][2][1][1]=st[1][8]
   t[1][1][1][2]=t[1][1][2][1]=t[1][2][1][1]=t[2][1][1][1]=st[1][9]
   t[1][1][0][2]=t[1][1][2][0]=t[0][2][1][1]=t[2][0][1][1]=st[1][10]
   t[1][1][0][1]=t[1][1][1][0]=t[0][1][1][1]=t[1][0][1][1]=st[1][11]
   t[2][2][2][2]=st[1][14]
   t[2][2][1][2]=t[2][2][2][1]=t[1][2][2][2]=t[2][1][2][2]=st[1][15]
   t[2][2][0][2]=t[2][2][2][0]=t[0][2][2][2]=t[2][0][2][2]=st[1][16]
   t[2][2][0][1]=t[2][2][1][0]=t[0][1][2][2]=t[1][0][2][2]=st[1][17]
   t[1][2][1][2]=t[1][2][2][1]=t[2][1][1][2]=t[2][1][2][1]=st[1][21]
   t[1][2][0][2]=t[1][2][2][0]=t[2][1][0][2]=t[2][1][2][0]=t[0][2][1][2]=t[2][0][1][2]=t[0][2][2][1]=t[2][0][2][1]=st[1][22]
   t[1][2][0][1]=t[1][2][1][0]=t[2][1][0][1]=t[2][1][1][0]=t[0][1][1][2]=t[1][0][1][2]=t[0][1][2][1]=t[1][0][2][1]=st[1][23]
   t[0][2][0][2]=t[2][0][0][2]=t[0][2][2][0]=t[2][0][2][0]=st[1][28]
   t[0][2][0][1]=t[0][2][1][0]=t[2][0][0][1]=t[2][0][1][0]=t[0][1][0][2]=t[1][0][0][2]=t[0][1][2][0]=t[1][0][2][0]=st[1][29]
   t[0][1][0][1]=t[0][1][1][0]=t[1][0][0][1]=t[1][0][1][0]=st[1][35]
   return [[s[0],s[1],s[2]],[s[3],s[4],s[5]],[s[6],s[7],s[8]]],t
   
def getVoidRatio2D(scene):
   Omega().stringToScene(scene)
   zSize = Omega().cell.hSize[2,2]
   return utils.voidratio2D(zlen=zSize)

def getVoidRatio3D(scene):
   Omega().stringToScene(scene)
   p = utils.porosity()
   return p/(1.0-p)

def avgRotation2D(scene):
   Omega().stringToScene(scene)
   rot = 0.0
   for b in Omega().bodies:
      rot += b.state.rot()[2]
   rot /= len(Omega().bodies)
   return rot
   
def avgRotation3D(scene):
   Omega().stringToScene(scene)
   rot = utils.Vector3.Zero
   for b in Omega().bodies:
      rot += b.state.rot()
   rot /= len(Omega().bodies)
   return [rot[0],rot[1],rot[2]]

# To record the debondding number
def DebondNum(scene):
   Omega().load('0.yade.gz')
   origBroNum = 0
   for inter in Omega().interactions:
      if inter.phys.cohesionBroken:
	     origBroNum += 1		 
   intrsOrigIds = [(i.id1,i.id2) for i in Omega().interactions]
   Omega().stringToScene(scene)
   num = 0
   for id1,id2 in intrsOrigIds:
      try:
        i = Omega().interactions[id1,id2]
        if not i.isReal:
           num +=1
           continue 
        if i.phys.cohesionBroken:
           num += 1
      except IndexError:
        num += 1
   return num - origBroNum

""" Used for perturbation method only (2D)   
def getTangentOperator(param):
   if len(param)!=2:
      raise RuntimeError,"No. of param. should be exactly 2. 0: RVE scene; 1: perturbation."
   Omega().stringToScene(param[0])
   pos = Omega().tags['id']
   perturbation = param[1]
   t=[[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]]
   ns=int(max(1e5*abs(perturbation),2))
   stress0 = utils.getStress()
   stress0 = .5*(stress0+stress0.transpose())
   #00
   strain=utils.Matrix3(perturbation,0,0, 0,0,0, 0,0,0)
   Omega().cell.velGrad=strain/(ns*Omega().dt)
   Omega().run(ns,True)
   stress1 = utils.getStress()
   stress1 = .5*(stress1+stress1.transpose())
   t[0][0][0][0]=(stress1[0]-stress0[0])/perturbation
   t[1][1][0][0]=(stress1[4]-stress0[4])/perturbation
   t[0][1][0][0]=t[1][0][0][0]=(stress1[1]-stress0[1])/perturbation
   Omega().stringToScene(param[0])
   #11
   strain=utils.Matrix3(0,0,0, 0,perturbation,0, 0,0,0)
   Omega().cell.velGrad=strain/(ns*Omega().dt)
   Omega().run(ns,True)
   stress1=utils.getStress()
   stress1=.5*(stress1+stress1.transpose())
   t[0][0][1][1]=(stress1[0]-stress0[0])/perturbation
   t[1][1][1][1]=(stress1[4]-stress0[4])/perturbation
   t[0][1][1][1]=t[1][0][1][1]=(stress1[1]-stress0[1])/perturbation
   Omega().stringToScene(param[0])
   #01
   strain=utils.Matrix3(0,perturbation,0, 0,0,0, 0,0,0)
   Omega().cell.velGrad=strain/(ns*Omega().dt)
   Omega().run(ns,True)
   stress1=utils.getStress()
   stress1=.5*(stress1+stress1.transpose())
   t[0][0][0][1]=(stress1[0]-stress0[0])/perturbation
   t[1][1][0][1]=(stress1[4]-stress0[4])/perturbation
   t[0][1][0][1]=t[1][0][0][1]=(stress1[1]-stress0[1])/perturbation
   Omega().stringToScene(param[0])
   #10
   strain=utils.Matrix3(0,0,0, perturbation,0,0, 0,0,0)
   Omega().cell.velGrad=strain/(ns*Omega().dt)
   Omega().run(ns,True)
   stress1=utils.getStress()
   stress1=.5*(stress1+stress1.transpose())
   t[0][0][1][0]=(stress1[0]-stress0[0])/perturbation
   t[1][1][1][0]=(stress1[4]-stress0[4])/perturbation
   t[0][1][1][0]=t[1][0][1][0]=(stress1[1]-stress0[1])/perturbation
   #symmetrize
   t[0][0][1][1]=t[1][1][0][0]=(t[0][0][1][1]+t[1][1][0][0])*0.5
   t[0][1][0][0]=t[1][0][0][0]=t[0][0][1][0]=t[0][0][0][1]=(t[0][1][0][0]*2+t[0][0][0][1]+t[0][0][1][0])*0.25
   t[0][1][0][1]=t[0][1][1][0]=t[1][0][0][1]=t[1][0][1][0]=(t[0][1][1][0]+t[0][1][0][1])*0.5
   t[0][1][1][1]=t[1][0][1][1]=t[1][1][0][1]=t[1][1][1][0]=(t[0][1][1][1]*2+t[1][1][0][1]+t[1][1][1][0])*0.25
   return t
"""

def finalize():
   Omega().exitNoBacktrace()

