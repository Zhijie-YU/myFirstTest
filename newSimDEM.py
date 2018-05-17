
# import yade modules
import sys
sys.path.append('/home/littlefish/yade/install/bin') # path where you have yadeimport.py!!! this should be changed under different conditions
from yadeimport import *



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

def getStressAndTangent2D(scene):
   Omega().stringToScene(scene)
   st = utils.getStressAndTangent(symmetry=True)  #getStressandTangent is from Yade.
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

# I add this function by myself to output stress tensor(4 values) and tangent operator(6 values)
def getStressAndTangent2DModi(scene):
   Omega().stringToScene(scene)
   st = utils.getStressAndTangent(symmetry=True)  #getStressandTangent is from Yade.
   s = st[0]
   s = .5*(s+s.transpose())
   t=[[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]]
   t[0][0][0][0]=st[1][0,0]
   t[1][1][0][0]=t[0][0][1][1]=st[1][0,1]
   t[0][1][0][0]=t[0][0][0][1]=t[1][0][0][0]=t[0][0][1][0]=st[1][0,5]
   t[1][1][1][1]=st[1][1,1]
   t[1][1][0][1]=t[0][1][1][1]=t[1][1][1][0]=t[1][0][1][1]=st[1][1,5]
   t[0][1][0][1]=t[0][1][1][0]=t[1][0][0][1]=t[1][0][1][0]=st[1][5,5]
   return [[s[0,0],s[0,1]],[s[1,0],s[1,1]]],t,[s[0,0],s[0,1],s[1,0],s[1,1]],[st[1][0,0],st[1][0,1],st[1][0,5],st[1][1,1],st[1][1,5],st[1][5,5]]

def getFabric2D(scene):
   Omega().stringToScene(scene)
   f = utils.fabricTensor(splitTensor=False)[0]
   return [[f[0,0],f[0,1]],[f[1,0],f[1,1]]]

def getVoidRatio2D(scene):
   Omega().stringToScene(scene)
   zSize = Omega().cell.hSize[2,2]
   return utils.voidratio2D(zlen=zSize)

def getStressTangentFabricAndVoid(out):
   #print type(out) ; out.shape==(1,15)
   stress = [[out[0][0],out[0][1]],[out[0][2],out[0][3]]]
   t = [[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]]
   t[0][0][0][0]=out[0][4]
   t[1][1][0][0]=t[0][0][1][1]=out[0][5]
   t[0][1][0][0]=t[0][0][0][1]=t[1][0][0][0]=t[0][0][1][0]=out[0][6]
   t[1][1][1][1]=out[0][7]
   t[1][1][0][1]=t[0][1][1][1]=t[1][1][1][0]=t[1][0][1][1]=out[0][8]
   t[0][1][0][1]=t[0][1][1][0]=t[1][0][0][1]=t[1][0][1][0]=out[0][9]  
   fabric = [[out[0][10],out[0][11]],[out[0][12],out[0][13]]]
   void = out[0][14]

   return [stress, t, fabric, void]


