from yade import pack
#this script is used to create a standard RVE '0.yade.gz',and should be executed as 'yade prepareRVE.py'
O.materials.append(FrictMat(young=6.e8,poisson=.8,frictionAngle=.0))

sp = pack.SpherePack()
size = .24
sp.makeCloud(minCorner=(0,0,.05),maxCorner=(size,size,.05),rMean=.005,rRelFuzz=.4,num=400,periodic=True,seed=1)
#makeCloud: create a pack of particles. min and max has the same z=0.05,meaning this is a 2D packing.
sp.toSimulation()
O.cell.hSize = Matrix3(size,0,0, 0,size,0, 0,0,.1)
#cell.hSize: Base cell vectors (columns of the matrix) cell is used because of periodic condition
#hSize need to be updated at every step from velGrad
print len(O.bodies)
for p in O.bodies:
   p.state.blockedDOFs = 'zXY'
   #blockedDOFs: block the degrees of freedom of the particle. String that may contain 'xyzXYZ' (translations and rotations).
   p.state.mass = 2650 * 0.1 * pi * p.shape.radius**2 # 0.1 = thickness of cylindrical particle
   inertia = 0.5 * p.state.mass * p.shape.radius**2
   p.state.inertia = (.5*inertia,.5*inertia,inertia)

O.dt = utils.PWaveTimeStep()
#Get timestep accoring to the velocity of P-Wave propagation; computed from sphere radii, rigidities and masses

O.engines = [
   ForceResetter(),
   InsertionSortCollider([Bo1_Sphere_Aabb()]),
   InteractionLoop(
      [Ig2_Sphere_Sphere_ScGeom()],  #collision geometry
      [Ip2_FrictMat_FrictMat_FrictPhys()],   #collision physics
      [Law2_ScGeom_FrictPhys_CundallStrack()]    #contact law
   ),
   PeriTriaxController(
      dynCell=True,
      goal=(-1.e5,-1.e5,0),
      stressMask=3,
      relStressTol=.001,
      maxUnbalanced=.001,
      maxStrainRate=(.5,.5,.0),
      doneHook='term()', #python command to be run when the desired state is reached 
      label='biax'
   ),
   #controller is used to define desired state
   NewtonIntegrator(damping=.1)  # Apply gravity force to particles(here the gravity is neglected). damping: numerical dissipation of energy.
]

def term():
   O.engines = O.engines[:3]+O.engines[4:]
   print getStress()
   print O.cell.hSize
   #print avgNumInteractions()
   setContactFriction(0.5) #Modify the friction angle (in radians) inside the material classes and existing contacts. 
   O.cell.trsf=Matrix3.Identity
   O.cell.velGrad=Matrix3.Zero
   for p in O.bodies:
      p.state.vel = Vector3.Zero
      p.state.angVel = Vector3.Zero
      p.state.refPos = p.state.pos
      p.state.refOri = p.state.ori
   O.save('0.yade.gz')
   O.pause()

O.run();O.wait()

