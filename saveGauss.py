"""
Function to output Gauss point tensorial values.
"""
def saveGauss2D(name='',pos=(),**kwargs):
   #**kwargs refers to a dictionary
   fout = file(name,'w')
   #key is the name of each dict item
   for key in kwargs: 
      data = kwargs[key].toListOfTuples()
      if len(pos)==0:
         fout.write('%s '%key+str(len(data))+'\n')
         for i in xrange(len(data)):
            fout.write(' '.join('%s %s'%x for x in data[i])+'\n')
      else:
         fout.write('%s '%key+str(len(pos))+'\n')
         for i in pos:
            fout.write(' '.join('%s %s'%x for x in data[i])+'\n')
   fout.close()
   
def saveGauss3D(name='',pos=(),**kwargs):
   fout = file(name,'w')
   for key in kwargs:
      data = kwargs[key].toListOfTuples()
      if len(pos)==0:
         fout.write('%s '%key+str(len(data))+'\n')
         for i in xrange(len(data)):
            fout.write(' '.join('%s %s %s'%x for x in data[i])+'\n')
      else:
         fout.write('%s '%key+str(len(pos))+'\n')
         for i in pos:
            fout.write(' '.join('%s %s %s'%x for x in data[i])+'\n')
   fout.close()
   
