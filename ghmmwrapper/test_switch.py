import ghmm
import ghmmwrapper
import ghmmhelper

A2 = [
      [[0.0,1.0,0.0],[0.0,0.0,1.0],[1.0,0.0,0.0]],
      [[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0]]
     ]
B2 = [[0.0,0.000001],[1.0,0.000001], [2.0,0.000001]]
pi2 = [1.0,0.0,0.0]
model2 = ghmm.HMMFromMatrices(ghmm.Float(),ghmm.GaussianDistribution(ghmm.Float), A2, B2, pi2)

ghmmwrapper.smodel_class_change_alloc(model2.cmodel)
ghmmwrapper.setPythonSwitching(model2.cmodel,"class_change","getClass")
#ghmmwrapper.setSwitchingFunction(model2.cmodel)

#print model2
seq = model2.sample(2,30)
