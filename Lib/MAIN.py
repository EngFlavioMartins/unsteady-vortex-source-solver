import numpy as np
import pickle
import os
import shutil

from Lib.Classes       import wakeVortex

from Lib.Subsolvers import updateInductionMatrix
from Lib.Subsolvers import moveBlade
from Lib.Subsolvers import wakeInducedVelocity
from Lib.Subsolvers import vortexRollUp
from Lib.Subsolvers import initializeInductionMatrix




def solveLinearSystem(A, WIV, myBladePanels, myTEVortex):

      npanels = len(myBladePanels)

      ## Construct the RHS:
      RHS = np.zeros(npanels+1)

      for p, iPanel in enumerate(myBladePanels):
            ## Normal kinetic velocity of the panel + wake-induced velocity
            
            RHS[p] =  -(iPanel.un + WIV[p]) 
            RHS[-1] += iPanel.Gamma 

      newGamma = np.linalg.solve(A, RHS)

      for p, panel in enumerate(myBladePanels):
            panel.Gamma = newGamma[p]

      myTEVortex.Gamma = newGamma[-1] ## update the circulation of the trailing edge vortex





def solveVortexPanels(timeStep, U, V, myBladePanels, Theta, myTEVortex, myWakeVortexes, saveToFile=True):



      ## Create a directory to store solution data:
      if saveToFile:
            myPath = './Plots/Animate/'
            if os.path.exists(myPath):
                  try:
                        shutil.rmtree(myPath)
                  except:
                        print('Folder not deleted')
            os.mkdir(myPath)
    

      ## (Induced velocity from blade panels onto blade panels)
      A = initializeInductionMatrix(myBladePanels, myTEVortex) # (npanels+1) x (npanels+1) - size

      solution = np.zeros(len(timeStep)) ## total circulation over the airfoil 

      #Uinf = np.array([-max(-U), -max(V)]) # Normalize by maximum velocity

      for t, dt in enumerate(timeStep):

            # Move the blade and the TE vortex 1 step forward in time
            moveBlade(t, dt, U, V, myBladePanels, Theta, myTEVortex)

            updateInductionMatrix(myBladePanels, myTEVortex, A)

            WIV = wakeInducedVelocity(myBladePanels, myWakeVortexes)
            
            solveLinearSystem(A, WIV, myBladePanels, myTEVortex)
            
            vortexRollUp(myBladePanels, myTEVortex, myWakeVortexes, dt)

            #logGamma(myBladePanels, myTEVortex, myWakeVortexes)

            myWakeVortexes[t] = wakeVortex(myTEVortex.Xv)
            myWakeVortexes[t].Gamma = myTEVortex.Gamma

            iterable = (panel.Gamma for panel in myBladePanels)
            solution[t] = np.sum( np.fromiter(iterable, float) )

            Uinf = np.array([-U[t], -V[t]]) ## define the reference velocity

            if saveToFile and t>0:
                  data = {'myBladePanels':myBladePanels, 'myWakeVortexes':myWakeVortexes, 'myTEVortex':myTEVortex, 'Uinf':Uinf, 'timeStep': dt}

                  with open('./Plots/animate/backup_' + str(t), 'wb') as outputFile:
                        pickle.dump(data, outputFile)

      return solution

