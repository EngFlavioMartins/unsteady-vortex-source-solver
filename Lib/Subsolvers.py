import numpy as np
from numba import njit

from Lib.Classes       import bladePanel
from Lib.Classes       import wakeVortex


def inducedVelocity(Point, jPanel, Gamma=1):
      """ Returns the  2D induced velocity from j-th onto i-th panel. The code is hard to read, but it is more efficient than declaring variables"""

      return Gamma / (2 * np.pi * ((Point[0] - jPanel.Xv[0])**2 + (Point[1] - jPanel.Xv[1])**2)) * np.array([[0,1],[-1,0]]).dot(np.array([(Point[0] - jPanel.Xv[0]),(Point[1] - jPanel.Xv[1])]))


def initializeInductionMatrix(myBladePanels, myTEVortex):
      """
      This is a (npanels+1)x(npaels+1) matrix that works as a LHS for the linear system
      A[i,j] = velocity normal to the i-th panel due to the j-th panel vortex
      """
      npanels = len(myBladePanels)
      A = np.ones([npanels+1,npanels+1])

      for i, iPanel in enumerate(myBladePanels):
            for j, jPanel in enumerate(myBladePanels):
                  A[i,j] = inducedVelocity(iPanel.Xcp, jPanel).dot(iPanel.n) # notice that Panel.Gamma = 1 at this stage
            A[i,-1] = inducedVelocity(iPanel.Xcp, myTEVortex).dot(iPanel.n)


      # After all operations are done, zero the circulations at all panels:
      myTEVortex.Gamma = 0

      for panel in myBladePanels: panel.Gamma=0

      return A


def updateInductionMatrix(myBladePanels, myTEVortex, A):
      """
      Updated the last column of the induction matrix
      """

      for i, iPanel in enumerate(myBladePanels):
            A[i,-1] = inducedVelocity(iPanel.Xcp, myTEVortex).dot(iPanel.n)



def moveBlade(t, dt, U, V, myBladePanels, Theta, myTEVortex):
      """ In this subroutine, we move the airfoil and advect the wake behind it"""

      TE = myBladePanels[-1].Xb  ## Current coordinate of the trailing edge
      Origin = myBladePanels[0].Xa

      ## Tranlate panels in space:
      for panel in myBladePanels:

            XcpOld = panel.Xcp                              ## get current location of control point

            panel.translate( U[t] * dt,  V[t] * dt)         ## Translate blade panels
            panel.rotate(Origin, Theta[t+1] - Theta[t])     ## Rotate points about TE coords.

            u = (panel.Xcp[0] - XcpOld[0])/dt               ## Update kinematic velocity of panels
            v = (panel.Xcp[1] - XcpOld[1])/dt

            panel.un = np.array([u,v]).dot(panel.n)         ## kinematic velocity in the normal and tan. dirs.
            panel.ut = -np.array([u,v]).dot(panel.t) 

      ## Update position of the trailing edge vortex to 20% of the distance travelled by the TE
      myTEVortex.Xv = TE - 0.3*(panel.Xb - TE) 



def _active_wake(myWakeVortexes):
      """The wake vortexes created so far (in shedding order)."""
      return [v for v in myWakeVortexes if v is not None]


@njit(cache=True)
def _normal_wake_velocity(cps, normals, wx, wg):
      """Velocity normal to each panel induced by point vortexes at wx with circulations wg."""
      n, m = cps.shape[0], wx.shape[0]
      out = np.zeros(n)
      for i in range(n):
            acc = 0.0
            for j in range(m):
                  dx = cps[i, 0] - wx[j, 0]
                  dy = cps[i, 1] - wx[j, 1]
                  c = wg[j] / (2 * np.pi * (dx * dx + dy * dy))
                  acc += c * dy * normals[i, 0] - c * dx * normals[i, 1]
            out[i] = acc
      return out


def wakeInducedVelocity(myBladePanels, myWakeVortexes):
      """
      This is the velocity the wake induced on the panels. W[i] = velocity normal to the i-th panel
      """
      wake = _active_wake(myWakeVortexes)
      if not wake:
            return np.zeros(len(myBladePanels))
      cps = np.array([p.Xcp for p in myBladePanels], float)
      normals = np.array([p.n for p in myBladePanels], float)
      wx = np.array([v.Xv for v in wake], float)
      wg = np.array([v.Gamma for v in wake], float)
      return _normal_wake_velocity(cps, normals, wx, wg)




def LEmotion(nTimeSteps, finalTime, xVelocity, yVelocity, rotSpeed, rotAmplitude, theta0, plungAmpl, plungFreq):

      """Simple function that defines the motion of the LE as constant in the X-direction"""

      rotSpeedRad = rotSpeed * 2 * np.pi
      plungFreqRad = plungFreq * 2 * np.pi
      
      Time = np.linspace(0, finalTime, num=nTimeSteps)
      dTime = Time[1:] - Time[0:-1]

      U = np.zeros_like(Time)
      V = np.zeros_like(Time)
      X = np.zeros_like(Time)
      Y = np.zeros_like(Time)
      Theta = np.zeros_like(Time)
      Omega = np.zeros_like(Time)

      for t, thisTime in enumerate(Time):
            U[t] = xVelocity       ## constant speed, m/s
            V[t] = yVelocity + plungAmpl * plungFreqRad * np.sin(plungFreqRad * thisTime) # m/s
            Omega[t] = rotAmplitude * rotSpeedRad * np.cos(rotSpeedRad * thisTime) ## rotational speed, rad/s

            X[t] = U[t] * thisTime ## position of the LE at t=thisTime, m
            Y[t] = V[t] * thisTime ## position of the LE at t=thisTime, m
            Theta[t] = rotAmplitude * np.sin(rotSpeedRad * thisTime) + theta0   ## angle of the thin airfoil, rad
            

      return Time, dTime, X, Y, U, V, Theta, Omega


def logGamma(myBladePanels, myTEVortex, myWakeVortexes):
      """ Compute total gamma present in the blades and wake vortexes """
      totalGammaBlade = 0
      for blade in myBladePanels: totalGammaBlade += blade.Gamma
      
      totalGammaWake = 0
      for vortex in myWakeVortexes:
            if vortex != None:
                  totalGammaWake += vortex.Gamma

      totalGamma = myTEVortex.Gamma + totalGammaWake + totalGammaBlade

      print('Total Gamma (Blade): %0.2f' % totalGammaBlade )
      print('Total Gamma (TE): %0.2f' % myTEVortex.Gamma )
      print('Total Gamma (Wake): %0.2f' % totalGammaWake )
      print('Total Gamma: %0.2f\n' % totalGamma )


@njit(cache=True)
def _roll_up(wx, wg, bx, bg, tex, teg, dt):
      """Move every wake vortex with the velocity induced by the trailing-edge vortex, the blade
      vortexes and the other wake vortexes. Vortexes are moved one after another in place, so
      each one sees the new positions of those moved before it, as in the original loop."""
      n = wx.shape[0]
      for i in range(n):
            dx = wx[i, 0] - tex[0]
            dy = wx[i, 1] - tex[1]
            c = teg / (2 * np.pi * (dx * dx + dy * dy))
            ux = c * dy
            uy = -c * dx
            for b in range(bx.shape[0]):
                  dx = wx[i, 0] - bx[b, 0]
                  dy = wx[i, 1] - bx[b, 1]
                  c = bg[b] / (2 * np.pi * (dx * dx + dy * dy))
                  ux += c * dy
                  uy += -c * dx
            for j in range(n):
                  if j != i:
                        dx = wx[i, 0] - wx[j, 0]
                        dy = wx[i, 1] - wx[j, 1]
                        c = wg[j] / (2 * np.pi * (dx * dx + dy * dy))
                        ux += c * dy
                        uy += -c * dx
            wx[i, 0] += ux * dt
            wx[i, 1] += uy * dt


def vortexRollUp(myBladePanels, myTEVortex, myWakeVortexes, dt):
      """
      Determine the new positions of the wake vortexes:
      """
      wake = _active_wake(myWakeVortexes)
      if not wake:
            return
      wx = np.array([v.Xv for v in wake], float)
      wg = np.array([v.Gamma for v in wake], float)
      bx = np.array([b.Xv for b in myBladePanels], float)
      bg = np.array([b.Gamma for b in myBladePanels], float)
      _roll_up(wx, wg, bx, bg, np.asarray(myTEVortex.Xv, float), float(myTEVortex.Gamma), dt)
      for v, x in zip(wake, wx):
            v.Xv[:] = x



def constructBladePanels(theta0, Npanels, nTimeSteps):
      """ Construct the blade panels given an AoA and number of panels. Notice that the chord length, c, is always 1 """
      
      x = np.linspace(0,1,num=Npanels+1)
      y = np.zeros_like(x)
      Rzz = np.array([[np.cos(-theta0), -np.sin(-theta0)],
                      [np.sin(-theta0),  np.cos(-theta0)]]) # rotation matrix


      ## Populate bladePanels list:
      bladePanels = np.empty(Npanels, dtype=object)

      for p in range(Npanels):
            
            Xa = Rzz.dot( np.array([ x[p],   y[p]   ]) )
            Xb = Rzz.dot( np.array([ x[p+1], y[p+1] ]) )

            bladePanels[p] = bladePanel(Xa, Xb)

      TEcoord = [1.2*np.cos(theta0), -1.2*np.sin(theta0)]
     
      TEVortex = wakeVortex(TEcoord)

      myWakeVortexes = np.empty(nTimeSteps-1, dtype=object) # constructWakeVortexes(nTimeSteps)
             
      return bladePanels, TEVortex, myWakeVortexes