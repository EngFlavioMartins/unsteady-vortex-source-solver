import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import pickle
import glob
from PIL import Image
import os
from matplotlib import colors
from Lib.Subsolvers import inducedVelocity

cm2in = 1/2.54


def showSystem(myBladePanels, myTEVortex, myWakeVortexes, Uinf, timeStep):
    
      fig, ax = plt.subplots(figsize=(15*cm2in, 5.0*cm2in), dpi=150)
      ax.set_xlim([-0.1, 3.0])
      ax.set_ylim([-0.5, 0.5])
      ax.set_xlabel(r'$x/c$')
      ax.set_ylabel(r'$y/c$')
      ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
      ax.grid(which='minor', color='#EEEEEE', linestyle=':', linewidth=0.5)
      ax.minorticks_on()
      #ax.set_title(r'Wake circulation strength, $\Gamma/U_\infty^2 \Delta t$')
      ax.set_aspect('1')

      LE = myBladePanels[0].Xa ## find leading edge coordinates
      TE = myBladePanels[-1].Xb-[LE[0],0] ## find the trailing edge cooedinates
      Umag2 = Uinf[0]**2 + Uinf[1]

      ## Lets color the shed vortexes using the strength (Gamma)
      wakeVortexColor = np.zeros(len(myWakeVortexes))
      wakeVortexX    = np.zeros_like(wakeVortexColor)
      wakeVortexY    = np.zeros_like(wakeVortexColor)
      for j, vortex in enumerate(myWakeVortexes):
            if vortex != None:
                  wakeVortexColor[j] = vortex.Gamma
                  wakeVortexX[j]     = vortex.Xv[0]-LE[0]
                  wakeVortexY[j]     = vortex.Xv[1]
      wakeVortexColorNorm            = wakeVortexColor / (Umag2 * timeStep)


      ## Lets color the airfoil vortexes using the strength (Gamma)
      airfoilVortexColor = np.zeros(len(myBladePanels))
      airfoilVortexX     = np.zeros_like(airfoilVortexColor)
      airfoilVortexY     = np.zeros_like(airfoilVortexColor)
      for j, panel in enumerate(myBladePanels):
            airfoilVortexColor[j] = panel.Gamma
            airfoilVortexX[j]     = panel.Xv[0]-LE[0]
            airfoilVortexY[j]     = panel.Xv[1]+0.1
      airfoilVortexColorNorm = airfoilVortexColor / (Umag2 * timeStep)


      ## Show the trailing edge vortex:
      """ 
      vmin = min(airfoilVortexColorNorm)
      vmax = max(airfoilVortexColorNorm)
      norm = colors.Normalize(vmin=vmin, vmax=vmax)
      AirfoilVortexPlot = ax.scatter(airfoilVortexX, airfoilVortexY, marker='^',
                                     c=airfoilVortexColorNorm, s=3, cmap='plasma', norm=norm)
      """
      
      ## Show the wake vortexes:
      vmin = min(wakeVortexColorNorm)
      vmax = max(wakeVortexColorNorm)
      norm = colors.Normalize(vmin=vmin, vmax=vmax)
      WakeVortexPlot = ax.scatter(wakeVortexX, wakeVortexY, 
                                  c=wakeVortexColorNorm, s=1, cmap='plasma', norm=norm)

      
      AirfoilPlotA = ax.plot([0, TE[0]], [LE[1], TE[1]], '-', color='black', linewidth=2.6)
      AirfoilPlotB = ax.plot([0, TE[0]], [LE[1], TE[1]], '-', color='slategrey', linewidth=2.0)

      cb1 = fig.colorbar(WakeVortexPlot, shrink=0.7, pad=0.05)
      #cb2 = fig.colorbar(AirfoilVortexPlot, shrink=0.7, pad=0.05)

      cb1.set_label(r'$\Gamma \, /\, U_\infty^2 \Delta t$')
      #cb2.set_label(r'Airfoil')
      

      plt.subplots_adjust(left=0.11, bottom=0.11, right=1.0, top=1.0)





def plotFields(myWakeVortexes, myBladePanels, myTEVortex, Uinf):


      LE = myBladePanels[0].Xa     ## find leading edge coordinates
      TE = myBladePanels[-1].Xb-LE ## find the trailing edge cooedinates

      nx, ny = [94, 34]
      x = np.linspace(LE[0]-1.2, LE[0]+8.2, num=nx)
      y = np.linspace(LE[1]-1.3, LE[1]+1.3, num=ny)
      xg, yg = np.meshgrid(x,y)
      ug, vg = Uinf[0]*np.ones((ny,nx)), Uinf[1]*np.zeros((ny,nx))   # (as before: no cross-stream free stream)

      ## Velocity induced by every point vortex (wake, blade and trailing edge) at every grid point
      vortexes = [v for v in myWakeVortexes if v is not None] + list(myBladePanels) + [myTEVortex]
      vx = np.array([v.Xv for v in vortexes], float)
      vgam = np.array([v.Gamma for v in vortexes], float)
      dx = xg[..., None] - vx[:, 0]
      dy = yg[..., None] - vx[:, 1]
      c = vgam / (2 * np.pi * (dx**2 + dy**2))
      ug += np.sum(c * dy, axis=-1)
      vg += np.sum(-c * dx, axis=-1)

      ## Get the fields:
      Umag  = np.sqrt(Uinf[0]**2 + Uinf[1]**2)  ## reference velocity
      uNorm = gaussian_filter(ug/Umag, sigma=0.8, order=0)
      vNorm = gaussian_filter(vg/Umag, sigma=0.8, order=0)
      dvdx  = np.diff(vNorm,axis=1)/(x[1]-x[0])
      dudy  = np.diff(uNorm,axis=0)/(y[1]-y[0])
      Omega = dvdx[1:,:] - dudy[:,1:]  ## vorticity field
      pc = 1.0 - (uNorm**2 + vNorm**2) ## Pressure coefficient:


      ## Move everything to LE=(0,0):
      xg = xg-LE[0]
      yg = yg-LE[1]


      fig, axs = plt.subplots(2, 1, figsize=(15*cm2in, 9*cm2in), dpi=150)

      plt.subplots_adjust(left=0.15, bottom=0.10, right=1.0, top=0.95)


      omgCont = axs[0].contourf(xg[1:,1:], yg[1:,1:], Omega/Umag, cmap = 'rainbow', levels = 15, extend='both')
      axs[0].set_title(r'Vorticity field, $\omega_z c/U_\infty$')


      pCont = axs[1].contourf(xg, yg, pc, cmap='Spectral_r', levels = 15, extend='both')
      axs[1].set_title(r'Pressure coefficient, $C_p$')


      for ax in axs:
            ax.streamplot(xg, yg, ug, vg, linewidth=0.5, color='#EEEEEE', arrowsize=1)
            ax.set_xlim([-1.0, 8.0])
            ax.set_ylim([-1.2, 1.2])
            ax.set_ylabel(r'$y/c$')
            ax.set_aspect('1')  
            ax.plot([0, TE[0]], [0, TE[1]], '-', color='black', linewidth=2.6)
            ax.plot([0, TE[0]], [0, TE[1]], '-', color='slategrey', linewidth=2.0)   
      ax.set_xlabel(r'$x/c$')
      
      fig.colorbar(omgCont, shrink=0.8, pad=0.05,format=lambda x, _: f"{x:.1f}")
      fig.colorbar(pCont, shrink=0.8, pad=0.05, format=lambda x, _: f"{x:0.1f}")




def makeGif(frame_folder):
    """ Makes a gif out of the animation files stored in the system"""
    frames = [Image.open(image) for image in  sorted(glob.glob(f"{frame_folder}/*.png"), key=os.path.getmtime)]
    frame_one = frames[0]
    frame_one.save("Animation.gif", format="GIF", append_images=frames,
               save_all=True, duration=70, loop=0)




def makeAnimation(nframes):
      """ creates animations of the vortex sheeding """

      for frame in range(1,nframes-2):

            file = open('./Plots/Animate/backup_' + str(frame),'rb')
            backupData = pickle.load(file)
            file.close()

            myBladePanels = backupData.get("myBladePanels")
            myTEVortex = backupData.get("myTEVortex")
            myWakeVortexes = backupData.get("myWakeVortexes")
            Uinf = backupData.get("Uinf")
            timeStep = backupData.get("timeStep")

            showSystem(myBladePanels, myTEVortex, myWakeVortexes, Uinf, timeStep)

            plt.savefig("./Plots/Animate/ %03d.png" % (frame-1), dpi=300)

            plt.close()
      

      

      makeGif('./Plots/Animate/')