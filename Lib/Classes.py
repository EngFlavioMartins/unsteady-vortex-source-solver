import numpy as np

class bladePanel:
      """
      Contains information of the vortexes that define a 2D thin airfoil model (2D flat plate)
      """
      def __init__(self, Xa, Xb):
            """
            Initializes the panel.

            Inputs
            ----------
            Xa:     2x1 float - Coordinates of the first end point
            Xb:     2x1 float - Coordinates of the second end point
            
            Parameters
            ----------
            Pc        float - pressure coefficient at control point
            length    float - length of the panel
            normalDir float - normal direction of the panel
            tanDir    float - tangent direction of the panel
            Ut        float - tangential velocity
            Gamma     float - source strength
            Xa        [x,y] - coordinates of the the starting point.
            Xa        [x,y] - coordinates of the the starting point.
            Xv        [x,y] - coordinates of the vortex-allocation point.
            Xcp       [x,y] - coordinates of the control point.
            """

            ## Set control point, normal and tangential directions of panel:
            self.Xa  = Xa
            self.Xb  = Xb
            self.length = np.sqrt( (Xb[0] - Xa[0])**2 + (Xb[1] - Xa[1])**2 )  # panel length
            self.Xv  = Xa + 0.25*(self.Xb - self.Xa)  # 1/4 chord length
            self.Xcp = Xa + 0.75*(self.Xb - self.Xa)  # 3/4 chord length

            self.un = 0.0  # x-velocity over time
            self.ut = 0.0  # y-velocity over time

            ## Orientation of panel (angle between x-axis and panel's normal)
            if Xb[0] - Xa[0] <= 0.0:
                angle = np.arccos( (Xb[1] - Xa[1]) / self.length)
            elif Xb[0] - Xa[0] > 0.0:
                angle = np.pi + np.arccos( -(Xb[1] - Xa[1]) / self.length)

            self.n = np.array([np.cos(angle), -np.sin(angle)]) ## tangential kinetic velocity of the panel
            self.t = np.array([np.sin(angle), np.cos(angle)])  ## normal kinetic velocity of the panel

            self.Gamma = 1.0  # source strength
            self.Lift  = 0.0  # lift force, N/m
            self.Pc    = 0.0  # pressure coefficient


      def rotate(self, origin, theta):
            """Rotate the panel by an angle theta about origin"""
            c, s = np.cos(theta), np.sin(theta)
            j = np.array([[c, s], [-s, c]])

            self.Xa = origin + np.dot(j, self.Xa-origin)  
            self.Xb = origin + np.dot(j, self.Xb-origin)  
            self.Xv = origin + np.dot(j, self.Xv-origin)
            self.Xcp = origin + np.dot(j, self.Xcp-origin) 

            ## Orientation of panel (angle between x-axis and panel's normal)
            if self.Xb[0] - self.Xa[0] <= 0.0:
                angle = np.arccos( (self.Xb[1] - self.Xa[1]) / self.length)
            elif self.Xb[0] - self.Xa[0] > 0.0:
                angle = np.pi + np.arccos( -(self.Xb[1] - self.Xa[1]) / self.length)

            self.n = np.array([np.cos(angle), -np.sin(angle)]) ## tangential kinetic velocity of the panel
            self.t = np.array([np.sin(angle), np.cos(angle)])  ## normal kinetic velocity of the panel            

            
      def translate(self, dx, dy):
            """Translate the points in a panel by distances (dx,dy)"""
            self.Xa  = self.Xa  + np.array([dx,dy]) 
            self.Xb  = self.Xb  + np.array([dx,dy])
            self.Xv  = self.Xv  + np.array([dx,dy]) 
            self.Xcp = self.Xcp + np.array([dx,dy]) 



class wakeVortex:
      """
      Contains information of the vortexes that define the 2D wake behind the airfoil model (2D flat plate)
      """
      def __init__(self, coord):
            """
            Initializes the panel.

            Inputs
            ----------
            coord: 2x1 float - Coordinates of the first end point
            
            Parameters
            ----------
            Xcp       [x,y] - coordinates of the control point.
            """

            ## Set control point, normal and tangential directions of panel:
            self.Xv = np.array(coord)
            self.Gamma = 1  # source strength
