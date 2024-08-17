import numpy as np
import matplotlib.pyplot as plt
#import torch 

class Beam:
    def __init__(self,type = 'gaussian',N = 2**10):
        # Physics parametrs
        self.wavelength = 633e-9                    # In meters 
        self.waist = 1e-3  
        L = 2*self.waist

        # Auxiliar Vector
        Aux = np.arange(-N/2 , N/2 + 1) * 2/N

        # Real space
        self.X, self.Y = np.meshgrid(L*Aux , L*Aux)
        self.R = np.sqrt(self.X**2 + self.Y**2)
        self.theta = np.arctan2(self.Y,self.X)

        # Frequency space
        k = 2*np.pi/self.wavelength
        zr = k*(self.waist**2)/2
        delta = 2*L/N
        kmax = np.pi/delta
        self.kx , self.ky = np.meshgrid(kmax*Aux,kmax*Aux)

        # Campos Iniciales
        if type =='gaussian':
            self.U0 = np.exp( -(self.R/self.waist)**2 )
        
        self.I0 = np.abs(self.U0)
        
    # Function to plot the initial Intensity
    def show_intensity0(self):
        fig, ax = plt.subplots(1,1)
        cax = ax.imshow(self.I0,cmap='hot')
        cbar = fig.colorbar(cax)
        plt.show()
        return fig


u0 = Beam('gaussian',N=2**10)
#u0.show_intensity0()