import numpy as np
import matplotlib.pyplot as plt
#import torch 

class Beam:
    def __init__(self,type = 'gaussian',N = 2**8):
        # Physics parametrs
        self.wavelength = 633e-9                    # In meters 
        self.waist = 1e-3  
        L = 2*self.waist

        # Auxiliar Vector
        if N % 2 == 0:
            Aux = np.arange(-N/2 , N/2) * 2/N
        else:
            Aux = np.linspace(-N/2 , N/2,N) * 2/N
        self.N = N                                  # Guardamos el tamaño

        # Real space
        self.x = L*Aux
        self.y = L*Aux
        self.X, self.Y = np.meshgrid(L*Aux , L*Aux)
        self.R = np.sqrt(self.X**2 + self.Y**2)
        self.theta = np.arctan2(self.Y,self.X)

        # Frequency space
        k = 2*np.pi/self.wavelength
        delta = 2*L/N
        kmax = np.pi/delta
        self.kx , self.ky = np.meshgrid(kmax*Aux,kmax*Aux)
        self.k0 = 2*np.pi/self.wavelength
        
        # Z-axis vector and some cases
        zr = k*(self.waist**2)/2
        z = np.linspace(0,zr,N)                 # Caso 1
        #z = Aux*zr*2                             # Caso 2
        self.z = z                              # Guardamos vector de 'z'



        # Select the initial field
        # More fields will be added in a near future
        # You can import you own field in this section
        if type =='gaussian':
            self.U0 = np.exp( -(self.R/self.waist)**2 )
        
        # Put the initial intensity as a property
        self.I0 = np.abs(self.U0)
        
    def paraxial_propagation(self):
        uz = np.zeros((self.N ,self.N,self.N)) + 0j
        #Fourier Transform of initial field 
        u0_fft_fftshift = np.fft.fftshift(np.fft.fft2(self.U0))
        for zi in range(self.N):
            kz = np.sqrt(self.k0**2 - (self.ky**2 + self.kx**2) + 0j)
            uz_fft = u0_fft_fftshift * np.exp(1j * self.z[zi] * kz)
            uz[:,:,zi] = np.fft.ifft2(np.fft.ifftshift(uz_fft))
            
        self.Uz =uz
        self.Iz = np.abs(uz)**2

        

        
    # Function to plot the initial Intensity
    def show_intensity0(self):
        plt.imshow(self.I0,cmap='gray')
        plt.colorbar()
        plt.show()

u0 = Beam('gaussian',N=2**8)
u0.show_intensity0()
u0.paraxial_propagation()
ax = plt.imshow(u0.Iz[:,128,:],cmap='gray')
plt.show()



