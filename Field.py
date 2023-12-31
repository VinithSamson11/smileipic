from smilei_project.extractdata import *
#from smilei_project.sim_params import *
import numpy as np

class Field:
    
    def __init__(self, h5item, field, moving, timesteps, modes,theta):
        
        self.h5item = h5item
        self.field = field
        self.moving = moving
        self.user_timesteps = timesteps
        self.theta = theta
        self.modes = modes
        
        self.window_time = (MovingWindow.time_start)/fs
        
        #Finding the simulation_timestep[sim_timesteps] and real time[real_time]
        timestep_list = [int(x) for x in list(h5item['data'].keys())]
        self.int_timesteps = min(timestep_list, key=lambda x: (abs(x - self.user_timesteps), x))
        self.sim_timesteps = str(self.int_timesteps).zfill(10)   #timestep format in the data dump
        self.real_time=int(self.sim_timesteps)*dt/fs    #timestep into femto-seconds unit
        
        self.AM_modes = Main.number_of_AM
        #self.field_name,self.field_array,self.real_field,self.imag_field=[],[],[],[]
        array_list=[]
        for mode in range(self.AM_modes):
        
            #self.field_name.append(self.field[1:]+"_mode_"+str(mode))
            #self.field_array.append(self.h5item['data'][self.sim_timesteps][self.field_name[mode]]*int(self.field[0]+"1"))
            #self.real_field.append(self.field_array[mode][:,::2])
            #self.imag_field.append(self.field_array[mode][:,1::2])

            field_name = self.field[1:]+"_mode_"+str(mode)
            field_array = self.h5item['data'][self.sim_timesteps][field_name]['data']*int(self.field[0]+"1")
            field_array_shape = field_array.shape
            reshaped_array = field_array.reshape(field_array_shape[0], field_array_shape[1]//2, 2)
            new_array = reshaped_array[:, :, 0] + 1j * reshaped_array[:, :, 1]
            array_list.append(new_array)
        
        self.mod_data= np.stack(array_list,axis=0)       #Modified array of shape (no.of.modes, Nx, Nr)
     
    def getNmodes(self):
        strings=list(self.h5item['data'][self.sim_timesteps].keys())
        max_suffix = float('-inf')
        max_suffix_string = None
    
        for s in strings:
            prefix, suffix = s.split('_mode_')
            try:
                suffix_int = int(suffix)
            except ValueError:
                continue
            if suffix_int > max_suffix:
                max_suffix = suffix_int
                max_suffix_string = s
    
        return int(max_suffix_string[-1])
    
    def mode_expansion_DFT(self):
        
        F = np.zeros_like(np.real(self.mod_data[0]))                       
        '''for m in self.modes:
            F += self.mod_data[m]*np.exp(-1j*m*self.theta)
        '''
        for m in self.modes:
            F += np.real(self.mod_data[m])*np.cos(m*self.theta)+np.imag(self.mod_data[m])*np.sin(m*self.theta)

        return F   
    
    def mode_expansion_FFT(self):

        rawdata=self.mod_data
        Nm,Nx,Nr = rawdata.shape
        Nth = (Nm+1)//2
        #if Ntheta is None or Ntheta < Nth:
        Ntheta = Nth
        fd = np.empty((Nx, Ntheta, Nr), dtype=np.complex128)

        fd[:, 0, :].real = rawdata[0, :, :]
        rawdatasw = np.swapaxes(rawdata, 0, 1)
        fd[:, 1:Nth, :].real = rawdatasw[:, 1::2, :]
        fd[:, 1:Nth, :].imag = rawdatasw[:, 2::2, :]

        fd = np.fft.fft(fd, axis=1).real
        mod_fd = np.swapaxes(fd,0,1)
        return mod_fd

    def mode_expansion(self):
        
        if self.theta is not None:
            F_total = self.mode_expansion_DFT()

        elif self.theta is None:
            F_total = self.mode_expansion_FFT()

        return F_total

    def getData(self):
        
        self.efield_names=["-Er","+Er","-El","+El","-Et","+Et"]
        self.rho_names=["-Rho","+Rho"]
        #conversion factors
        density_factor = (ncrit*e)*(1e12/1e6) #density normalized units --> pC/cm^3
        efield_factor = me*c*omega0/e/1e12 #E field normalized units --> TV/m           
        
        if self.field in self.rho_names:
            scalar_field_array = self.mode_expansion() * density_factor
       
        if self.field in self.efield_names:
            scalar_field_array = self.mode_expansion() * efield_factor
        
        return scalar_field_array
        
    def plot(self,vmin=None, vmax=None, cmap=None,saveas=None):
        

        if self.real_time>self.window_time:
            xmin = 1e6*(window_vx*c*(self.real_time - self.window_time + 5)*fs/omega0)
        else:
            xmin = 0
        xmax = xmin + (Lx*c_over_omega0*1e6)
    
        ymin = 0
        ymax = Lr*c_over_omega0*1e6
    
        fig, ax = plt.subplots()
    
        im=ax.imshow(self.getData().T ,cmap=cmap, origin='lower', 
           extent=(xmin, xmax, ymin, ymax),vmin=vmin,vmax=vmax,aspect="auto")
    
        cbar = plt.colorbar(im, ax=ax, shrink=1, aspect=20)
        plt.xlabel('x [um]')
        plt.ylabel('r [um]')
        if self.field in self.efield_names:
            plt.title(self.field[1:]+' [TV/m]  t='+str(self.real_time)+" fs")
        elif self.field in self.rho_names:
            plt.title(self.field[1:]+' [pC/cm^3]  t='+str(self.real_time)+" fs") 
        if saveas!=None:
            plt.savefig(saveas)
        plt.show()

    

#field = Field()