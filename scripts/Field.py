from extractdata import *
from smilei import *
import numpy as np

class Field:
    
    def _validate(self):
        
        for im in self.modes:
            if im<0 or im > (self.AM_modes-1):
                raise IOError('Invalid mode : {} '.format(self.modes))

        if os.path.isfile(self.h5path) is False:
            raise IOError('{} file/path does not exist'.format(self.h5path))
        
        if self.field[1:] not in ['El','Er','Rho']:
            raise IOError('Field {} not available in simulation dump'.format(self.field)) 
    
    def grp(self,dir,flag=None): #grp stands for group in a HDF5 file
        
        if flag=='ar':
            return np.array(self.h5item.get(dir))
        elif flag=='l':
            return list(self.h5item.get(dir))
        else:
            return self.h5item.get(dir)
    
    def __init__(self, h5item, field, moving=None, timesteps=None, modes=None,theta=None):
        
        self.h5path = h5item
        if os.path.isfile(self.h5path) is False:
            raise IOError('{} file or path does not exist'.format(self.h5path))
        
        self.field = field
        if self.field[1:] not in ['El','Er','Rho']:
            raise IOError('Field {} not available in simulation dump'.format(self.field)) 
        
        self.moving = moving
        self.user_timesteps = timesteps

        self.window_time = (MovingWindow.time_start)/fs
        self.AM_modes = Main.number_of_AM
        
        #self._validate()  #to check the validity of user arguments 
        self.h5item = gethdf5(h5item)
        self.timestep_list = [int(x) for x in list(self.grp(dir='data',flag='ar'))]

        if modes is not None:
            self.modes = modes
        else:
            self.modes = [i for i in range(self.AM_modes)]
        for im in self.modes:
            if im<0 or im > (self.AM_modes-1):
                raise IOError('Invalid mode : {} '.format(self.modes))
            
        if theta is not None:
            self.theta = theta
        else:
            self.theta = 0 

        if len(self.modes) == 0:
             raise IOError('Invalid mode : {} '.format(self.modes))
        for im in self.modes:
            if im<0 or im > (self.AM_modes-1):
                raise IOError('Invalid mode : {} '.format(self.modes))

    def getTime(self,usertime):    
        
        #Finding the simulation_timestep[sim_timesteps] and real time[real_time]
        if usertime is not None:
            int_timesteps = min(self.timestep_list, key=lambda x: (abs(x - usertime), x))
        else:
            int_timesteps = self.timestep_list[-1]
        
        sim_timesteps = str(int_timesteps).zfill(10)   #timestep format in the data dump
        real_time=int(sim_timesteps)*dt/fs    #timestep into femto-seconds unit
        
        return [int_timesteps,sim_timesteps,real_time]

    def getArray(self, timestep=None):

        if timestep is None:
            timestep = self.getTime()
        #self.field_name,self.field_array,self.real_field,self.imag_field=[],[],[],[]
        array_list=[]
        for mode in range(self.AM_modes):
        
            #self.field_name.append(self.field[1:]+"_mode_"+str(mode))
            #self.field_array.append(self.h5item['data'][self.sim_timesteps][self.field_name[mode]]*int(self.field[0]+"1"))
            #self.real_field.append(self.field_array[mode][:,::2])
            #self.imag_field.append(self.field_array[mode][:,1::2])

            field_name = self.field[1:]+"_mode_"+str(mode)
            field_array = self.grp(dir='data/'+timestep+'/'+field_name,flag='ar')*int(self.field[0]+"1")
            #field_array = self.h5item['data'][self.sim_timesteps][field_name]['data']*int(self.field[0]+"1")
            field_array_shape = field_array.shape
            reshaped_array = field_array.reshape(field_array_shape[0], field_array_shape[1]//2, 2)
            new_array = reshaped_array[:, :, 0] + 1j * reshaped_array[:, :, 1]
            array_list.append(new_array)
        
        mod_data= np.stack(array_list,axis=0)     #Modified array of shape (no.of.modes, Nx, Nr)
        return mod_data             

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
    
    def mode_expansion_DFT(self,timestep):
        
        mod_data=self.getArray(timestep=timestep)
        F = np.zeros_like(np.real(mod_data[0]))                       
        
        '''for m in self.modes:
            F += self.mod_data[m]*np.exp(-1j*m*self.theta)
        '''
        for m in self.modes:
            F += np.real(mod_data[m])*np.cos(m*self.theta)+np.imag(mod_data[m])*np.sin(m*self.theta)

        return F   
    
    def mode_expansion_FFT(self,timestep):

        rawdata=self.getArray(timestep=timestep)
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

    def mode_expansion(self,timestep):
        
        if self.theta is not None:
            F_total = self.mode_expansion_DFT(timestep=timestep)

        elif self.theta is None:
            F_total = self.mode_expansion_FFT(timestep=timestep)

        return F_total

    def getAxis(self, axis,timestep=None):
        
        timestep=self.getTime(usertime=timestep)[1]
        Nx,Nr = self.getArray(timestep=timestep)[0].shape

        if axis is None:
            raise IOError("Invalid axis")
        elif len(axis)>1:
            raise IOError("Only one axis at a time")
        
        if axis == "x":
            x_min = np.array(self.grp(dir='data/'+timestep).attrs['x_moved'])*c_over_omega0*1e6  
            x_max = x_min+Lx*c_over_omega0*1e6
            x_axis = np.linspace(x_min, x_max, Nx-1)
            return x_axis    
        elif axis == "r":
            r_max = Lr*c_over_omega0*1e6
            r_axis = np.linspace(0, r_max, Nr-1)
            return r_axis
        else:
            raise IOError("Invalid axis")
    
    def getData(self,timestep=None):
        
        timestep = self.getTime(usertime=timestep)[1]
        #self.efield_names=["-Er","+Er","-El","+El","-Et","+Et"]
        #self.rho_names=["-Rho","+Rho"]
        #conversion factors
        density_factor = (ncrit*e)*(1e12/1e6) #density normalized units --> pC/cm^3
        efield_factor = me*c*omega0/e/1e12 #E field normalized units --> TV/m           
        
        if self.field[1:] == "Rho":
            scalar_field_array = self.mode_expansion(timestep=timestep) * density_factor
       
        if self.field[1:2] == "E":
            scalar_field_array = self.mode_expansion(timestep=timestep) * efield_factor
        
        return scalar_field_array

    
    def plot(self,vmin=None, vmax=None, cmap=None,saveas=None):

        if self.real_time>self.window_time:
            xmin = self.h5item.get('data/'+self.sim_timesteps).attrs['x_moved']*c_over_omega0*1e6
        else:
            xmin = 0
        xmax = xmin + (Lx*c_over_omega0*1e6)
    
        ymin = 0
        ymax = Lr*c_over_omega0*1e6
    
        fig, ax = plt.subplots()
    
        im=ax.imshow(self.getData().T ,cmap=cmap, origin='lower', 
           extent=(xmin, xmax, ymin, ymax),vmin=vmin,vmax=vmax,aspect="auto")
    
        plt.colorbar(im, ax=ax, shrink=1, aspect=20)
        plt.xlabel('x [um]')
        plt.ylabel('r [um]')
        if self.field in self.efield_names:
            plt.title(self.field[1:]+' [TV/m]  t='+str(self.real_time)+" fs")
        elif self.field in self.rho_names:
            plt.title(self.field[1:]+' [pC/cm^3]  t='+str(self.real_time)+" fs") 
        if saveas!=None:
            plt.savefig(saveas)
        plt.show()

