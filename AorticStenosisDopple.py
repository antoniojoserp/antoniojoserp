#################################################################
#                                                               #
#   CW Doppler in customized aortic valve models                #
#                                                               #
#         coded by Antonio Romero, MD                           #
#   to be presented at:                                         #
#   European Society of Cardiology Congress 2022, Barcelona     #
#     http://shorturl.at/erLN9                                #
#                                                               #
#################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.patches as patches


class Valve:
    '''Valve object contains information relative to valves
    '''
    
    def ava_to_rad (ava):
        '''translates aortic valve area (cm2) into orifice radius (mm)
        d=Sqr(AVA/Pi)'''
        return (  (ava / np.pi)**0.5 ) * 10
    
    
    def __init__(self, l1=10, l2=10, length=5, rad=10, ava=0.8, 
        thick=1, shape='cusp', wall_thick=2, trim_factor=0):
        '''l1 = LVOT length, l2 = distal aortic wall, length = valve length
        (from base to tip)
        '''
        
        self.l1=float(l1) 
        self.l2=float(l2)
        if not self.l2 > 0: self.l2 = self.l1
        self.length=float(length)
        self.rad=float(rad)
        self.orif=float( Valve.ava_to_rad(ava) )
        self.thick=float(thick)
        self.shape=str(shape)
        self.wall_thick = float(wall_thick)
        self.trim_factor = float(trim_factor)
        self.y_range = float(self.length + self.l1)
        self.v_max = float()
        self.alfa = float()
        return
       
    def line (self):
        '''returns inner contour points  (cordinates) of the valve'''
        puntos=10  # number of ponts

    
        valv_xlen = self.rad - self.orif # valve width from base to orifice
        puntos_x=np.empty(puntos, dtype='float')

        for celda in enumerate ( puntos_x):
            i = celda[0]
            interval = i/(puntos-1)
            puntos_x[i]= self.rad - (valv_xlen)*interval

        if self.shape=='cusp':
            
            factor = self.trim_factor
            #trims a % of the centre of thr cusp; as indicted in trim_factor


            _x_len = valv_xlen/(1-factor)
            _x0 = self.rad - _x_len #punto teorico en que _x = 0 
            _x = lambda x:  x - _x0
            _y_len = (_x_len)**2  #amplitud de la curva (restar el valor en la punta)
            _x1=_x(self.rad-valv_xlen) #punta de la valva
            y_len =  _y_len - _x1**2   #amplitud de la valva (para el ajuste) entre base y punta

            puntos_y = np.array([ ( _x(x) )**2 for x in puntos_x], dtype='float')
            puntos_y = puntos_y*self.length/y_len -self.length  #ajusta la amplitud de los puntos para
            #que coincidan con la longitud de la estenosis

        if self.shape =='line':
            puntos_y= - (self.rad-puntos_x)*(self.length/valv_xlen)

        return np.array([puntos_x, puntos_y])

    def contour (self):

        line_1 = self.line () #valve upper (proximal) contour
        line_2 = line_1[:, ::-1].copy() # valve lower (distal) contour
        line_2[1] = line_2[1] - self.thick
        return np.concatenate([line_1,line_2], axis=1)

    def right_wall (self):
        '''creates right wall and valve contour coordinates, the left will be
        created as a mirror 
        '''
        
        l1=self.l1
        l2=self.l2
        
        valve_points=self.contour()
        valve_points[1]=valve_points[1] - l1
        roof = np.array([[self.rad+self.wall_thick, self.rad],
            [0,0]], dtype='float')
        base= - (l1 + self.length + l2) # base y coordinates
        basement = np.array([[self.rad, self.rad+self.wall_thick],
            [base,base]], dtype='float')
        

        polygon = np.concatenate([roof, valve_points, basement], axis=1)
        return polygon

    def plot(self):
        '''plots the valve'''
        polygon= self.right_wall()
        plt.fill(polygon[0], polygon[1], 'k')
        plt.fill(-polygon[0], polygon[1], 'k')
        plt.xlim(-25,25)
        plt.ylim ( -(self.l1 + self.l2 +self.length + 5), 10)
        _ax = plt.gca()
        plt.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = False, bottom = False)
        for key, spine in _ax.spines.items():
            spine.set_visible(False)
        
        return
    
    def radio (self, y):
        ''' for a given section or level (height, y coordinates) from LVOT
        base (y=0) to valve cusp (y= -(l1+length), returns the radius'''

        l1=self.l1
        radio=float()
        if y > 0: radio= np.nan
        elif y >= -l1: radio = self.rad
        # watch out, negative (<0) y values. While on LVOT, radius = self.rad
        elif y >= - (l1 + self.length):
            contour = self.line().T # needs to be trasposed
            y_ = y - (-l1)   # y_: _y coordinates referenced to valve base
            for i in range( len(contour)-1 ):
                '''this loop selects the correct rangeof y points; when finds
                the correct interval, interpolates y and radius values'''

                x1,y1 = contour[i]
                x2,y2 = contour[i+1]
                if y_ < y2: continue
                elif y_ >= y2:  
                    radio = x1 + ( (y_-y1)/(y2-y1) )*(x2-x1)
                    break
        
        #Note:
        # The following condition (y_ >= y2) would theooretically not be
        #neccessary, as the program should not get y values beyond valve tip. In
        #practice, float point migration errors may produce inapropriate values
        #of the radius'''
        
        
        elif y >= - (l1 + self.length + 2): 
            radio = self.orif
                
        else: radio = np.nan
            
        return radio

    
class Wave:
    '''Created the wave object, which returns the a and y coordinates of the
    higher velocities Doppler tracing during entire systole.
    Lower velocities would be inferred  proportionally.
    Some parameters of the Doppler tracing, sucha as aceleration time and
    ejection time (wave width), will be inferred from the max velocity'''
    
    def calc_t_acel (self):
        if self.peak <= 700:  result = 70 + 2*(self.peak/100)**2
        if self.peak > 700:   result = 168
        return result
    
    def calc_t_eyec (self):
        if self.peak <700: result = 240 + 0.143*self.peak
        if self.peak >= 700: result = 340
        return result
        
    def __init__(self, peak=450, time=600):
        self.peak=float(peak)
        self.t_acel=Wave.calc_t_acel(self)
        self.t_eyec=Wave.calc_t_eyec(self)
        self.t_decel= self.t_eyec - self.t_acel
        self.time=float(time)
        
        #determinar las agrupaciones de bandas (necesario para la animacion)
        
    
    def wavepoint (self, ms, q_eyec=80):
        ''' returns velocity for a given milisecond'''
        v=float()
        t = ms - (q_eyec + self.t_acel) # =0 at the curve peak
        k_1 = self.peak/self.t_acel**2
        k_2 = self.peak/self.t_decel**2
        if t< (- self.t_acel) : v = np.nan
        elif t <= 0: v=  k_1*t**2 - self.peak
        elif t < self.t_decel: v= k_2*t**2 - self.peak
        else: v = np.nan
        return v
    
    def get_outer_wave (self, sep_ms=1):
        '''returns an np array with the outer (max velocity) Doppler
        coordinates, by default, spaced 1 ms in time'''
        
        timeline=np.arange(0,self.time+1, sep_ms, dtype='float')
        velocidades=np.array([self.wavepoint(ms) for ms in timeline], dtype='float')
        
        return velocidades
    
    def plot (self):
        '''plots outer contour (maximum velocities) of doppler tracing'''
        plt.plot(self.get_timeline(), self.get_outer_wave())
        plt.xlim(0, 500)
        return
    
    def get_timeline (self, sep_ms=1):
        ''' return an one-dimensional array with only  time coordinates'''

        return np.arange(0,self.time+1, sep_ms)
        
    def wave (self, veloc, sep_ms=1):
        '''returns Doppler tracing produced by a given velocity (in some way
        analog to pulsed-Doppler). Trace is inferred from the outer contour
        in proportion'''
        escala = veloc/self.peak
        return self.get_outer_wave(sep_ms)*escala
    
  
class Simulation:
    '''Simulation object stores main simulation parameters, to be used
    in the plot of one or varius valves.
    -bands = number of equally spaced lines (levels, sections) which are
    considered in the simullation
    - groups: used in some features, such as animations or pulsed Doppler
    simualtion. If <> 1, the bands are distributed into the specified number
    of groups
    '''
    
    
    def band_grouping(self):
        ''' Distributes bands between groups; in case of non exact division,
        remains are distributed among the first groups.
        Returns list (group_limits) which contains the cutoff limits of groups'''

        b = self.bands
        g = self.band_groups
        n = b // g # bands per group (initial approach)
        resto = int( b % g)
        group_limits=list()
        ''' The list group_limits stores the value of the index [0/bands-1] of
        the first band of the NECT group. eg: group_limits[0] indicates the
        first band to be included in grouop 1
        '''

        limit=0
        for i in range(g):
            limit = limit + n # adds n to the previous group limit
            if resto > 0:
                limit += 1
                resto -= 1
            group_limits.append(limit)

        return group_limits
    
    def __init__(self, bands=80, basal_velocity=110, band_groups=10, method='scatter', valves=[Valve()],
                 seed=None, dispersion=0.05, sim_alfa=None):
        self.bands=int(bands)
        self.band_groups=int(band_groups)
        self.basal_velocity=float(basal_velocity)
        self.group_limits= list( self.band_grouping() )
        self.method=method
        self.valves=valves
        self.dispersion = float(dispersion)
        self.sep_ms = 1
           
        # seed parameter for np.random engine
        if seed != None: self.seed=seed
        else:
            self.seed = int (
                np.random.randint(0,9999,1 )
                )
            
        for valve in self.valves:
            valve.v_max = Simulation.v_conteq (valve.orif, valve, basal_velocity)
            # bright adjustement (alfa)

            if sim_alfa == None:
                valve.alfa = 12 / self.bands
            else: valve.alfa = sim_alfa
                        
        self.v_max_scale = max([valve.v_max for valve in self.valves])
        
        
        
        
    def v_conteq (radio, valve=Valve(), basal_velocity=110):
            '''for a given radius and valve object anda basal velocity,
            returns velocity based on continuity equation.
            '''
            basal_v = basal_velocity
            rad_LVOT = valve.rad
            
            return (basal_v * rad_LVOT**2) / radio**2
        

    def create_fig (self): 
        
        valves=self.valves
        
        if len(valves) ==1:      
            # when only a valve is given, the figure is arranged horizontally
            fig , ax = plt.subplots(1,2)
            ax = ax[:,np.newaxis]  
       
        else:
            # else, pairs of valve schemes and Doppler
            #tracings are arranged vertically
            fig , ax = plt.subplots( 2, len(valves) , sharey='row')
            
        for i, ax_ in enumerate(ax[1,:]):
            ax_.set_facecolor('k')
            plt.sca(ax_)
            plt.ylim(-self.v_max_scale*1.2, 0)
            plt.xlim(50,450)
            if i == 0:
                plt.xlabel ('Time from QRS onset (ms)', loc = 'left')
                plt.ylabel ('Velocity (cm/seg)', loc = 'center')
            ax_.tick_params(axis='both', which='major', labelsize=8)
            ax_.set_xticks([200, 400])

        for i, valve in enumerate(valves):
            plt.sca(ax[0,i])
            valve.plot()
            if len (valves)>1: # inserts figure letter if there are more 
                abecedario = "abcdefghijklmnopqrstuvwxyz"
                plt.title(abecedario[i].upper(), loc='left', fontsize=20)
        
      
        return fig, ax
    
    
    def calc_band(self, n, valve ):
        ''' for a given band, returns the ROI coordinates in the valve scheme
        and the scatter points for CW tracing'''
        interval= valve.y_range/self.bands
        wave = Wave(valve.v_max)
                        
        y_1 = - n*interval      #upper limit coords
        y_mid = - (0.5 + n*interval) #medium point
        y_2 = - (n+1)*interval  #lower point
            

        
        radio1= valve.radio(y_1)
        radio2 = valve.radio(y_2)
        radio = radio2   # used for velocity calculations
            
            
        vel = Simulation.v_conteq(radio, valve, self.basal_velocity)
                            
        # right coordinates of the region of interest in the scheme
        x_coords = np.array([radio1, radio2], dtype='float' )
        y_coords = np.array([y_1, y_2])
        
        contour = np.array([x_coords, y_coords], dtype='float')
        
        # doppler scatter points
        sep_ms=self.sep_ms  # spacing (timing) between points
        x_points = wave.get_timeline(sep_ms)
        x_points = x_points [~(np.isnan(x_points) )]
        y_points = wave.wave(vel, sep_ms) 
        
        # adds random dispersion
        aleatory_dispersion = np.random.RandomState(seed=self.seed + n).normal(1,
                self.dispersion, len(y_points))
        y_points = y_points * aleatory_dispersion
        points=np.array([x_points, y_points], dtype='float')
        
        return contour, points

            
                                 
    
    def calc_group (self, group, valve=None):
        '''returns data (rois, velocity points) for the specified group.
        If the group number equals or surpass actual number groups, returns
        complete data for the figure
        '''
        if valve == None: valve = self.valves[0]
        bands=self.bands
        limits=list (self.group_limits )
                
        v_max = Simulation.v_conteq (valve.orif, valve, self.basal_velocity)
        
        if group >= self.band_groups:
            band_0 = 0
            band_n = self.bands
            
        else:
            band_0 = 0 if group == 0 else limits[group-1]
            band_n = limits[group]
        
        cont= np.array([[],[]])   # ROI (polygon) contours 
        points = np.array([[],[]])
        
                 
        for n in range (band_0, band_n):
            new_cont, new_points = self.calc_band(n, valve)
            
            #for bands 1 and following, the first point is redundant and is
            # removed
            if n > band_0: new_cont = new_cont[:,1:]
            
            cont = np.append(cont, new_cont, axis=1)
            points = np.append(points, new_points, axis=1)
                        

        
        # completes the ROI coordinates
        left_cont = cont[:,::-1].copy() # inverts the columns order
        left_cont[0,:] = left_cont[0,:]*(-1)  #create negative x values (mirror)
        polig = np.append(cont, left_cont, axis=1)

                
        return polig, points
    
    def plot_group(self, group, line=False, fondo_negro=True):
        
        bands=self.bands
        basal_v=self.basal_velocity
        columnas = len(self.valves)
                        
        fig, ax = self.create_fig()
        for i , valve in enumerate (self.valves):
            wave = Wave(valve.v_max)
            polig, points = self.calc_group(group, valve)
                      
            plt.sca(ax[0,i])
            plt.fill(polig[0], polig[1], color = '#F39D9A')
            
            plt.sca(ax[1,i])
            point_number =  len (points[0])
            
            
            if fondo_negro==False:
                ax[1,i].set_facecolor('w')
                color_puntos='k'
            else: color_puntos = 'w'
                
            plt.scatter(points[0], points[1], c=color_puntos, marker = '.', alpha = valve.alfa)
            
            # show LVOT velocity as reference (dotted line)
            if line == True and group >= self.band_groups:
                for v, valve in enumerate (self.valves):
                    
                    plt.sca(ax[1,i])    
                    plt.plot(wave.get_timeline(1), wave.wave(basal_v, 1),
                        'r--', linewidth=3, alpha=0.2)          
        
        return
                
               
    def plot (self, line=False, fondo_negro=True):
        '''plots complete figure'''
        self.plot_group(self.band_groups +1, line, fondo_negro=fondo_negro)
              
        return
        

def ask_bands():

    try:
        bands = input('\t- number of sections (levels) per simulation (100-800)? (default=600):')
        if bands == '': bands = 600
        bands = int(bands)
        if bands <100 or bands > 800:
            print ('\t- sorry, this code only accepts between 100 and 800 sections')
            bands = ask_bands()
        
    except:
        print('\n\t-not a valid value \n')
        bands= ask_bands()
    return bands

def ask_basal_vel():
    try:
        b = input('\n\t- LVOT peak velocity (cm/s)? (default=110 cm/s):')
        if b == '': b = 110
        b = float(b)
    except:
        print('\n\t-not a valid value \n')
        b= ask_basal_vel()
    return b

def ask_alfa():
    try:
        a = input('\n\t- force brightness (0-1)? (default: auto)')
        if a == '': a = None
        else:
            a = float(a)/5
            if a <0 or a > 0.2:
                print('\t - only accepted values between 0 and 1')
                a = ask_alfa()
    except:
        print('\n\t-not a valid value \n')
        a= ask_alfa()
    return a

def ask_line():
    try:
        _ = input('\n\t- Show LVOT velocity line (Y/N)? (default=Y)')
        _ = str(_) 
        if _ == '' or _ .lower in ['y', 'yes', 'si', 's', 'sí', 'True']:
             _ = True
        else: _ = False
    except:
        print('\n- not a valid value \n')
        _ = ask_line()
    return _
    
def ask_valves():
    try:
        _ = input ('\n\t- number of figures plotted simultaneously (1-3)? (default:1):')
        if _ == '': _=1
        else:
            _ = int(_)
            if _ <1 or _ > 3:
                print('\t - sorry, accepted values between 1 and 3')
                _ = ask_valves()
    except:
        print('\n\t-not a valid value \n')
        _ = ask_valves()
    return _

def ask_range(linf, lsup, text=None, default=0):
    
    if text == None: text = f'- enter value ({linf}, {lsup})'
    try:
        _ = input(text)
        if _ == '': _ = default
        else:
            _ = float(_)
            if _ < linf or _ > lsup:
                print(f'\t - sorry, accepted values between {linf} and {lsup}')
                _ = ask_range(linf, lsup, text, default)
    except:
        print('\n\t-not a valid value \n')
        _ = ask_range(linf, lsup, text, default)
    return _

def ask_shape():
    try:
        _ = input('\n - Shape of the cusps? 1=cusp (default), 2 = line')
        if _ == '' or _ == 'cusp': _ = 1
        elif _ == 'line': _ = 2
        elif int(_) in [1,2]: _ = int(_)
        else: raise Exception
    except:
        print('\n- not a valid value \n')
        _ = ask_shape()
    dic={1:'cusp', 2:'line'}
    
    return dic[_]




    
            
#    
#  MAIN    
#

print('CW Doppler simulation:\n\n\-Basic parameters:')


bands = ask_bands()
basal_velocity= ask_basal_vel()
basal_line = ask_line()
alfa = ask_alfa()
v = ask_valves()
disp = ask_range(0,100, '\n\t- random dispersion standard deviation (as % of actual values) (0-100) (default:5%)', 5)/100

valves = list()
for i in range(v):
    print (f'######################\n VALVE NUMBER: {i+1}\n')
    v_length = ask_range(1,10,'\n- Valve length (mm) (1-10)? (default = 5)',5)
    l1 = ask_range(5,20, '\n - LVOT length (mm) (5-20): (default=10)', 10)
    r = ask_range(5,15, '\n - LVOT radius (mm) (5-15): (default=10)', 10)
    ava = ask_range(0.2,2, '\n - Aortic valve area (cm2) (0.2-2): (default=0.8)', 0.8)
    shape = ask_shape()
    val = Valve(l1=l1, length = v_length, rad=r, ava=ava, shape= shape)
    valves.append(val)
    
sim = Simulation(bands=bands, basal_velocity = basal_velocity, band_groups=1, valves= valves, sim_alfa=alfa, dispersion=disp)
sim.plot_group(1, line=basal_line)
plt.show()
    
