import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation


class Gardner_burger_solver:
    def __init__(self, N, u_0, x, dt, alpha, rhi):
        self.x = x
        dx = x[1] - x[0]
        self.k = 2*np.pi*np.fft.fftfreq(N, d = dx)
        self.u_0f = np.fft.fft(u_0) # the solution on the frequency domain
        self.t = 0.0
        self.dt = dt
        #parameters
        self.alpha = alpha
        self.rhi = rhi
    def RK4(self, u):
        k = self.k
        alpha = self.alpha
        u1 = 2j*alpha*k*np.fft.fft(np.real(np.fft.ifft(u))**3)
        return  u1 
    def get_u(self):
        return np.real(np.fft.ifft(self.u_0f))
    def next_step(self):
        #parameters
        x = self.x
        rhi = self.rhi
        k = self.k
        dt = self.dt
        u = self.u_0f
        t = self.t
        self.t = t + dt
        # split step fourier
        u1 = u*np.exp((1j*k**3 - rhi*k**2)*dt)
        k1 = self.RK4(u1)
        k2 = self.RK4(u1 + dt*k1/2)
        k3 = self.RK4(u1 + dt*k2/2)
        k4 = self.RK4(u1 + dt*k3)
        u = u1 + dt*(k1 + 2*k2 + 2*k3 + k4)/6
        self.u_0f = u

def step_function(ul,ur):
    u_0 = np.zeros(N)  
    for i in range(len(x)):
        if x[i] < 0:
    	    u_0[i] = ul
        else:
    	    u_0[i] = ur
    return u_0

def plot_the_graph(new_wave, t_f):
    dt = new_wave.dt
    for i in range(round(t_f/dt)):
        new_wave.next_step()
        print('t = ', (i+1)*dt)
    np.savetxt('region_1.txt', (new_wave.get_u(),new_wave.x))
    plt.plot(new_wave.x, new_wave.get_u(), 'b')
    plt.show()

#initiate the simulation
def init():
    line.set_ydata([])
    return line,

# The function of animation
def animation_frame(i):
    new_wave.next_step()
    line.set_ydata(new_wave.get_u())
    print("t=", new_wave.t)
    return line,

# definition of the coordinate domain
N = 2**17
xmax = 10000
xmin = -10000
dx = (xmax - xmin)/N
x = np.linspace(xmin, xmax, N)

# The definition of the parameters
alpha = 0.20
rhi = 0.01

# mKDV initial condition 1 soliton
ul = 1.0
ur = 0.2  #-ul + 0.09
u_0 = step_function(ul, ur)


# The time simulation (just for the plot)
t_f = 4000
dt = 0.01

new_wave = Gardner_burger_solver(N, u_0, x, dt, alpha, rhi)

print("V= ", -2*alpha*(ul**2+ur*ul+ur**2))

#plot_the_graph(new_wave, t_f)

#Create the animation
#graph configuration

fig, ax = plt.subplots()
ax.set_ylabel("u(x,t)")
ax.set_xlabel("x")
ax.set_ylim(np.min(u_0)- 0.02, 4.0)
#ax.set_xlim(np.min(x), np.max(x))

line, = ax.plot(x, new_wave.get_u())
animation = FuncAnimation(fig, func=animation_frame, init_func=init, frames=np.arange(0,t_f, dt), interval=5, blit = True)

plt.show()
