import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from matplotlib.animation import FuncAnimation

# Constants
m1, m2 = 6.4, 6.4 #mass of the pendulums
l1, l2 = 1, 1 #length of pendulums in metres
g = 9.81 #acceleration due to gravity in m/s^2

# Initial conditions
theta1_0, theta2_0 = np.pi/3, np.pi/2 # initial angle of the pendulums
omega1_0, omega2_0 = 0, 0 #initial angular velocity of pendulums (radians)

t0 = 0.0 #initial time in seconds
tf = 20.0 #final time in seconds
h = 0.01 #time step size in seconds

# Define time array
t = np.arange(t0, tf + h, h)

# Function for the double pendulum system
def double_pendulum(t, y):
    theta1, theta2, omega1, omega2 = y
    
    dtheta1_dt = omega1
    dtheta2_dt = omega2
    
    #domega1_dt = (-g*(2*m1 + m2)*np.sin(theta1) - m2*g*np.sin(theta1 - 2*theta2) - 2*np.sin(theta1 - theta2)*m2*(omega2**2*l2 + omega1**2*l1*np.cos(theta1 - theta2))) / (l1*(2*m1 + m2 - m2*np.cos(2*theta1 - 2*theta2)))
    #domega2_dt = (2*np.sin(theta1 - theta2) * (omega1**2*l1*(m1 + m2) + g*(m1 + m2)*np.cos(theta1) + omega2**2*l2*m2*np.cos(theta1 - theta2))) / (l2*(2*m1 + m2 - m2*np.cos(2*theta1 - 2*theta2)))
    domega1_dt = ( -g/l1*(2*np.sin(theta1)-np.cos(theta1-theta2)*np.sin(theta2)) - np.sin(theta1-theta2)*(omega2**2 + np.cos(theta1-theta2)*omega1**2)  )/(2 - (np.cos(theta1-theta2))**2)
    domega2_dt = ( -g/l1*(np.sin(theta2)-2*np.cos(theta1-theta2)*np.sin(theta1)) + np.sin(theta1-theta2)*(omega1**2 + np.cos(theta1-theta2)*omega2**2)  )/(2 - (np.cos(theta1-theta2))**2)
    
    
    return dtheta1_dt, dtheta2_dt, domega1_dt, domega2_dt



# Solving the system using Runge Kutta method
sol = solve_ivp(double_pendulum, [t[0], t[-1]], [theta1_0, theta2_0, omega1_0, omega2_0], method='RK45',t_eval=t)

# Interpolating to get smoother animation
theta1_interp = interp1d(t, sol.y[0])
theta2_interp = interp1d(t, sol.y[1])

# Creating the figure and axis
fig, ax = plt.subplots(figsize=(5,5))
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)
ax.set_aspect('equal')
ax.grid()

# Creating the lines for the double pendulum
line1, = ax.plot([], [], 'o-', lw=2, color='b')
line2, = ax.plot([], [], 'o-', lw=2, color='r')

# Function to update the lines for the double pendulum at each time step
def update(i):
    theta1 = theta1_interp(t[i])
    theta2 = theta2_interp(t[i])
    
    x1 = l1*np.sin(theta1)
    y1 = -l1*np.cos(theta1)
    x2 = x1 + l2*np.sin(theta2)
    y2 = y1 - l2*np.cos(theta2)
    
    line1.set_data([0, x1], [0, y1])
    line2.set_data([x1, x2], [y1, y2])
    
    return line1, line2,

# Creating the animation
anim = FuncAnimation(fig, update, frames=len(t), interval=20, blit=True)

# Show the animation
plt.show()