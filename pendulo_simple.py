import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def update(num, x, y, line):
    line.set_data([0, x[num]], [0, y[num]])
    return line,

def show_animation(L, t, x, y, title):
    fig, ax = plt.subplots()
    line, = ax.plot([], [], 'o-', lw=2)
    ax.set_xlim(-L, L)
    ax.set_ylim(-1.1*L, 1.1*L)
    ax.set_aspect('equal')
    ax.grid()
    
    ani = animation.FuncAnimation(fig, update, len(t), fargs=[x, y, line], interval=10, blit=True)
    #show animation
    plt.title(title)
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)")
    plt.show()

def E (mass, length, gravity, theta, omega):
    return T(mass, length, omega) + V(mass, gravity, length, theta)

def T (mass, length, omega):
    return (1/2)*mass*(length**2)*(omega**2)

def V (mass, gravity, length, theta):
    return mass*gravity*length - (mass*gravity*length)*(np.cos(theta))

def r (length, x):
    return np.array([length*math.sin(x), -length*math.cos(x)])

#cambio de variables
def z (theta_prime): 
    return theta_prime

def z_prime (x, w0):
    return -(w0**2)*math.sin(x)

def plotEnergy (t, mass, length, gravity, theta, omega, title):
    plt.plot(t, E(mass, length, gravity, theta, omega), label = "Total energy")
    plt.plot(t, T(mass, length, omega), label = "Kinetic energy")
    plt.plot(t, V(mass, gravity, length, theta), label="potential energy")
    plt.title(title)
    plt.legend()
    plt.show()

def plotTrayectory (x, y, title):
    plt.plot(x, y)
    plt.title(title)
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)")
    plt.show()

def euler_semi_implicito (y0, x0, t0, tn, step, gravity, length, mass): #f = z' = theta''
    w0 = math.sqrt(gravity/length)
    t_start = t0
    t_end = tn
    t = np.arange(t_start, t_end, step)
    positions = np.zeros(len(t)) #z'
    velocity = np.zeros(len(t)) #z
    positions[0] = y0 
    velocity[0] = x0 #velocidad

    for i in range(1, len(t)):
        velocity[i] = velocity[i-1] + step*z_prime(positions[i-1], w0)
        positions[i] = positions[i-1] + step*z(velocity[i])

    plt.plot (t, positions, 'b', label='trayectoria')
    plt.xlabel('t')
    plt.ylabel('theta')
    plt.show()
    
    _x = length*np.sin(positions)
    _y = -length*np.cos(positions)
    plotTrayectory(_x, _y, "Trayectory of simple pendulum - Euler")

    plotEnergy(t, mass, length, gravity, positions, velocity, "energy in Euler")
    

    #show_animation(length, t, _x, _y, "Motion of simple pendulum - Euler")
    return positions
    
# RUNGE KUTTA 4
g = 9.81 #gravitational acceleration (m/s^2)
L = 9.8 #length of the pendulum (m)
theta0 = math.pi/3  #initial angle (radians)
omega0 = 0.0  #initial angular velocity (radians)
m = 1.0 #mass of the pendulum
t0, tn = 0.0, 40.0
h = 0.1


def f(ti, yi):
    theta, omega = yi
    theta_dt = omega
    #omega_dt = -(g/L)*np.sin(theta)
    #omega_dt = z_prime(theta, math.sqrt(g/L))
    omega_dt = -(g/L)*np.sin(theta) 
    return np.array([theta_dt, omega_dt])

def rungeKutta4_method(gravity, length, theta0, omega0, mass, t0, tn, h):
    y0 = np.array([theta0, omega0])
    
    t = np.arange(t0, tn+h, h)
    y = np.zeros((len(t), 2)) # vamos a tener un arreglo de len(t) arreglos de dos elementos
    y[0] = y0
    
    for i in range(len(t)-1):
        k1 = h*f(t[i], y[i])
        k2 = h*f(t[i] + 0.5*h, y[i] + 0.5*k1)
        k3 = h*f(t[i] + 0.5*h, y[i] + 0.5*k2)
        k4 = h*f(t[i] + h, y[i] + k3)
        
        y[i+1] = y[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    
    #plot solution
    theta = y[:,0] 
    _x = L*np.sin(theta)
    _y = -L*np.cos(theta)

    #plotTrayectory(_x, _y, "Trayectory of simple pendulum - Runge Kutta 4")
    
    """
    plt.plot(t, theta, 'b')
    plt.xlabel('Time (s)')
    plt.ylabel('Angle (radians)')
    plt.title('Pendulum motion using Runge Kutta 4 method')
    plt.show()
    """
    
    omega = np.zeros(len(t))
    for i in range(len(t)):
        omega[i] = y[i,1]
    
    plotEnergy(t, m, L, g, theta, omega, "energy in RK4")
    
    #kinetic_energy = 0.5*m*(L*L)*(omega*omega)
    #potential_energy = -m*g*L*np.cos(theta) + m*g*L
    #total_energy = kinetic_energy + potential_energy
    """
    print("Total energy:\n")
    print(total_energy)
    
    print("Kinetic energy:\n")
    print(kinetic_energy)
    
    print("Potential energy:\n")
    print(potential_energy)
    """
    #return t, theta, _x, _y
    
    show_animation(L, t, _x, _y, "Motion of simple pendulum - RK4")
    
    return theta

def estabilidad_euler ():
    return

if __name__ == '__main__':
    euler_semi_implicito(math.pi/3, 0, 0, 10, 0.01, 9.81, 1, 10)
    rungeKutta4_method(9.81, 1, math.pi/6, 0, 10, 0, 100, 0.01)


    