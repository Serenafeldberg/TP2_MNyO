import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def update(num, x, y, line):
    line.set_data([0, x[num]], [0, y[num]])
    return line,

def show_animation(L, t, x, y):
    fig, ax = plt.subplots()
    line, = ax.plot([], [], 'o-', lw=2)
    ax.set_xlim(-L, L)
    ax.set_ylim(-1.1*L, 1.1*L)
    ax.set_aspect('equal')
    ax.grid()
    
    ani = animation.FuncAnimation(fig, update, len(t), fargs=[x, y, line], interval=10, blit=True)
    #show animation
    plt.title("Motion of simple pendulum")
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)")
    plt.show()

def E (mass, length, gravity, theta):
    return T(mass, length, theta) + V(mass, gravity, length, theta)

def T (mass, length, theta):
    return (1/2)*mass*(math.pow(length,2))*(math.pow(theta, 2))

def V (mass, gravity, length, theta):
    return mass*gravity*length - (mass*gravity*length)*(math.cos(theta))

def r (length, x):
    return np.array([length*math.sin(x), -length*math.cos(x)])

#cambio de variables
def z (theta_prime): 
    return theta_prime

def z_prime (x, w0):
    return -(w0**2)*math.sin(x)

def euler_explicito (y0, x0, t0, step, N_iter, f, gravity, length): #f = z' = theta''
    mass = 3
    w0 = math.sqrt(gravity/length)
    t_start = t0
    t_end = N_iter*step
    t = np.arange(t_start, t_end, step)
    positions = np.zeros(len(t)) #z'
    velocity = np.zeros(len(t)) #z
    positions[0] = y0 
    velocity[0] = x0 #velocidad

    for i in range(1, len(t)):
        velocity[i] = velocity[i-1] + step*f(positions[i-1], w0)
        positions[i] = positions[i-1] + step*velocity[i]

    plt.plot (t, positions, 'b', label='trayectoria')
    plt.xlabel('t')
    plt.ylabel('theta')
    plt.show()
    
    _x = length*np.sin(positions)
    _y = -length*np.cos(positions)
    plt.plot(_x, _y)
    plt.title("Trayectory of simple pendulum with Euler")
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)")
    plt.show()
    
    show_animation(length, t, _x, _y)


if __name__ == '__main__':
    euler_explicito(math.pi/3, 0, 0, 0.01, 1000, z_prime, 9.81, 1)
    