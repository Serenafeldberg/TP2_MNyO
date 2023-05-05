import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.animation import FuncAnimation
import math
from scipy.interpolate import interp1d

def E (mass, length, gravity, theta, omega):
    return T(mass, length, omega) + V(mass, gravity, length, theta)

def T (mass, length, omega):
    return (1/2)*mass*(length**2)*(omega**2)

def V (mass, gravity, length, theta):
    return mass*gravity*length - (mass*gravity*length)*(np.cos(theta))

def solve_pendulum(theta1_0, theta2_0, theta1_p0, theta2_p0, delta, t0, tf, N, length, gravity):
    # Definir las ecuaciones diferenciales del péndulo doble
    def f(t, y):
        theta1, theta1_p, theta2, theta2_p = y
        cos_delta = np.cos(delta)
        sin_delta = np.sin(delta)
        theta1_pp = ( -gravity/length*(2*np.sin(theta1)-np.cos(theta1-theta2)*np.sin(theta2)) - np.sin(theta1-theta2)*(theta2_p**2 + np.cos(theta1-theta2)*theta1_p**2)  )/(2 - (np.cos(theta1-theta2))**2)
        theta2_pp = ( -gravity/length*(np.sin(theta2)-2*np.cos(theta1-theta2)*np.sin(theta1)) + np.sin(theta1-theta2)*(theta1_p**2 + np.cos(theta1-theta2)*theta2_p**2)  )/(2 - (np.cos(theta1-theta2))**2)
        return [theta1_p, theta1_pp, theta2_p, theta2_pp]

    # Calcular el intervalo de tiempo
    dt = (tf - t0) / N

    # Crear las matrices para almacenar las soluciones
    theta1 = np.zeros(N)
    theta2 = np.zeros(N)

    theta1_p = np.zeros(N)
    theta2_p = np.zeros(N)

    # Establecer las condiciones iniciales
    theta1[0] = theta1_0
    theta2[0] = theta2_0

    theta1_p[0] = theta1_p0
    theta2_p[0] = theta2_p0

    # Usar el método de Runge-Kutta de cuarto orden para resolver las ecuaciones diferenciales
    for i in range(N-1):
        y = [theta1[i], theta1_p[i], theta2[i], theta2_p[i]]
        k1 = dt * np.array(f(t0 + i*dt, y))
        k2 = dt * np.array(f(t0 + i*dt + dt/2, y + k1/2))
        k3 = dt * np.array(f(t0 + i*dt + dt/2, y + k2/2))
        k4 = dt * np.array(f(t0 + i*dt + dt, y + k3))
        y = y + (k1 + 2*k2 + 2*k3 + k4) / 6
        theta1[i+1] = y[0]
        theta1_p[i+1] = y[1]
        theta2[i+1] = y[2]
        theta2_p[i+1] = y[3]

    return theta1, theta2, theta1_p, theta2_p

def anim (t, theta1, theta2, l1, l2, title ):
    fig, ax = plt.subplots(figsize=(5,5))
    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    ax.set_aspect('equal')
    ax.grid()

    # Creating the lines for the double pendulum
    line1, = ax.plot([], [], 'o-', lw=2, color='b')
    line2, = ax.plot([], [], 'o-', lw=2, color='r')

    theta1_interp = interp1d(t, theta1)
    theta2_interp = interp1d(t, theta2)

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
    anim = FuncAnimation(fig, update, frames=len(t), interval=10, blit=True, repeat=False)

    # Show the animation
    plt.title(title)
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)")
    plt.show()

def pendulum_trayectory (t, theta, label_):
    plt.plot(t, theta, label = label_)
    plt.legend()
    plt.title ("Trayectory of the pendulum")
    plt.xlabel("Time (s)")
    plt.ylabel("Angle (rad)")
    plt.show()

def pendulum_position (x, y, label_):
    plt.plot(x, y, label = label_)
    plt.legend()
    plt.title("Position of the pendulum")
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)")
    plt.show()

def plot_energy (t, theta, theta_p, length, gravity, label1, label2, label3):
    plt.plot(t, E(10, length, gravity, theta, theta_p), label = label1)
    plt.plot(t, T(10, length, theta_p), label = label2)
    plt.plot(t, V(10, 9.81, length, theta), label = label3)
    plt.title("Energy of the pendulum")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.legend()
    plt.show()
    

if __name__ == "__main__":
    
    # Parámetros
    length = 1.0
    theta1_0 = np.pi/6
    theta2_0 = np.pi/4
    theta1_p0 = 0.0
    theta2_p0 = 0.0
    delta = theta1_0 - theta2_0
    t0 = 0.0
    tf = 20
    N = 2000
    t = np.arange(t0, tf, (tf-t0)/N)

    # Resolver las ecuaciones diferenciales
    theta1, theta2, theta1_p, theta2_p = solve_pendulum(theta1_0, theta2_0, theta1_p0, theta2_p0, delta, t0, tf, N, length, 9.81)

    anim(t, theta1, theta2, length, length, "Double Pendulum animation")
    #animate_pendulum(np.linspace(t0, tf, N+1), theta1, theta2, length, length)

    x1_ = length*np.sin(theta1)
    y1_ = -length*np.cos(theta1)

    x2_ =  x1_ + length*np.sin(theta2)
    y2_ = y1_ - length*np.cos(theta2)

    pendulum_trayectory(t, theta1, "First Pendulum trayectory")

    pendulum_position(x1_, y1_, "First Pendulum Position")

    plot_energy(t, theta1, theta1_p, length, 9.81, "First Pendulum Total Energy", "First Pendulum Cinetic Energy", "First Pendulum Potential Energy")

    pendulum_trayectory(t, theta2, "Second Pendulum trayectory")

    pendulum_position(x2_, y2_, "Second Pendulum Position")

    plot_energy(t, theta2, theta2_p, length, 9.81, "Second Pendulum Total Energy", "Second Pendulum Cinetic Energy", "Second Pendulum Potential Energy")
    

