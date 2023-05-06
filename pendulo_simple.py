import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
#import sympy as sp

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
    plt.xlabel("Time (s)")
    plt.ylabel("Energy")
    plt.title(title)
    plt.legend()
    plt.show()

def plotTrayectory (x, y, title):
    plt.plot(x, y)
    plt.title(title)
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)")
    plt.show()

def plotAngles(t, theta, title):
    plt.plot(t, theta)
    plt.title(title)
    plt.xlabel("Time (s)")
    plt.ylabel("Angle (rad)")
    plt.show()

def euler_semi_implicito (y0, x0, t0, tn, step, gravity, length, mass): #f = z' = theta''
    w0 = math.sqrt(gravity/length)
    t_start = t0
    t_end = tn
    t = np.arange(t_start, t_end+step, step)
    positions = np.zeros(len(t)) #z'
    velocity = np.zeros(len(t)) #z
    positions[0] = y0 
    velocity[0] = x0 #velocidad

    for i in range(1, len(t)):
        velocity[i] = velocity[i-1] + step*z_prime(positions[i-1], w0)
        positions[i] = positions[i-1] + step*z(velocity[i])
    
    _x = length*np.sin(positions)
    _y = -length*np.cos(positions)
    
    #plotTrayectory(_x, _y, "Trayectory of simple pendulum - Semi-implicit Euler method")
    #plotAngles(t, positions, 'Pendulum motion - Semi-implicit Euler method')
    #plotEnergy(t, mass, length, gravity, positions, velocity, "Energy - Semi-implicit Euler method")
    

    #show_animation(length, t, _x, _y, "Motion of simple pendulum - Euler")
    return t, positions, _x, _y, velocity

def euler_explicito(theta0, omega0, t0, tn, h, gravity, length, mass):
    y0 = np.array([theta0, omega0])
    t = np.arange(t0, tn+h, h)
    y = np.zeros((len(t), 2)) # vamos a tener un arreglo de len(t) arreglos de dos elementos
    y[0] = y0

    for i in range(len(t)-1):  
        y[i+1] = y[i] + h*f(t[i], y[i], gravity, length)

    #plot solution

    theta = y[:,0] 
    _x = length*np.sin(theta)
    _y = -1*length*np.cos(theta)
    
    omega = np.zeros(len(t))
    for i in range(len(t)):
        omega[i] = y[i,1]
    
    #plotTrayectory(_x, _y, 'Trayectory of pendulum - Explicit Euler method')
    #plotEnergy(t, mass, length, gravity, theta, omega, 'Energy - Explicit Euler method')
    #plotAngles(t, theta, 'Pendulum motion - Explicit Euler method')
    
    return t, theta, _x, _y, omega
    
    
'''
# RUNGE KUTTA 4
g = 9.81 #gravitational acceleration (m/s^2)
L = 9.8 #length of the pendulum (m)
theta0 = math.pi/3  #initial angle (radians)
omega0 = 0.0  #initial angular velocity (radians)
m = 1.0 #mass of the pendulum
t0, tn = 0.0, 40.0
h = 0.1
'''


def f(ti, yi, gravity, length):
    theta, omega = yi
    theta_dt = omega
    #omega_dt = -(g/L)*np.sin(theta)
    #omega_dt = z_prime(theta, math.sqrt(g/L))
    omega_dt = -(gravity/length)*np.sin(theta) 
    return np.array([theta_dt, omega_dt])

def rungeKutta4_method(gravity, length, theta0, omega0, mass, t0, tn, h):
    y0 = np.array([theta0, omega0])
    
    t = np.arange(t0, tn+h, h)
    y = np.zeros((len(t), 2)) 
    y[0] = y0
    
    for i in range(len(t)-1):
        k1 = h*f(t[i], y[i], gravity, length)
        k2 = h*f(t[i] + 0.5*h, y[i] + 0.5*k1, gravity, length)
        k3 = h*f(t[i] + 0.5*h, y[i] + 0.5*k2, gravity, length)
        k4 = h*f(t[i] + h, y[i] + k3, gravity, length)
        
        y[i+1] = y[i] + (1/6)*(k1 + 2*k2 + 2*k3 + k4)
    
    theta = y[:,0] 
    _x = length*np.sin(theta)
    _y = -length*np.cos(theta)
    
    omega = np.zeros(len(t))
    for i in range(len(t)):
        omega[i] = y[i,1]
    
    #plotTrayectory(_x, _y, 'Trayectory of pendulum - RK4 method')
    #plotEnergy(t, mass, length, gravity, theta, omega, 'Energy - RK4 method')
    #plotAngles(t, theta, 'Pendulum motion - RK4 method')
    
    #show_animation(length, t, _x, _y, "Motion of simple pendulum - RK4")
    
    return t, theta, _x, _y, omega

def jacobiana (theta, omega, gravity, length):

    # Definimos las variables simbólicas
    x1, x2 = sp.symbols('x1 x2')

    # Definimos las funciones simbólicas
    f1 = x2                              # z = theta'
    f2 = -(gravity/length)*sp.sin(x1)    # z' = theta''

    # Calculamos las derivadas parciales
    df1_dx1 = sp.diff(f1, x1)
    df1_dx2 = sp.diff(f1, x2)
    df2_dx1 = sp.diff(f2, x1)
    df2_dx2 = sp.diff(f2, x2)

    # Construimos la matriz Jacobiana general
    J = sp.Matrix([[df1_dx1, df1_dx2], [df2_dx1, df2_dx2]])
    print('Matriz Jacobiana:')
    print(J)

    # Definimos el punto en el que queremos evaluar la matriz Jacobiana
    x1_0 = omega
    x2_0 = theta

    # Sustituimos los valores del punto en las derivadas parciales
    df1_dx1_val = df1_dx1.subs([(x1, x1_0), (x2, x2_0)])
    df1_dx2_val = df1_dx2.subs([(x1, x1_0), (x2, x2_0)])
    df2_dx1_val = df2_dx1.subs([(x1, x1_0), (x2, x2_0)])
    df2_dx2_val = df2_dx2.subs([(x1, x1_0), (x2, x2_0)])

    # Construimos la matriz Jacobiana
    J = np.array([[df1_dx1_val, df1_dx2_val], [df2_dx1_val, df2_dx2_val]])

    # Mostramos la matriz Jacobiana
    print('Matriz Jacobiana:')
    print(J)

    return J

def estabilidad (J, theta, omega):
    # Calculamos los autovalores de la matriz Jacobiana
    eigenvalues = np.linalg.eigvals(J)

    # Mostramos los autovalores
    print('Autovalores:')
    print(eigenvalues)

    # Determinamos la estabilidad del punto de equilibrio
    if np.all(np.real(eigenvalues) < 0):
        print('El punto de equilibrio es estable.')
    elif np.all(np.real(eigenvalues) > 0):
        print('El punto de equilibrio es inestable.')
    else:
        print('Se necesitan técnicas adicionales para determinar la estabilidad del punto de equilibrio.')
        
def plot_trayectories(x_ee, x_ei, x_rk, y_ee, y_ei, y_rk, title):
    plt.plot(x_ee, y_ee, label = 'Explicit Euler method')
    plt.plot(x_ei, y_ei, label = 'Semi-implicit Euler method')
    plt.plot(x_rk, y_rk, label = 'RK4 method')
    plt.title(title)
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)")
    plt.legend()
    plt.show()

def plotAllAngles(t, theta_ee, label1, theta_ei, label2, theta_rk, label3, title):
    plt.plot(t, theta_ee, label = label1)
    plt.plot(t, theta_ei, label = label2)
    plt.plot(t, theta_rk, label = label3)
    plt.title(title)
    plt.xlabel("Time (s)")
    plt.ylabel("Angle (rad)")
    plt.legend()
    plt.show()
    
def plotAllTotalEnergies (t, mass, length, gravity, theta1, omega1, label1, theta2, omega2, label2, theta3, omega3, label3,  title):
    plt.plot(t, E(mass, length, gravity, theta1, omega1), label = label1)
    plt.plot(t, E(mass, length, gravity, theta2, omega2), label = label2)
    plt.plot(t, E(mass, length, gravity, theta3, omega3), label = label3)
    plt.xlabel("Time (s)")
    plt.ylabel("Energy")
    plt.title(title)
    plt.legend()
    plt.show()
    
    
if __name__ == '__main__':
    #t_ei, theta_ei, x_ei, y_ei = euler_semi_implicito(math.pi/3, 0, 0, 10, 0.01, 9.81, 1, 10)
    #t_rk, theta_rk, x_rk, y_rk = rungeKutta4_method(9.81, 1, math.pi/3, 0, 10, 0, 10, 0.01)
    #t_ee, theta_ee1, x_ee1, y_ee1, omega_ee1 = euler_explicito(math.pi/3, 0, 0, 10, 0.01, 9.81, 1, 10)

    """
    print("PRIMER PUNTO")
    j1 = jacobiana(0, (3/2)*math.pi, 9.81, 9.8) # punto de equilibrio (0,0)
    print("SEGUNDO PUNTO")
    j2 = jacobiana(0, math.pi, 9.81, 9.8) # punto de equilibrio (0, pi)
    """