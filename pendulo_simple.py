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

def euler_semi_implicito (y0, x0, t0, tn, step, gravity, length, mass,ani_title): 
    w0 = math.sqrt(gravity/length)
    t_start = t0
    t_end = tn
    t = np.arange(t_start, t_end+step, step)
    positions = np.zeros(len(t)) 
    velocity = np.zeros(len(t)) 
    positions[0] = y0 
    velocity[0] = x0 

    for i in range(1, len(t)):
        velocity[i] = velocity[i-1] + step*z_prime(positions[i-1], w0)
        positions[i] = positions[i-1] + step*z(velocity[i])
    
    _x = length*np.sin(positions)
    _y = -length*np.cos(positions)
    
    #plotTrayectory(_x, _y, "Trayectory of simple pendulum - Semi-implicit Euler method")
    #plotAngles(t, positions, 'Pendulum motion - Semi-implicit Euler method')
    #plotEnergy(t, mass, length, gravity, positions, velocity, "Energy - Semi-implicit Euler method")

    show_animation(length, t, _x, _y, ani_title)
    return t, positions, _x, _y, velocity

def euler_explicito(theta0, omega0, t0, tn, h, gravity, length, mass, f, ani_title):
    y0 = np.array([theta0, omega0])
    t = np.arange(t0, tn+h, h)
    y = np.zeros((len(t), 2)) 
    y[0] = y0

    for i in range(len(t)-1):  
        y[i+1] = y[i] + h*f(t[i], y[i], gravity, length)

    theta = y[:,0] 
    _x = length*np.sin(theta)
    _y = -1*length*np.cos(theta)
    
    omega = np.zeros(len(t))
    for i in range(len(t)):
        omega[i] = y[i,1]
    
    #plotTrayectory(_x, _y, 'Trayectory of pendulum - Explicit Euler method')
    #plotEnergy(t, mass, length, gravity, theta, omega, 'Energy - Explicit Euler method')
    #plotAngles(t, theta, 'Pendulum motion - Explicit Euler method')
    
    show_animation(length, t, _x, _y, ani_title)
    
    return t, theta, _x, _y, omega
    


def f(ti, yi, gravity, length):
    theta, omega = yi
    theta_dt = omega
    omega_dt = -(gravity/length)*np.sin(theta) 
    return np.array([theta_dt, omega_dt])

def rungeKutta4_method(gravity, length, theta0, omega0, mass, t0, tn, h, f):
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
        
def plot_trayectories(x_ee, x_ee_l, y_ee, y_ee_l, title):
    plt.plot(x_ee, y_ee, label = 'Original ODE')
    plt.plot(x_ee_l, y_ee_l, label = 'Linearized ODE')
    plt.title(title)
    plt.xlabel("Horizontal position (m)")
    plt.ylabel("Vertical position (m)")
    plt.legend()
    plt.show()

def plotAllAngles(t, theta_ee, label1, theta_ee_l, label2, title):
    plt.plot(t, theta_ee, label = label1)
    plt.plot(t, theta_ee_l, label = label2)
    plt.title(title)
    plt.xlabel("Time (s)")
    plt.ylabel("Angle (rad)")
    plt.legend()
    plt.show()
    
def plotAllTotalEnergies (t1, t2, t3, mass, length, gravity, theta1, omega1, label1, theta2, omega2, label2, theta3, omega3, label3, title):
    plt.plot(t1, E(mass, length, gravity, theta1, omega1), label = label1)
    plt.plot(t2, E(mass, length, gravity, theta2, omega2), label = label2)
    plt.plot(t3, E(mass, length, gravity, theta3, omega3), label = label3)
    plt.xlabel("Time (s)")
    plt.ylabel("Energy")
    plt.title(title)
    plt.legend()
    plt.show()

def linearized_f(ti, yi, gravity, length):
    theta, omega = yi
    theta_dt = omega
    omega_dt = -(gravity/length)*theta 
    return np.array([theta_dt, omega_dt])

def relativeError_totalEnergy(mass, length, gravity, theta1, omega1, label1, theta2, omega2, label2, theta3, omega3, label3, title, t1, t2, t3):
    
    energy1 = E(mass, length, gravity, theta1, omega1)
    energy2 = E(mass, length, gravity, theta2, omega2)
    energy3 = E(mass, length, gravity, theta3, omega3)
    
    ground_truth1 = energy1[0]
    ground_truth2 = energy2[0]
    ground_truth3 = energy3[0]
    
    error1 = np.zeros(len(energy1))
    for i in range(len(energy1)):
        error1[i] = abs(energy1[i]-ground_truth1)/ground_truth1 
    
    error2 = np.zeros(len(energy2))
    for i in range(len(energy2)):
        error2[i] = abs(energy2[i]-ground_truth2)/ground_truth2
    
    error3 = np.zeros(len(energy3))
    for i in range(len(energy3)):
        error3[i] = abs(energy3[i]-ground_truth3)/ground_truth3 
    
    
    plt.plot(t1, error1, label = label1)
    plt.plot(t2, error2, label = label2)
    plt.plot(t3, error3, label = label3)
    plt.xlabel("Time (s)")
    plt.ylabel("Error")
    plt.title(title)
    plt.legend()
    plt.show()
    
    
if __name__ == '__main__':
    #t_ei, theta_ei, x_ei, y_ei, omega_ei = euler_semi_implicito(math.pi/3, 0, 0, 50, 0.01, 9.81, 1, 10)
    t_rk, theta_rk, x_rk, y_rk, omega_rk = rungeKutta4_method(9.81, 1, math.pi/3, 0, 10, 0, 10, 0.01, f)
    t_rk_l, theta_rk_l, x_rk_l, y_rk_l, omega_rk_l = rungeKutta4_method(9.81, 1, math.pi/3, 0, 10, 0, 10, 0.01, linearized_f)
    #t_ee, theta_ee, x_ee, y_ee, omega_ee = euler_explicito(math.pi/3, 0, 0, 10, 0.01, 9.81, 1, 10, f)
    #t_ee_l, theta_ee_l, x_ee_l, y_ee_l, omega_ee_l = euler_explicito(math.pi/3, 0, 0, 10, 0.01, 9.81, 1, 10, linearized_f)
    
    energia_rk4 = E(1, 1, 9.81, theta_rk, omega_rk)
    energia_rk4l = E(1, 1, 9.81, theta_rk_l, omega_rk_l)
    

    '''
    print("Energia inicial: ", energia_rk4[0])
    print("Energia final: ", energia_rk4[-1])
    print("Energia media: ", np.mean(energia_rk4))
    print("desviacion estandar: ", np.std(energia_rk4))

    print("Energia inicial linealizada: ", energia_rk4l[0])
    print("Energia final linealizada: ", energia_rk4l[-1])
    print("Energia media linealizada: ", np.mean(energia_rk4l))
    print("desviacion estandar linealizada: ", np.std(energia_rk4l))
    '''

    """
    print("PRIMER PUNTO")
    j1 = jacobiana(0, (3/2)*math.pi, 9.81, 9.8) # punto de equilibrio (0,0)
    print("SEGUNDO PUNTO")
    j2 = jacobiana(0, math.pi, 9.81, 9.8) # punto de equilibrio (0, pi)
    """