import pendulo_doble as pd
import pendulo_simple as ps 

def animations_pd1 (t):
    # Definir las condiciones iniciales
    theta1_p0 = 0
    theta2_p0 = 0
    delta = 0
    t0 = 0
    tf = 20
    N = 2000
    length = 1
    gravity = 9.8

    #Angulos iniciales pi/6 y pi/4
    theta1_0 = pd.math.pi/6
    theta2_0 = pd.math.pi/4

    theta1, theta2, theta1_p, theta2_p = pd.solve_pendulum(theta1_0, theta2_0, theta1_p0, theta2_p0, delta, t0, tf, N, length, gravity)

    pd.anim(t, theta1, theta2, length, length, 'Motion of double pendulum - Figure 29')

    #Angulos iniciales pi/6 y pi/3.9
    theta1_0 = pd.math.pi/6
    theta2_0 = pd.math.pi/3.9

    theta1, theta2, theta1_p, theta2_p = pd.solve_pendulum(theta1_0, theta2_0, theta1_p0, theta2_p0, delta, t0, tf, N, length, gravity)

    pd.anim(t, theta1, theta2, length, length, 'Motion of double pendulum - Figure 33')

    #Angulos iniciales pi y pi/1.1 -> angulos grandes
    theta1_0 = pd.math.pi
    theta2_0 = pd.math.pi/1.1

    theta1, theta2, theta1_p, theta2_p = pd.solve_pendulum(theta1_0, theta2_0, theta1_p0, theta2_p0, delta, t0, tf, N, length, gravity)

    pd.anim(t, theta1, theta2, length, length, 'Motion of double pendulum - Figure 37')

    #Angulos iniciales pi/2.25 y pi/2.25 -> angulos medianos
    theta1_0 = pd.math.pi/2.25
    theta2_0 = pd.math.pi/2.25

    theta1, theta2, theta1_p, theta2_p = pd.solve_pendulum(theta1_0, theta2_0, theta1_p0, theta2_p0, delta, t0, tf, N, length, gravity)

    pd.anim(t, theta1, theta2, length, length, 'Motion of double pendulum - Figure 41')

    #Angulos iniciales pi/30 y pi/30 -> angulos pequeños
    theta1_0 = pd.math.pi/30
    theta2_0 = pd.math.pi/30

    theta1, theta2, theta1_p, theta2_p = pd.solve_pendulum(theta1_0, theta2_0, theta1_p0, theta2_p0, delta, t0, tf, N, length, gravity)

    pd.anim(t, theta1, theta2, length, length, 'Motion of double pendulum - Figure 45')


def animations_ps():
    mass = 10
    length = 1
    gravity = 9.81
    t0, tn = 0, 10
    h = 0.01
    theta0, omega0 = 0, ps.math.pi/3
    
    #Figura 2
    ps.euler_explicito(theta0, omega0, t0, tn, h, gravity, length, mass, ps.f, "Motion of simple pendulum - Figure 2")
    ps.euler_semi_implicito(theta0, omega0, t0, tn, h, gravity, length, mass, "Motion of simple pendulum - Figure 5")
    
    

if __name__ == "__main__":
    #t = pd.np.linspace(0, 20, 2000)
    #animations_pd1(t)
    animations_ps()
