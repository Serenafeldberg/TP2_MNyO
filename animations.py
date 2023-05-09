import pendulo_doble as pd

#PENDULO DOBLE

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

    pd.anim(t, theta1, theta2, length, length, 'Péndulo doble theta1 = pi/6, theta2 = pi/4')

    #Angulos iniciales pi/6 y pi/3.9
    theta1_0 = pd.math.pi/6
    theta2_0 = pd.math.pi/3.9

    theta1, theta2, theta1_p, theta2_p = pd.solve_pendulum(theta1_0, theta2_0, theta1_p0, theta2_p0, delta, t0, tf, N, length, gravity)

    pd.anim(t, theta1, theta2, length, length, 'Péndulo doble theta1 = pi/6, theta2 = pi/3.9')

    #Angulos iniciales pi y pi/1.1 -> angulos grandes
    theta1_0 = pd.math.pi
    theta2_0 = pd.math.pi/1.1

    theta1, theta2, theta1_p, theta2_p = pd.solve_pendulum(theta1_0, theta2_0, theta1_p0, theta2_p0, delta, t0, tf, N, length, gravity)

    pd.anim(t, theta1, theta2, length, length, 'Péndulo doble theta1 = pi, theta2 = pi/1.1')

    #Angulos iniciales pi/2.25 y pi/2.25 -> angulos medianos
    theta1_0 = pd.math.pi/2.25
    theta2_0 = pd.math.pi/2.25

    theta1, theta2, theta1_p, theta2_p = pd.solve_pendulum(theta1_0, theta2_0, theta1_p0, theta2_p0, delta, t0, tf, N, length, gravity)

    pd.anim(t, theta1, theta2, length, length, 'Péndulo doble theta1 = pi/2.25, theta2 = pi/2.25')

    #Angulos iniciales pi/30 y pi/30 -> angulos pequeños
    theta1_0 = pd.math.pi/30
    theta2_0 = pd.math.pi/30

    theta1, theta2, theta1_p, theta2_p = pd.solve_pendulum(theta1_0, theta2_0, theta1_p0, theta2_p0, delta, t0, tf, N, length, gravity)

    pd.anim(t, theta1, theta2, length, length, 'Péndulo doble theta1 = pi/30, theta2 = pi/30')




if __name__ == "__main__":
    t = pd.np.linspace(0, 20, 2000)
    animations_pd1(t)