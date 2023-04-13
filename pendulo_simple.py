import math


def E (t, v):
    return t+v

def T (m, l, theta):
    return (1/2)*m*(math.pow(l,2))*(math.pow(theta, 2))

def V (m, g, l, theta):
    return m*g*l - (m*g*l)*(math.cos(theta))