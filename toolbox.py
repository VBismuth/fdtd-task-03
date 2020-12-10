# -*- encoding: utf-8 -*-

from numpy import pi, sqrt, exp, zeros, fft
from typing import List
from matplotlib import pyplot as plt

class GaussianPlaneWave:

    def __init__(self, delay, width, Sc=1.0, eps=1.0, mu=1.0):
        self.d = delay
        self.w = width
        self.Sc = Sc
        self.eps = eps
        self.mu = mu

    def getE(self, m, q):
        return exp(-(((q - m * sqrt(self.eps * self.mu) / \
        self.Sc) - self.d) / self.w) ** 2)

class Probe:
    def __init__(self, pos: int, T_max: int):
        self.pos = pos
        self.max = T_max
        
        self.E = zeros(T_max)
        self.H = zeros(T_max)
        
        self.time = 0

    def addData(self, E: List[float], H: List[float]):

        if self.time<self.max:
            self.E[self.time] = E[self.pos]
            self.H[self.time] = H[self.pos]
            self.time += 1
        else:
            print("Can't add more data")

def DispInProc(E: List[float],
            H: List[float], 
            p_pos: int, s_pos: int):

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)

    ax = axs[0]
    ax.set_xlim([0, len(E)])
    ax.set_ylim([min(E)-0.2, max(E)+0.2])
    ax.set_ylabel('Ez, В/м')
    ax.grid()
    
    
    ax.plot(E)
    ax.plot(s_pos, 0, 'ok', p_pos, 0, 'rx')

    ax = axs[1]
    ax.set_xlim([0, len(H)])
    ax.set_ylim([min(H)-0.02, max(H)+0.02])
    ax.set_xlabel('x, отсчет')
    ax.set_ylabel('Hy, А/м')
    ax.grid()
    
    ax.plot(H)
    ax.plot(s_pos, 0, 'ok', p_pos, 0, 'rx')
    fig.subtitle("Процесс распостранения \
                  электромагнитной волны")
    plt.show()

def DispGraphs(probe: Probe, s_pos: int):
    E = probe.E
    H = probe.H
    p_pos = probe.pos
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)

    ax = axs[0]
    ax.set_ylim([min(E)-0.2, max(E)+0.2])
    ax.set_ylabel('Ez, В/м')
    ax.grid()
    
    
    ax.plot(E)

    ax = axs[1]
    ax.set_xlim([0, len(H)])
    ax.set_ylim([min(H)-0.002, max(H)+0.002])
    ax.set_xlabel('q, отсчет')
    ax.set_ylabel('Hy, А/м')
    ax.grid()
    
    ax.plot(H)
    fig.subtitle("Временной сигнал, \
                  зарегистрированный датчиком")
    plt.show()



