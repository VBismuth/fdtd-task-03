# -*- encoding: utf-8 -*-

from numpy import pi, sqrt, exp, zeros
import numpy as np
from scipy.fft import fft, fftfreq
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

def DispInProc(E: List[float], H: List[float],
               dx: float, Xmax: int,
               p_pos: int, s_pos: int, dur):

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)
    x = np.linspace(0.0, Xmax*dx, Xmax)

    ax = axs[0]
    ax.set_ylim([min(E)-0.2, max(E)+0.2])
    ax.set_ylabel('Ez, В/м')
    ax.grid()
    
    
    ax.plot(x, E)
    ax.plot(s_pos*dx, 0, 'ok', p_pos*dx, 0, 'rx')

    ax = axs[1]
    ax.set_xlim([0, Xmax*dx])
    ax.set_ylim([min(H)-0.02, max(H)+0.02])
    ax.set_xlabel('x, м')
    ax.set_ylabel('Hy, А/м')
    ax.grid()
    
    ax.plot(x, H)
    ax.plot(s_pos*dx, 0, 'ok', p_pos*dx, 0, 'rx')
    ax.plot(s_pos, 0, 'ok', p_pos, 0, 'rx')
    fig.suptitle("Распостранение ЭМ волны" + \
                 " через\n "+ str(dur*1e6) + \
                 " мкс от начала моделирования")
    plt.show()

def DispGraphs(probe: Probe, s_pos: int, dt: float, maxTime):
    E = probe.E
    H = probe.H
    p_pos = probe.pos
    t = np.linspace(0.0, len(E)*dt, len(E))
    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)

    ax = axs[0]
    ax.set_ylim([min(E)-0.2, max(E)+0.2])
    ax.set_ylabel('Ez, В/м')
    ax.grid()
    
    
    ax.plot(t*1e6, E)

    ax = axs[1]
    ax.set_xlim([0, max(t)*1e6])
    ax.set_ylim([min(H)-0.002, max(H)+0.002])
    ax.set_xlabel('q, мкс')
    ax.set_ylabel('Hy, А/м')
    ax.grid()
    
    ax.plot(t*1e6, H)
    fig.suptitle("Временной сигнал,"+ \
    " зарегистрированный датчиком")
    plt.show()


def DispSpectrum(probe: Probe, dt, N):
    E = probe.E
    H = probe.H
    
    Ef = fft(E)
    E_max = np.abs(max(Ef))
    Hf = fft(H)

    H_max = np.abs(max(Hf))

    xf = fftfreq(N, dt)

    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True)

    ax = axs[0]
    ax.set_ylim([-0.2, 1.2])
    ax.set_ylabel('$|{E_m}/{E_{m max}}|$')
    ax.grid()
    ax.plot(xf/1e9, np.abs(Ef/E_max))

    ax = axs[1]
    ax.set_ylim([-0.2, 1.2])
    ax.set_xlim([0, max(xf)/1e9])
    ax.set_ylabel('$|{H_m}/{H_{m max}}|$')
    ax.set_xlabel('f, ГГц')
    ax.grid() 
    ax.plot(xf/1e9, np.abs(Hf/H_max))
    

    fig.suptitle("Нормированный амплитудный спектр\n" \
    "сигнала, зарегистрированный датчиком")
    plt.show()
