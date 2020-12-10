# -*- encoding: utf-8 -*-

import numpy as np
import toolbox as tlbx

if __name__ == "__main__":
    
    # Параметры по заданию Вар15
    X_size = 7.5   # Размер области моделирования, м
    dx = 1e-3      # Размер ячейки разбиения, м
    Df = 0.5e+9    # Ширина спектра сигнала, Гц

    sig = 10.0     # Ширина спектра сигнала в отсчетах
    dt = 1/sig/Df  # Дискрет по времени, с
   
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * np.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета
    maxTime = 5000

    # Размер области моделирования в отсчетах
    maxSize = int(X_size // dx)

    # Положение источника в отсчетах
    sourcePos = 1000

    # Датчики для регистрации поля в отсчетах
    p_pos = 3000
    probe = tlbx.Probe(p_pos, maxTime)

    Ez = np.zeros(maxSize)
    Hy = np.zeros(maxSize)

    source = tlbx.GaussianPlaneWave(120.0, sig, Sc)

    probe.addData(Ez, Hy)

    for q in range(1, maxTime):
        # Граничные условия для поля H
        Hy[-1] = Hy[-2]

        # Расчет компоненты поля H
        Ez_shift = Ez[1:]
        Hy[:-1] = Hy[:-1] + (Ez_shift - Ez[:-1]) * Sc / W0

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Hy[sourcePos - 1] -= (Sc / W0) * source.getE(sourcePos, q)
        # Hy[sourcePos - 1] -= (Sc / W0) * source.getE(0, q)
        
        # Граничные условия для поля E
        Ez[0] = Ez[1]

        # Расчет компоненты поля E
        Hy_shift = Hy[:-1]
        Ez[1:] = Ez[1:] + (Hy[1:] - Hy_shift) * Sc * W0

        # Источник возбуждения с использованием метода
        # Total Field / Scattered Field
        Ez[sourcePos] += Sc * source.getE(sourcePos - 0.5, q + 0.5)
        # Ez[sourcePos] += Sc * source.getE(-0.5, q + 0.5)
        
        # Регистрация поля в датчиках
        probe.addData(Ez, Hy)
        
        if q == int(maxTime/5 * 4):
            dur = q*dt
            tlbx.DispInProc(Ez, Hy, dx, maxSize, p_pos, sourcePos, dur)
   
    tlbx.DispGraphs(probe, sourcePos, dt, maxTime)
    tlbx.DispSpectrum(probe, dt, maxTime)
    
