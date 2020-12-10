# -*- encoding: utf-8 -*-

import numpy as np
import toolbox as tlbx

if __name__ == "__main__":
    
    # Волновое сопротивление свободного пространства
    W0 = 120.0 * np.pi

    # Число Куранта
    Sc = 1.0

    # Время расчета в отсчетах
    maxTime = 400

    # Размер области моделирования в отсчетах
    maxSize = 200

    # Положение источника в отсчетах
    sourcePos = 50

    # Датчики для регистрации поля
    p_pos = 75
    probe = tlbx.Probe(p_pos, maxTime)

    Ez = np.zeros(maxSize)
    Hy = np.zeros(maxSize)

    source = tlbx.GaussianPlaneWave(30.0, 10.0, Sc)

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
        
        if q == 120:
            tlbx.DispInProc(Ez, Hy, p_pos, sourcePos)
   
    tlbx.DispGraphs(probe, sourcePos)
    
