from numpy import array
p = array([125., 125., 125., 125., 125., 125., 125., 125., 124., 121., 118., 112., 106, 98.3, 88.3, 78.5, 68.3, 58.0, 47.6, 37.6, 28.3, 20.0, 12.9, 7.6, 4.9, 3.7, 3.1, 2.7, 2.3, 2.2, 2.0, 1.8, 1.7, 1.6, 1.5, 1.4, 1.4])
powerRatio = p/360.
DCoffsetZ = array([3600., 3400., 3200., 3000., 2800., 2600., 2400., 2200., 2000., 1800., 1600., 1400., 1200., 1000., 800., 600., 400., 200., 0., -200., -400., -600., -800., -1000., -1200., -1400., -1600., -1800., -2000., -2200., -2400., -2600., -2800., -3000., -3200., -3400., -3600.,])

def getPowerDict(inputPower):
    """powerDict(offset) returns the power given by the offset"""
    power = inputPower*powerRatio
    return dict(zip(DCoffsetZ, power))