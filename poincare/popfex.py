"""
popfex : Poincare Plot Feature Extaction Library
License : GPL3 
programmed by Birol Kuyumcu    
"""

import numpy as np
import matplotlib.pyplot as plt

def showTimeSeries(ts):
	x = np.arange(len(ts))
	plt.plot(x,ts)
	return

def showSMatrix(ts):
	angle,mag = PoincareVect(ts)
	img = aMatrix(angle,delay=1)	
	plt.imshow(img)
	return

def timeS2Poincare(ts,ndelay=1):
    x = ts[0:-2*ndelay]
    y = ts[ndelay:-ndelay]
    return x,y

def showPoincare(ts):
    x,y = timeS2Poincare(ts)
    plt.plot(x,y,'r.')
    return

def PoincareVect(ts,ndelay=1):
    x1 = ts[0:-2*ndelay]
    x2 = y1 = ts[ndelay:-ndelay]
    y2 = ts[2*ndelay:]
    deltax = x2-x1
    deltay = y2-y1
    angle = np.arctan2(deltax, deltay) * 180 / np.pi    
    angle = np.where(angle < 0, angle+360,angle)
    angle = angle / 18
    mag = (deltax**2 + deltay**2)**0.5    
    return angle,mag

def showPoincareVect(ts):
	angle,mag = PoincareVect(ts)
	sm = sMatrix(angle,mag)
	x_pos = np.arange(20)
	plt.bar(x_pos, sm)

	return

def showPoincareVect2(ts):
    angle,mag = PoincareVect(ts)
    plt.plot(mag,'r')

    return

def sMatrix(angle,mag):
    smat = np.zeros(20)
    for j in range (angle.shape[0]):
        smat[int(angle[j])] += mag[j]
        
    return (smat / smat.max())

def aMatrix(angle,delay=1):
    amat = np.zeros((20,20))
    for i in range(angle.shape[0]-delay):
        amat[int(angle[i]),int(angle[i+delay])] +=1
    return amat / amat.max()