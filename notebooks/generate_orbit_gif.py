from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import sys
import imageio
from config import *

orbit = [(xn,0)]

for i in range(n):
    
    orbit.append((xn, xn))
    orbit.append((xn, F(xn)))
    xn = F(xn)
    
    if abs(xn) > 100:
        break
    
orbit_df = pd.DataFrame(orbit, columns= ["x", "fx"])

mnx, mxx = min(-0.5, orbit_df.x.min()),  max(0.5, orbit_df.x.max())
mny, mxy = min(-0.5, orbit_df.fx.min()), max(0.5, orbit_df.fx.max())

rx = mxx - mnx
ry = mxy - mny

mnx -= (0.1 *rx)
mxx += (0.1 *ry)
mny -= (0.1 *rx)
mxy += (0.1 *ry)


def plot_orbit(F, orbit_df):
    
    x = np.arange(mnx, mxx, (mxx - mnx) / 30)
    y = np.apply_along_axis(F, 0, x)
    y2 = x

    fig, ax = plt.subplots(figsize=(8, 8), dpi=80)

    # x axis
    ax.plot([mnx - (0.1 *rx), mxx + (0.1 *rx)], [0,0], color="#00000033", linewidth=1, linestyle='dashed')
    #y axis
    ax.plot([0,0], [mny - (0.1 *rx), mxy + (0.1 *ry)], color="#00000033", linewidth=1, linestyle='dashed')

    # function
    ax.plot(x, y)
    # identity
    ax.plot(x, y2)
    
    ax.scatter(orbit_df.x[0], orbit_df.fx[0], color='g')

    ax.plot(orbit_df.x[:-1], orbit_df.fx[:-1], linewidth=2, linestyle=':', color='g')
    ax.plot(orbit_df.x[-2:], orbit_df.fx[-2:], linewidth=2, linestyle='-', color='g')

    ax.set_xlim(mnx, mxx )
    ax.set_ylim(mny, mxy )
    
    # Used to return the plot as an image rray
    fig.canvas.draw()       # draw the canvas, cache the renderer
    image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
    image  = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

    return image


if __name__ == '__main__':
    
    imageio.mimsave('./%s.gif'% (sys.argv[1]), [plot_orbit(F, orbit_df[:i]) for i in range(2,orbit_df.shape[0])], fps=3)