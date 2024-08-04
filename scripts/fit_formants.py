import json
from matplotlib.widgets import Slider
import numpy as np 
from matplotlib import pyplot as plt
import scipy
from scipy.signal import find_peaks

def norm(array):
  return array / max(abs(array))

data = json.load(open('plot_data.txt'));
peaks = data["peaks"]
spectrum = np.array(data["spectrum"])
dspectrum = np.gradient(spectrum)
ddspectrum = np.gradient(dspectrum)
# dddspectrum = np.gradient(ddspectrum)
# dspectrum = np.array([spectrum[i+1] - spectrum[i] for i in range(0,len(spectrum)-1)])
# ddspectrum = np.array([dspectrum[i+1] - dspectrum[i] for i in range(0,len(spectrum)-2)])
# dddspectrum = np.array([ddspectrum[i+1] - ddspectrum[i] for i in range(0,len(spectrum)-3)])

fig, ax = plt.subplots()
plt.axhline()
plt.vlines(peaks, -1, 1, color="black")
plt.plot(np.arange(0,len(spectrum)), norm(spectrum), color="blue")
plt.plot(np.arange(0,len(spectrum)), 2 * norm(dspectrum), color="green")
plt.plot(np.arange(0,len(spectrum)), 4 * norm(ddspectrum), color="red")
# plt.plot(np.arange(0,len(spectrum)), 8 * norm(dddspectrum), color="purple")
# plt.plot(np.arange(0.5, len(spectrum) - 0.5), 2 * norm(dspectrum), color="green")
# plt.plot(np.arange(1.0, len(spectrum) - 1.0), 4 * norm(ddspectrum), color="red")
# plt.plot(np.arange(1.5, len(spectrum) - 1.5), 8 * norm(dddspectrum), color="purple")
plt.show()
exit()


r = round(peaks[0])

def get_points(peak_idx):
  peak = peaks[peak_idx]
  i = round(peak)
  xf = peak - i
  yf = norm(np.array([spectrum[i-2], spectrum[i-1], spectrum[i], spectrum[i+1], spectrum[i+2]]))
  dyf = norm(np.array([yf[1] - yf[0], yf[2] - yf[1], yf[3] - yf[2], yf[4] - yf[3]]))
  ddyf = norm(np.array([dyf[1] - dyf[0], dyf[2] - dyf[1], dyf[3] - dyf[2]]))
  dddyf = norm(np.array([ddyf[1] - ddyf[0], ddyf[2] - ddyf[1]]))
  x = xf + r
  y = norm(spectrum[i - r : i + r])
  dy = norm(dspectrum[i - r : i + r])
  ddy = norm(ddspectrum[i - r : i + r])
  dddy = norm(dddspectrum[i - r : i + r])
  return xf, yf, dyf, ddyf, dddyf, x, y, dy, ddy, dddy

fig, (ax1, ax2, ax3) = plt.subplots(1, 3, width_ratios=[8,1,8])
fig.tight_layout()
ax1.plot(np.array([-2, 2]), [0,0], color="black")
ax3.plot(np.array([0, 2*r-1]), [0,0], color="black")

xf0, yf0, dyf0, ddyf0, dddyf0, x0, y0, dy0, ddy0, dddy0 = get_points(0)
lxf, = ax1.plot([xf0, xf0], [-1, 1], color = 'black')
lyf, = ax1.plot(np.array([-2, -1, 0, 1, 2]), yf0, marker='o', color="blue")
ldyf, = ax1.plot(np.array([-1.5, -0.5, 0.5, 1.5]), dyf0, marker='o', color="green")
lddyf, = ax1.plot(np.array([-1, 0, 1]), ddyf0, marker='o', color="red")
ldddyf, = ax1.plot(np.array([-0.5, 0.5]), dddyf0, marker='o', color="purple")
lx, = ax3.plot([x0, x0], [0, 1], color="black")
ly, = ax3.plot(np.arange(0, 2 * r), y0, color="blue")
ldy, = ax3.plot(np.arange(0.5, 2 * r + 0.5), dy0, color="green")
lddy, = ax3.plot(np.arange(1.0, 2 * r + 1.0), ddy0, color="red")
ldddy, = ax3.plot(np.arange(1.5, 2 * r + 1.5), dddy0, color="purple")

slider = Slider(
  ax=ax2,
  label='Peak',
  valmin=0,
  valmax=len(data["peaks"]),
  valstep=1,
  valinit=0,
  orientation="vertical"
)

def update(val):
  xf, yf, dyf, ddyf, dddyf, x, y, dy, ddy, dddy = get_points(slider.val)
  lxf.set_xdata([xf, xf])
  lyf.set_ydata(yf)
  ldyf.set_ydata(dyf)
  lddyf.set_ydata(ddyf)
  ldddyf.set_ydata(dddyf)
  lx.set_xdata([x, x])
  ly.set_ydata(y)
  ldy.set_ydata(dy)
  lddy.set_ydata(ddy)
  ldddy.set_ydata(dddy)
  fig.canvas.draw_idle()
slider.on_changed(update)

plt.show()



exit()

M = 4096
F = 11
C = (np.pi ** 2) / 3

y = np.array(data["spectrum"])
x = np.arange(0, len(y))

peaks = find_peaks(y)[0]

def f(args = np.ones(4*F)):
  tot = np.zeros(x.shape)
  for i in range(0, F):
    a, b, c, d = args[i], args[F+1], args[F+F+i], args[F+F+F+i]
    a *= peaks[i]
    b *= 2 * F / M
    c *= y[peaks[i]] - 0.08
    d *= 0.125 - 0.08
    x1 = C * (b * (x - a)) ** 2
    tot += c * (1 + x1) / (1 + (3*c/d - 2) * x1)
  return tot

def lsq(args):
  return np.sum((y - f(args))**2)

# guesses = np.ones(4*F)
# res = scipy.optimize.minimize(lsq, guesses)
# print(res)

plt.plot(x, y)
# plt.plot(x, f())
# plt.plot(x, f(res.x))
# plt.plot(x, (y - f(res.x))**2)
plt.show()
