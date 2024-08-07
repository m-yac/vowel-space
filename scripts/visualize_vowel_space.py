import json
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.widgets import Slider
import numpy as np

data = json.load(open(f'formantData.json'))
formantData = np.array(data['formantData'])
f1, f2, f3, f4 = formantData[0::4], formantData[1::4], formantData[2::4], formantData[3::4]

f1Min, f1Max = min(f1), max(f1)
f2Min, f2Max = min(f2), max(f2)
f3Min, f3Max = min(f3), max(f3)
f4Min, f4Max = min(f4), max(f4)

fig = plt.figure(figsize=(12, 12))

ax = fig.add_subplot(projection='3d')
ax.set_xlabel('F1')
ax.set_ylabel('F2')
ax.set_zlabel('F3')

def log_tick_formatter(val, pos=None):
    return f'{round(2**val)}'
ax.xaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
ax.yaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))
ax.zaxis.set_major_formatter(mticker.FuncFormatter(log_tick_formatter))

c_by_boundary = [(1.0 if i == 0 or i == 15 else 0.2,
                  1.0 if j == 0 or j == 15 else 0.2,
                  1.0 if k == 0 or k == 15 else 0.2)
                  for i in range(0,16) for j in range(0,16) for k in range(0,16)]

pts = ax.scatter(np.log2(f1), np.log2(f2), np.log2(f3), c=np.log2(f4), label='F4')

fig.colorbar(pts).set_label("F4")

grid_x, grid_y, grid_z = np.meshgrid(np.linspace(np.log2(f3Min), np.log2(f3Max), 16), np.linspace(np.log2(f2Min), np.log2(f2Max), 16), np.linspace(np.log2(f1Max), np.log2(f1Min), 16))
grid_x, grid_y, grid_z = grid_x.flatten(), grid_y.flatten(), grid_z.flatten()

axSlider = fig.add_axes([0.2, 0.1, 0.65, 0.03])   
slider = Slider(
    ax=axSlider,
    label='Projection',
    valmin=0,
    valmax=1,
    valinit=1,
)

def update(val):
    x = (1 - val) * grid_z + val * np.log2(f1)
    y = (1 - val) * grid_y + val * np.log2(f2)
    z = (1 - val) * grid_x + val * np.log2(f3)
    pts._offsets3d = (x, y, z)
    fig.canvas.draw()
slider.on_changed(update)

fig.tight_layout()
plt.subplots_adjust(bottom=0.2)

plt.show()

