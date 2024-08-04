import json
import numpy as np
import pyvista as pv

# Source: https://github.com/pyvista/pyvista/discussions/3211#discussioncomment-3476299
def reorder_edges(edge_data):
    edges = edge_data.lines.reshape(-1, 3)[:, 1:]  # n_edges x 2

    u, v = edges.T

    adj = np.empty(u.shape, dtype=u.dtype)
    adj[u] = v

    v = edges[0, 0]  # Starting vertex
    vert_idxs = [v]
    for _ in range(len(edges)):
        v = adj[v]
        vert_idxs.append(v)

    new_lines = np.empty((edges.shape[0], 3), dtype=edges.dtype)
    new_lines[:, 0] = 2
    new_lines[:, 1] = range(edges.shape[0])
    new_lines[:, 2] = range(1, edges.shape[0] + 1)
    return pv.PolyData(edge_data.points[vert_idxs], lines=new_lines)

data = json.load(open(f'formantData.json'))
formantData = np.array(data['formantData'])
f1, f2, f3, f4 = formantData[0::4], formantData[1::4], formantData[2::4], formantData[3::4]

f1Min, f1Max = min(f1), max(f1)
f2Min, f2Max = min(f2), max(f2)
f3Min, f3Max = min(f3), max(f3)

cloud = pv.PolyData(np.c_[np.log2(f1), np.log2(f2), np.zeros(len(f1))])
surf = cloud.delaunay_2d(alpha=1/8)
boundary = reorder_edges(surf.extract_feature_edges(boundary_edges=True, non_manifold_edges=False, manifold_edges=False))

boundary_points = []
for pt in boundary.points:
  for i in range(0, len(f1)):
    if np.log2(f1[i]) == pt[0] and np.log2(f2[i]) == pt[1]:
      boundary_points.append(i)

print(boundary_points)
with open(f'formantDataBoundary.js', 'w') as f:
  f.write(f'var formantDataBoundary = {boundary_points};\n')
  f.write(f'var minF1 = {f1Min};\n')
  f.write(f'var maxF1 = {f1Max};\n')
  f.write(f'var minF2 = {f2Min};\n')
  f.write(f'var maxF2 = {f2Max};\n')
  f.write(f'var minF3 = {f3Min};\n')
  f.write(f'var maxF3 = {f3Max};\n')

p = pv.Plotter()
p.add_mesh(cloud, point_size=10, render_points_as_spheres=True)
p.add_mesh(surf, opacity=0.5)
p.add_mesh(boundary, scalars=range(boundary.n_points), render_lines_as_tubes=True, line_width=5, show_scalar_bar=False)
p.show()



