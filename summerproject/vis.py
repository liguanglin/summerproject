import os
import numpy as np
from open3d import *
import re

f=open('points.txt','r')
points = re.split(' |\n',f.read())[:-1]
points = list(map(int,points))
points = np.asarray(points)
points = points.reshape((-1,3))
colors = np.zeros(points.shape)
point_cloud = PointCloud()
point_cloud.points = Vector3dVector(points)
point_cloud.colors = Vector3dVector(colors)
draw_geometries([point_cloud])
