import generateBlock as GB
import numpy as np

################################
# For block B
l = 6  # y/2
w = 4  # x/2
alpha = 2/180*3.1415926
h = w*np.tan(alpha)  # z
vertexes = []
vertexes.append([-w, l/2, h])
vertexes.append([-w, -l/2, h])
vertexes.append([0, -l/2, 0])
vertexes.append([0, l/2, 0])
vertexes.append([-w, l/2, 0])
vertexes.append([-w, -l/2, 0])
edges = [[0, 1], [1, 2], [2, 3], [3, 0], [
    3, 4], [4, 5], [5, 2], [1, 5], [0, 4]]
faces = [[0, 1, 2, 3], [0, 3, 4], [1, 5, 2], [2, 5, 4, 3], [0, 4, 5, 1]]

block1 = GB.Block(vertexes, edges, faces, fixed=1)

#######################################
# For block A
l = 4  # y/2
w = 3  # x/2
alpha = 2/180*3.1415926
h = w*np.tan(alpha)  # z
vertexes = []
vertexes.append([0, -l/2, 0])
vertexes.append([0, l/2, 0])
vertexes.append([w, l/2, 0])
vertexes.append([w, -l/2, 0])
vertexes.append([w, -l/2, -h])
vertexes.append([w, l/2, -h])
edges = [[0, 1], [1, 2], [2, 3], [3, 0], [
    3, 4], [4, 5], [5, 2], [1, 5], [0, 4]]
faces = [[0, 3, 2, 1], [0, 4, 3], [1, 2, 5], [2, 3, 4, 5], [0, 1, 5, 4]]

offset = 1
print(offset*np.tan(alpha))
GB.move(vertexes, [-offset, 0, offset*np.tan(alpha)])

block2 = GB.Block(vertexes, edges, faces)

GB.genBlk([block1, block2], "smallAngle.dda")
