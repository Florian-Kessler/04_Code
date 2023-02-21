M = np.array([[-7.93273539e-01, -1.31847220e-02,  6.08722643e-01,  1.13000000e+03],
              [ 2.21123987e-02,  9.98482113e-01,  5.04431538e-02,  4.53000000e+02],
              [-6.08463749e-01,  5.34755369e-02, -7.91777893e-01,  7.64000000e+02],
              [ 0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  1.00000000e+00]])
[t1, t2, t3] = rotation_angles_from_matrix(M[:3, :3], 'zyx')
M2 = rotation_matrix_from_angles(t1, t2+3.14159, t3, 'zyx')
print(round(t1*180/np.pi))
print(round(t2*180/np.pi))
print(round(t3*180/np.pi))
print(M2)
origin = [0, 0, 0]
ori = M[:3, 3]
center = [5.75, 8.75, 0]
v1 = M[0, :3]
v2 = M[1, :3]
v3 = M[2, :3]
p1P = [54, 474, 820]
plt.figure()
ax = plt.axes(projection='3d')
fact = 100
for i in range(-0, 5):
    ax.scatter3D(ori[0]+i*fact*v1[0], ori[1]+i*fact*v1[1], ori[2]+i*fact*v1[2], c='r', alpha=1)
    ax.scatter3D(ori[0]+i*fact*v2[0], ori[1]+i*fact*v2[1], ori[2]+i*fact*v2[2], c='g', alpha=1)
    ax.scatter3D(ori[0]+i*fact*v3[0], ori[1]+i*fact*v3[1], ori[2]+i*fact*v3[2], c='b', alpha=1)

    ax.scatter3D(origin[0]+i*fact*1, origin[0]+i*fact*0, origin[0]+i*fact*0, c='r')
    ax.scatter3D(origin[0]+i*fact*0, origin[0]+i*fact*1, origin[0]+i*fact*0, c='g')
    ax.scatter3D(origin[0]+i*fact*0, origin[0]+i*fact*0, origin[0]+i*fact*1, c='b')

ax.scatter3D(p1P[0], p1P[1], p1P[2], c='c')
ax.scatter3D([0, 1500], [0, 1500], [0, 1500], alpha=0)
ax.scatter3D(686, 480, 177, c='c', alpha=1)  # Ti tip
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()