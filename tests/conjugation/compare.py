import numpy as np
import matplotlib.pyplot as plt

N = 10
L = 100.0
dq = 2*np.pi/L
dx = L/N


r0 = np.array([dx*4,-2*dx,13*dx])

qxs = np.linspace(0,N/2,num=N//2+1,endpoint=True)*dq

qfwd = np.linspace(0,N/2,num=N//2+1,endpoint=True)*dq
qbwd = np.linspace(-N/2+1,0,num=N//2-1,endpoint=False)*dq

qys = np.concatenate((qfwd,qbwd)) 
qzs = np.concatenate((qfwd,qbwd))

real = np.empty([N,N,N//2+1],float)+100
imag = np.empty([N,N,N//2+1],float)+100


realfac = np.sqrt(2)
imagfac = 0

for k,qz in enumerate(qzs):
    for j,qy in enumerate(qys):
        for i,qx in enumerate(qxs):

            tmp = r0[0]*qx+r0[1]*qy+r0[2]*qz
            
            real[k,j,i] = np.cos(tmp)
            imag[k,j,i] = np.sin(tmp)


            if (i == 0 or i == N//2):
                if (j == 0 or j == N//2):
                    if (k == 0 or k == N//2):
                        real[k,j,i] *= realfac



real_c = real*0
imag_c = imag*0

for i in range(4):
    chunk = np.loadtxt(f"real_{i}.txt")
    if ( i != 3):
        real_c[i*3:(i+1)*3,:,:] = chunk.reshape((3,N,N//2+1))
    else:
        real_c[i*3:,:,:] = chunk.reshape((1,N,N//2+1))


        
    chunk = np.loadtxt(f"imag_{i}.txt")
    if ( i != 3):
        imag_c[i*3:(i+1)*3,:,:] = chunk.reshape((3,N,N//2+1))
    else:
        imag_c[i*3:,:,:] = chunk.reshape((1,N,N//2+1))



print(real[~np.isclose(real,real_c)])
print(imag[~np.isclose(imag,imag_c)])




"""
for out in [imag,imag_c]:

    squarechunk0 = out[:N//2+1,:N//2+1,0]
    squarechunk1 = out[:N//2+1,N//2+1:,0]
    squarechunk2 = out[N//2+1:,:N//2+1,0]
    squarechunk3 = out[N//2+1:,N//2+1:,0]


    mat = np.zeros([N,N],float)+10000

    mat[:N//2-1,:N//2-1] = squarechunk3
    mat[:N//2-1,N//2-1:] = squarechunk2
    mat[N//2-1:,:N//2-1] = squarechunk1
    mat[N//2-1:,N//2-1:] = squarechunk0

    plt.contourf(mat)
    plt.colorbar()
    plt.show()

"""
