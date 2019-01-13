#原子彈熱傳導
import math
from numpy import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
seterr(all='ignore')

#模擬的最大時間及間隔
maxIter = 80000
tPlot = 100

#畫面大小
ny = 401
ratio = 0.5	
nx =nz = int(ny*ratio)

'''dx = 1/(ny-2)
Pr = 1 #nulb/kapa
Ra = 20000'''
#buoyancy = [0,gr]
gr = 0.0005

#沒用到
cxO = nx//2
cyO = 2*ny//3
r = ny//9

#溫度
Thot = 1#500
Tcold = 0#300
T0 = (Thot+Tcold)/2

#dt = sqrt(gr*dx)
#空氣參數，熱擴散、黏滯係數
nulb = 0.018#dt/(dx**2)*sqrt(Pr/Ra)
kapa = 0.026#dt/(dx**2)/sqrt(Pr*Ra)
omegaN = 1/(3*nulb+0.5)
omegaT1 = 1/(3*kapa+0.5)#1.2#
omegaT2 = omegaT1*0.6 #障礙物的熱擴散係數

#D2Q9, D2Q5
v = array([[1,1],[1,0],[1,-1],[0,1],[0,0],[0,-1],[-1,1],[-1,0],[-1,-1]])
'''v = array([0,0,0],[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1],
	[1,1,0],[-1,1,0],[1,-1,0],[-1,-1,0],
	[1,0,1],[-1,0,1],[1,0,-1],[-1,0,-1],
	[0,1,1],[0,-1,1],[0,1,-1],[0,-1,-1])
weightN = array([1/3,1/18,1/18,1/18,1/18,1/18,1/18,
	1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36,1/36]) '''#collision的權重
weightN = array([1/36,1/9,1/36,1/9,4/9,1/9,1/36,1/9,1/36])
weightT= array([0, 1/6, 0, 1/6,1/3,1/6, 0 , 1/6,0])
col1 = array([0,1,2])
col2 = array([3,4,5])
col3 = array([6,7,8])
colup = array([2,5,8])

#初始物理量
rho = 1
def inivel(d, x, y):
	return 0.0#(1-d)*0.04*(1+1E-4*sin(y/(ny-1)*2*pi))#+1
vel = fromfunction(inivel, (2,nx,ny))
u = zeros((2,nx,ny))#fromfunction(inivel, (2,nx,ny)) #u0
cx, cy = nx//2, ny//3 #熱源中心
T = Tcold*ones((nx,ny))
T[cx-5:cx+5,cy-5:cy+5] = Thot*10
force = zeros((nx,ny))

#設置障礙物
def obstacle_fun(x, y):
	return False#(x-cxO)**2+(y-cyO)**2<r**2
obstacle = zeros((nx,ny),dtype=bool)#fromfunction(obstacle_fun, (nx,ny))
#print(type(obstacle))
boundaryUD = True
boundaryLR = False
if boundaryUD:
	obstacle[:,-1]=0
	obstacle[:,0]=1
if boundaryLR:
	obstacle[0,:]=True
	obstacle[-1,:]=True
omegaT = ones((nx,ny))*omegaT1
omegaT[obstacle] = omegaT2

def equilibrium():
	usqr = 3/2 * (u[0]**2+u[1]**2)
	feq = zeros((9,nx,ny))
	for i in range(9):
		vu = 3 * (v[i,0]*u[0,:,:]+v[i,1]*u[1,:,:])
		feq[i,:,:] = rho*weightN[i] * (1+vu+0.5*vu**2-usqr)
	return feq
	
def equilibriumT():
	#usqr = 3/2 * (u[0]**2+u[1]**2)
	Teq = zeros((9,nx,ny))
	for i in range(9):
		vu = 3 * (v[i,0]*u[0,:,:]+v[i,1]*u[1,:,:])
		Teq[i,:,:] = T*weightT[i]*(1+vu)#+0.5*vu**2-usqr)
	return Teq

#初始溫度場、速度場
fin = equilibrium()
Tin = equilibriumT()

for time in range(maxIter):
	#walls
	fin[col3,-1,:] = fin[col3,-2,:]
	fin[colup,:,-1] = fin[colup,:,-2]
	fin[col1,0,:] = fin[col1,1,:]

	#macro
	rho = sum(fin,axis=0)
	T = sum(Tin, axis=0)
	u = zeros((2,nx,ny))
	for i in range(9):
		u[0,:,:] += v[i,0] * fin[i,:,:]
		u[1,:,:] += v[i,1] * fin[i,:,:]
	u /= rho
	
	#設置熱源、冷源
	T[cx,cy] = Thot*5 
	
	#left wall
	#u[:,0,:] = vel[:,0,:]
	#rho[0,:] = 1/(1-u[0,0,:]) * (sum(fin[col2,0,:],axis=0)+2*sum(fin[col3,0,:],axis=0))
	
	#eq
	feq = equilibrium()
	Teq = equilibriumT()
	#fin[[0,1,2],0,:] = feq[[0,1,2],0,:]+fin[[8,7,6],0,:]-feq[[8,7,6],0,:]
	
	#collision
	force = zeros((9,nx,ny))
	avgT = average(T)
	for i in range(9):
		force[i] = 3*weightN[i]*rho*(T-avgT)*(v[i,0]*0+v[i,1]*gr)/(Thot-Tcold)
		force[i,obstacle] = 0.0
	fout = fin - omegaN * (fin - feq) + force
	
	Tout = Tin-omegaT*(Tin-Teq)
	
	#boundary
	for i in range(9):
		fout[i,obstacle] = fin[8-i,obstacle]
		
	#streaming
	for i in range(9):
		fin[i,:,:] = roll(roll(fout[i,:,:],v[i,0],axis=0),v[i,1],axis=1)
		Tin[i,:,:] = roll(roll(Tout[i,:,:],v[i,0],axis=0),v[i,1],axis=1)
		
	#boundary for T
	heatup = True
	if heatup:
		Tin[5,:,-1]=Tcold - Tin[3,:,-1] - Tin[7,:,-1] - Tin[1,:,-1] - Tin[4,:,-1]
		Tin[3,:,0] =Tcold  - Tin[7,:,0] - Tin[1,:,0] - Tin[5,:,0] - Tin[4,:,0]
	elif True:
		Tin[3,:,0]=Tcold - Tin[5,:,0] - Tin[7,:,0] - Tin[1,:,0] - Tin[4,:,0]
		Tin[5,:,-1] =Thot  - Tin[7,:,-1] - Tin[1,:,-1] - Tin[3,:,-1] - Tin[4,:,-1]
	if True:
		Tin[1,0,:]=Tcold - Tin[5,0,:] - Tin[7,0,:] - Tin[3,0,:] - Tin[4,0,:]
		Tin[7,-1,:] =Tcold  - Tin[1,-1,:] - Tin[5,-1,:] - Tin[3,-1,:] - Tin[4,-1,:]
		'''Tin[:,0,:] =(4*Tin[:,1,:]-Tin[:,2,:])/3
		Tin[:,-1,:]=(4*Tin[:,-2,:]-Tin[:,-3,:])/3'''
	
	#visualize
	if (time%tPlot==0):
		set_printoptions(threshold=nan)
		f = open("image%d.txt" %time ,"w")
		T_T = array_repr(T)
		f.write(str(T_T))
		f.close()
		"""
		print(T.max(), T.min())
		image = T#-(u[1,2:,1:-1]-u[1,:-2,1:-1])/2+(u[0,1:-1,2:]-u[0,1:-1,:-2])/2#u[0]**2+u[1]**2)#
		image = fliplr(image)
		v1max = T[:,cy+3:].max()
		plt.clf()
		plt.imshow(image.transpose(), cmap=cm.jet, vmax=Thot, vmin=Tcold)
		plt.savefig('thermo2.{0:03d}.png'.format(time//tPlot))"""