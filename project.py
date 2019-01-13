from vpython import*
import numpy as np
from numpy import array
import time
t_range = (0,1)

f = open("image0.txt" , "r")
tem = eval(f.read())
f.close()

width = tem.shape
size=vec(width[0],width[1],1)

size_l=(int(size.x),int(size.y),int(size.z))
scene = canvas(width=800, height=800, background=vector(1,1,1),align = 'left',center=size/2)
class boxx(box):
	def __init__(self,pos):
		size=vec(1,1,1)
		#opacity=0.5
		color=vec(0,0.5,1)
		box.__init__(self,pos=pos,size=size,color=color)

boxxx = []
def create(s):
	for x in range(int(s.x)):
		for y in range(int(s.y)):
				boxxx.append(boxx(vec(x,y,1)))

create(size)
t = 0
while True:
	rate(1)
	time.sleep(10)
	T = t*100
	f = open("image%d.txt" %T , "r")
	tem = eval(f.read())
	f.close()
	i = 0 
	for x in range(width[0]):
		for y in range(width[1]):
				boxxx[i].color = vec(tem[x,y],0.5,1-tem[x,y])
				i+=1
	print(t)
	t += 1

