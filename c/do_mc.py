#!/usr/bin/env python
import subprocess

import matplotlib.pyplot as plt

cmdline="./a.out 1000"
output=subprocess.check_output(cmdline,shell=True)

hxdata=[]
hydata=[]

pxdata=[]
pydata=[]

xdata=[]
ydata=[]

energy=0

for iline in output.split("\n"):
    ilist=iline.strip().split()
    if ilist==[]:
        continue
    if ilist[0]==">>":
        energy=ilist[-1]
    if ilist[0]==">>>":
        xdata.append(int(ilist[-2]))
        ydata.append(int(ilist[-1]))
        if ilist[1]=="1":
            hxdata.append(int(ilist[-2]))
            hydata.append(int(ilist[-1]))
        else:
            pxdata.append(int(ilist[-2]))
            pydata.append(int(ilist[-1]))

            
fig=plt.figure()
ax=fig.add_subplot(111)
lim=10
ax.set_xlim([-lim,lim])
ax.set_ylim([-lim,lim])
ax.scatter(hxdata,hydata,color="g",alpha=0.5,s=300)
ax.scatter(pxdata,pydata,color="b",alpha=0.5,s=300)
ax.plot(xdata,ydata,lw=2,alpha=0.3)
ax.set_title("Energy="+energy)
plt.show()

