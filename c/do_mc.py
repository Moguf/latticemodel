#!/usr/bin/env python
import subprocess
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def ploter(seed):
    cmdline="./a.out 100000000 " + str(seed)
    output=""
    output=subprocess.check_output(cmdline,shell=True)

    hxdata=[]
    hydata=[]

    pxdata=[]
    pydata=[]

    xdata=[]
    ydata=[]

    energy=0
    time=0
    for iline in output.split("\n"):
        ilist=iline.strip().split()
        if ilist==[]:
            continue
        if ilist[0]==">>":
            energy=ilist[-1]
        if ilist[0]==">":
            time = ilist[-1]
        if ilist[0]==">>>":
            print "%3s %3s %3s %3s %3s"%(ilist[1],ilist[2],ilist[3],ilist[4],ilist[5])
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
    lim=13
    ax.set_xlim([-lim,lim])
    ax.set_ylim([-lim,lim])
    ax.set_aspect('equal')
    ax.scatter(hxdata,hydata,color="g",alpha=0.5,s=200)
    ax.scatter(pxdata,pydata,color="b",alpha=0.5,s=200)
    ax.plot(xdata,ydata,lw=2,alpha=0.3)
    ax.set_title("Energy="+energy+":time="+time)
    plt.savefig("test.png")

    
ploter(3)
