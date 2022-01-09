import numpy as np
import random as rd
import time
import os
from tkinter.messagebox import *

nstep=0
iout=0
TTT=0.0

var_space_group = ''
var_KMC_element1 = ''
var_KMC_element2 = ''
var_KMC_step = ''
var_KMC_gap = ''
var_KMC_temperature = ''

def get_sg(sg):
    global var_space_group
    var_space_group = sg

def get_element1(e1):
    global var_KMC_element1
    var_KMC_element1 = e1

def get_element1(e2):
    global var_KMC_element2
    var_KMC_element2 = e2

def get_step(step):
    global var_KMC_step
    var_KMC_step = step

def get_gap(gap):
    global var_KMC_gap
    var_KMC_gap = gap

def get_temperature(t):
    global var_KMC_temperature
    var_KMC_temperature = t
##############################################

def KMC_calculation_():
     global nstep
     global iout
     global TTT
     def make_bulk(dx, dy, dz, cl):
         nli=int(dx/cl)
         nlj=int(dy/cl)
         nlz=int(dz/cl)
         xyz=[]
         nbulk=0
         for i in range(-nli, nli):
             for j in range(-nlj, nlj):
                 for k in range(-nlz, nlz):
                     xyz.append([i*cl, j*cl, k*cl])
                     xyz.append([(0.5+i)*cl, (0.5+j)*cl, k*cl])
                     xyz.append([i*cl, (0.5+j)*cl, (0.5+k)*cl])
                     xyz.append([(0.5+i)*cl, j*cl, (0.5+k)*cl])
                     nbulk=nbulk+4
         return nbulk,xyz

     def readstru():
         with open('ini_stru.xyz') as f:
             lines=f.readlines()
             npar=int(lines[0])
             dx,dy,dz,cl=list(map(float,lines[1].split()))
             nli=int(dx/cl)
             nlj=int(dy/cl)
             nlz=int(dz/cl)
             xyz=[]
             nbulk=0
             for i in range(-nli,nli):
                 for j in range(-nlj,nlj):
                     for k in range(-nlz,nlz):
                         xyz.append([i*cl,j*cl,k*cl])
                         xyz.append([(0.5+i)*cl,(0.5+j)*cl,k*cl])
                         xyz.append([i*cl,(0.5+j)*cl,(0.5+k)*cl])
                         xyz.append([(0.5+i)*cl,j*cl,(0.5+k)*cl])
                         nbulk=nbulk+4        
             elem=[]
             nelem=0 
             nt=[100]*nbulk
             minz=float(lines[2].split()[3])
             for i in range(npar):           
                 el=lines[2+i].split()[0]
                 if (all(elem[nnn]!=el for nnn in range(nelem))):
                     elem.append(el)
                     nelem=nelem+1
                 x=float(lines[2+i].split()[1])
                 y=float(lines[2+i].split()[2])
                 z=float(lines[2+i].split()[3])
                 if (z<minz):
                     minz=z
                 for j in range(nbulk):
                     dx=xyz[j][0]-x
                     dy=xyz[j][1]-y
                     dz=xyz[j][2]-z
                     if (abs(dx)<0.001 and abs(dy)<0.001 and abs(dz)<0.001):
                         for k in range(nelem):
                             if (el==elem[k]):
                                 nt[j]=k
         return (nelem,elem,npar,nt,minz,nbulk,xyz)

     def NN_cal(nbulk,ntype,xyz,rcut,nelem):
         a=np.array(xyz)
         cn=np.array([0]*nbulk)
         effcn=np.array([[0]*nelem]*nbulk)
         nn=np.array([[0]*12]*nbulk)
         for i in range(nbulk):
             for j in range(i+1,nbulk):
                 r=np.sum((a[i]-a[j])**2)
                 r=np.sqrt(r)
                 if (r<rcut):
                     cn[i]=cn[i]+1
                     cn[j]=cn[j]+1
                     nn[i][cn[i]-1]=j
                     nn[j][cn[j]-1]=i
                     if (ntype[i]<nelem):
                         effcn[j][ntype[i]]=effcn[j][ntype[i]]+1
                     if (ntype[j]<nelem):
                         effcn[i][ntype[j]]=effcn[i][ntype[j]]+1
         return cn,effcn,nn
                 
     def Hij(ai,cn,ecn,nn,nt,nelem,elem):
         NNbond={'Au-Au':-0.2,'Pd-Pd':-0.36,'Pt-Pt':-0.40,
                 'Au-Pd':-0.304,'Pd-Au':-0.304,
                 'Pd-Pt':-0.38,'Pt-Pd':-0.38}
         Eai=0.0
         Eaj=np.zeros(12)
         eli=elem[nt[ai]]
         for iel in range(nelem):
             el=eli+'-'+elem[iel]
             Ebond=ecn[ai][iel]*NNbond[el]
             Eai=Eai+Ebond
         for icn in range(cn[ai]):
             aj=nn[ai,icn]
             if (nt[aj]>nelem):
                 for iel in range(nelem):
                     el=eli+'-'+elem[iel]
                     Ebond=ecn[aj][iel]*NNbond[el]
                     Eaj[icn]=Eaj[icn]+Ebond
             Eaj[icn]=Eaj[icn]-NNbond[eli+'-'+eli]
         return (Eai,Eaj)

     def Hev(ai,xyz,minz,cn,ecn,nn,nt,nelem,Eai,Eaj,a,b):
         if (xyz[ai][2]==minz):
             si=sum(ecn[ai])
             Eev=si*a+b
             Eai=Eai+Eev
         for icn in range(cn[ai]):
             aj=nn[ai,icn]
             if (xyz[aj][2]==minz):
                 if (nt[aj]>nelem):
                     sj=sum(ecn[aj])-1
                     Eev=sj*a+b
                     Eaj[icn]=Eaj[icn]+Eev
         return (Eai,Eaj)

     def ratei(ai,cn,ecn,nn,nt,nelem,Eai,Eaj,T,Ebase):
         kb=8.62e-5
         rij=np.zeros(12)
         ri=0.0
         si=sum(ecn[ai])
         for icn in range(cn[ai]):
             aj=nn[ai,icn]
             if (nt[aj]>nelem):
                 sj=sum(ecn[aj])-1
                 if (si*sj!=0.0):
                     Ediff=Eaj[icn]-Eai-Eaj[icn]**2/(Eai+Eaj[icn])
                 if (si!=0.0 and sj==0.0):
                     Ediff=Eaj[icn]-Eai-Ebase-Eaj[icn]**2/(Eai+Ebase+Eaj[icn])
                 if (si==0.0 and sj!=0.0):
                     Ediff=Eaj[icn]+Ebase-Eai-(Eaj[icn]+Ebase)**2/(Eai+Ebase+Eaj[icn])
                 if (si+sj==0.0):
                     Ediff=0.0
                 rij[icn]=np.exp(-Ediff/kb/T)
                 if (rij[icn]>1.0):
                     rij[icn]=1.0
                 ri=ri+rij[icn]
         return (rij,ri)

     def pickij(rall,cycle,ri):
         rpoint=rd.uniform(0.0,rall)
         rsum=0.0
         for i in range(cycle):
             rsum=rsum+ri[i]
             if (rsum>=rpoint):
                 break
         return (i)

    #  main body    
    ## with open('KMC_input') as f1:
    # lines=f1.readlines()
     time_start=time.time()
     key_sp=''
     nstep=int(var_KMC_step )
     iout=int(var_KMC_gap )
     TTT=float(var_KMC_temperature )
    ##     nstep,iout=int(lines[1].split()[0]),int(lines[1].split()[1])
    ##     TTT=float(lines[2])
    ##     print nstep
    ##     print iout
    ##     print TTT
     Ebase,Epa,Epb=0,0,0

     rcut=3.0

     nelem,elem,npar,ntype,minz,nbulk,xyz=readstru()
     print(minz,Epa,Epb)
     print(nelem,elem)
     print(nbulk)

     nsp=0
     xyzsp=[]
     if (key_sp=='sp'):
         for i in range(nbulk):
             if (xyz[i][2]<minz):
                 nsp=nsp+1
                 xyzsp.append(xyz[i])
                 ntype[i]=nelem       
     cn,effcn,nnsite=NN_cal(nbulk,ntype,xyz,rcut,nelem)
     print(npar)

     Ei=np.zeros(nbulk)
     Ej=np.zeros((nbulk,12))
     rij=np.zeros((nbulk,12))
     rsite=np.zeros(nbulk)
     for ai in range(nbulk):
         if (ntype[ai]<nelem and sum(effcn[ai])<12):
             Ei[ai],Ej[ai]=Hij(ai,cn,effcn,nnsite,ntype,nelem,elem)
             if (key_sp=='sp'):
                 Ei[ai],Ej[ai]=Hev(ai,xyz,minz,cn,effcn,nnsite,ntype,nelem,Ei[ai],Ej[ai],Epa,Epb)       
             rij[ai],rsite[ai]=ratei(ai,cn,effcn,nnsite,ntype,nelem,Ei[ai],Ej[ai],TTT,Ebase) 
     rtot=np.sum(rsite)
     print(rtot)
     f=open('cluster.xyz','w')
     for istep in range(nstep):
         if ((istep%iout)==0):
             f.write(str(npar+nsp)+'\n')
             f.write(str(istep)+'\n') 
             for n in range(nelem):
                 for i in range(nbulk):
                     if (ntype[i]==n):
                         el=elem[ntype[i]]
                         f.write(el+'  '+str(xyz[i][0])+' '+str(xyz[i][1])+' '+str(xyz[i][2])+'\n')    
             for i in range(nsp):
                 f.write('H'+'  '+str(xyzsp[i][0])+'  '+str(xyzsp[i][1])+'  '+str(xyzsp[i][2])+'\n')

         picki=pickij(rtot,nbulk,rsite)
         pickicn=pickij(rsite[picki],12,rij[picki])
         pickj=nnsite[picki,pickicn]

         temp_ntype=ntype[picki]
         ntype[picki]=ntype[pickj]
         ntype[pickj]=temp_ntype
         for icn in range(cn[picki]):
             aj=nnsite[picki,icn]
             effcn[aj][temp_ntype]=effcn[aj][temp_ntype]-1
         for jcn in range(cn[pickj]):
             ak=nnsite[pickj,jcn]
             effcn[ak][temp_ntype]=effcn[ak][temp_ntype]+1

         rij[picki]=0.0
         rsite[picki]=0.0

         for icn in range(cn[picki]):
             aj=nnsite[picki,icn]
             if (ntype[aj]<nelem):
                 Ei[aj],Ej[aj]=Hij(aj,cn,effcn,nnsite,ntype,nelem,elem)
                 if (key_sp=='sp'):
                     Ei[ai],Ej[ai]=Hev(ai,xyz,minz,cn,effcn,nnsite,ntype,nelem,Ei[ai],Ej[ai],Epa,Epb)       
                 rij[aj],rsite[aj]=ratei(aj,cn,effcn,nnsite,ntype,nelem,Ei[aj],Ej[aj],TTT,Ebase)            

         for jcn in range(cn[pickj]):
             ak=nnsite[pickj,jcn]
             if (ntype[ak]<nelem):
                 Ei[ak],Ej[ak]=Hij(ak,cn,effcn,nnsite,ntype,nelem,elem)
                 if (key_sp=='sp'):
                     Ei[ai],Ej[ai]=Hev(ai,xyz,minz,cn,effcn,nnsite,ntype,nelem,Ei[ai],Ej[ai],Epa,Epb)       
                 rij[ak],rsite[ak]=ratei(ak,cn,effcn,nnsite,ntype,nelem,Ei[ak],Ej[ak],TTT,Ebase)
         rtot=np.sum(rsite)

     f.write(str(npar+nsp)+'\n')
     f.write(str(istep)+'\n')
     for n in range(nelem):
         for i in range(nbulk):
             if (ntype[i]==n):
                 el=elem[ntype[i]]
                 f.write(el+'  '+str(xyz[i][0])+' '+str(xyz[i][1])+' '+str(xyz[i][2])+'\n')
     for i in range(nsp):
         f.write('H'+'  '+str(xyzsp[i][0])+'  '+str(xyzsp[i][1])+'  '+str(xyzsp[i][2])+'\n')
     f.close()
     print(nsp)
     with open('last_one.xyz','w') as f:
         f.write(str(npar+nsp)+'\n')
         f.write(str(istep)+'\n')
         for n in range(nelem):
             for i in range(nbulk):
                 if (ntype[i]==n):
                     el=elem[ntype[i]]
                     f.write(el+'  '+str(xyz[i][0])+' '+str(xyz[i][1])+' '+str(xyz[i][2])+'\n')
         for i in range(nsp):
             f.write('H'+'  '+str(xyzsp[i][0])+'  '+str(xyzsp[i][1])+'  '+str(xyzsp[i][2])+'\n')
     time_end=time.time()
     total_time=round(time_end-time_start , 4)
   #  import tkMessageBox as mBox
     showinfo('Calculation Message Info Box','Job Completed.'+'\n'+'Total Cost About: '+str(total_time)+' Seconds')