import numpy as np
import time
import os
from functools import reduce
from scipy.optimize import fsolve
from tkinter.messagebox import *
##def DJF_update():
##    if var_support =='Cu111-ZnO0001':
##         element_chosen.current(1)
##         environment1_chosen.current(6)
##         element_chosen.config(state=DISABLED)
##         environment1_chosen.config(state=DISABLED)
##         element_chosen.after(1000,DJF_update)
##         environment1_chosen.after(1000,DJF_update)
##    else:
##         element_chosen.config(state=NORMAL)
##         environment1_chosen.config(state=NORMAL)       
##
##         element_chosen.after(1000,DJF_update)
##         environment1_chosen.after(1000,DJF_update)
         
#_thread.start_new_thread( DJF_update,() )

##t1 = threading.Thread(target=DJF_update)
##t1.start()
##t1.join()
#DJF_update()

########################################################
NPs = ''
latt_para = 0
A_atom = np.array([])
gamma = np.array([])
Gas = np.array([])
GasPP = np.array([])
E_ads = np.array([[]])
w = np.array([[[]]])
theta = np.array([[]])
pressure = 0
temperature = 0
d = 0
Flag = False
R = 8.314
unit_coversion = 0.0000103643
k_b = 0.000086173303
bond_length = 3.0
filename_xyz=''
filename_record = 'record.txt'  # filename for record coverage file
#############################################################
var_support = ''
var_element = ''
var_radius = ''
var_pressure = ''
var_temperature = ''
var_environment1 = ''
var_environment2 = ''
var_environment3 = ''
var_environment4 = ''
var_environment5 = ''
var_PP1 = ''
var_PP2 = ''
var_PP3 = ''
var_PP4 = ''
var_PP5 = ''

def get_support(support):
    global var_support
    var_support = support

def get_element(element):
    global var_element
    var_element = element

def get_radiu(radius):
    global var_radius
    var_radius = radius

def get_pressure(pressure):
    global var_pressure
    var_pressure = pressure

def get_temperature(t):
    global var_temperature
    var_temperature = t

def get_environment1(e1):
    global var_environment1
    var_environment1 = e1

def get_environment2(e2):
    global var_environment2
    var_environment2 = e2

def get_environment3(e3):
    global var_environment3
    var_environment3 = e3

def get_environment4(e4):
    global var_environment4
    var_environment4 = e4

def get_environment5(e5):
    global var_environment5
    var_environment5 = e5

def get_PP1(pp1):
    global var_PP1
    var_PP1 = pp1

def get_PP2(pp2):
    global var_PP2
    var_PP2 = pp2

def get_PP3(pp3):
    global var_PP3
    var_PP3 = pp3

def get_PP4(pp4):
    global var_PP4
    var_PP4 = pp4
    
def get_PP5(pp5):
    global var_PP5
    var_PP5 = pp5
################################################################

def geometry_1():  
     if (var_support == 'Cu111-ZnO000-1' and (var_element != 'None' or var_environment1 != 'H2O-gas')) :
          showinfo('Support','The only supported nanoparticle in this version is: \nCu(111)-ZnO(000-1)-H2O-gas\nIf you want to calculate the structure of supported metal NP, Please choose :\n Element : None\n Environment1 : H2O-gas')
          return  
     if ((var_element != 'Pd' ) & (var_element != 'Pt' ) & (var_element  != 'Cu' ) & (var_element != 'Au' ) & (var_element != 'Rh' )) and (var_support =='None' ):
          showinfo('Element','Please choose the given metal element as a selection')
          return
     if (var_radius =='') or (float(var_radius )<=15):
          showinfo('Radius','Please enter the radius greater than 15 angstrom')
          return
     if var_pressure =='':
          showinfo('Pressure','Please enter the pressure')
          return
     if var_temperature =='':
          showinfo('Pressure','Please enter the temperature')
          return
     if (var_environment1 !='' and var_environment1 !='None') and (var_PP1 =='' and
                                                                             (float(var_PP1 )<=0 or float(var_PP1 )>100)):
          showinfo('Partial Pressure 1','Please enter the partial pressure 1 (from 0-100) of environment1')
          return
     if (var_environment2 !='' and var_environment2 !='None') and (var_PP2 =='' and
                                                                             (float(var_PP2 )<=0 or float(var_PP2 )>100)):
          showinfo('Partial Pressure 2','Please enter the partial pressure 2 (from 0-100) of environment2')
          return
     if (var_environment3 !='' and var_environment3 !='None') and (var_PP3 =='' and
                                                                             (float(var_PP3 )<=0 or float(var_PP3 )>100)):
          showinfo('Partial Pressure 3','Please enter the partial pressure 3 (from 0-100) of environment3')
          return
     if (var_environment4 !='' and var_environment4 !='None') and (var_PP4 =='' and
                                                                             (float(var_PP4 )<=0 or float(var_PP4 )>100)):
          showinfo('Partial Pressure 4','Please enter the partial pressure 4 (from 0-100) of environment4')
          return
     if (var_environment5 !='' and var_environment5 !='None') and (var_PP5 =='' and
                                                                             (float(var_PP5 )<=0 or float(var_PP5 )>100)):
          showinfo('Partial Pressure 5','Please enter the partial pressure 5 (from 0-100) of environment5')
          return
     # Lattice parameter
     L={'Pd':3.891,'Pt':3.924,'Cu':3.615,'Au':4.078,'Rh':3.804}
     #Area per surface atom[100,101,111]
     G={'Pd':[7.569,10.704,6.555],'Pt':[7.698,10.887,6.667],'Cu':[6.553,9.239,5.658],
                  'Au':[8.316,11.761,7.202],'Rh':[7.237,10.234,6.627]}
     #surface tension
     H={'Pd':[0.094,0.099,0.081],'Pt':[0.111,0.117,0.088],'Cu':[0.091,0.095,0.081],
                  'Au':[0.050,0.057,0.040],'Rh':[0.148,0.146,0.125]}
     # Adsorption energy[100,101,111]
     M={'Pd-O2':[-0.96,-0.95,-0.39],'Pd-H2':[-0.35,-0.43,-0.25],'Pd-N2':[-0.460,-0.654,-0.233],
        'Pd-CO2':[-0.47,-0.54,-0.02], 'Pd-NO2':[-1.14,-1.33,-0.88],'Pd-CO':[-1.564,-1.611,-1.591],
        'Pd-NO':[-1.75,-1.90,-1.79],
        'Pt-O2':[-0.83,-1.19,-0.21],'Pt-H2':[-0.69,-0.62,-0.48],'Pt-N2':[-0.24,-0.42,-0.01],
        'Pt-CO2':[0,0,-0.02],'Pt-NO2':[-1.35,-1.65,-0.88],'Pt-CO':[-1.84,-1.91,-1.48],
        'Pt-NO':[-1.92,-2.07,-1.47],     
        'Cu-O2':[-1.28,-1.02,-0.16],'Cu-H2':[-0.02,0,-0.01],'Cu-N2':[-0.02,0,-0.01],
        'Cu-CO2':[-0.01,-0.01,-0.03],'Cu-NO2':[-1.22,-1.42,-0.89],'Cu-CO':[-0.58,-0.64,-0.47],
        'Cu-NO':[-0.99,-0.95,-0.70],
        'Au-O2':[-0.22,-0.34,0],'Au-H2':[-0.09,-0.04,0],'Au-N2':[-0.08,-0.03,0],
        'Au-CO2':[-0.10,-0.06,-0.01],'Au-NO2':[-0.78,-0.85,-0.50],'Au-CO':[-0.60,-0.62,-0.22],
        'Au-NO':[-0.62,-0.59,-0.31],
        'Rh-CO':[-1.638,-1.724,-1.624],'Rh-NO':[-2.17,-2.26,-2.11],'Rh-O2':[0,0,0],
        'Cu-H2O-gas':[-0.40,-0.47,-0.29],"Pd-H2O-gas":[-0.39,-0.50,-0.39],"Pt-H2O-gas":[-0.47,-0.56,-0.38],"Au-H2O-gas":[-0.3,-0.3,-0.27],
        'Cu-H2O-liquid':[-0.40,-0.47,-0.29],"Pd-H2O-liquid":[-0.39,-0.50,-0.39],"Pt-H2O-liquid":[-0.47,-0.56,-0.38],"Au-H2O-liquid":[-0.3,-0.3,-0.27]}
     # Interaction [100,101,111]
     N={'Pd-O2':[-0.33,-0.32,-0.21],'Pd-H2':[-0.06,-0.08,-0.06],'Pd-N2':[-0.189,-0.222,-0.121],
        'Pd-CO2':[-0.31,-0.44,-0.19], 'Pd-NO2':[-0.45,-0.53,-0.31],
        'Pt-O2':[-0.36,-0.40,-0.15],'Pt-H2':[-0.01,0,-0.02],'Pt-N2':[-0.18,-0.18,-0.13],
        'Pt-CO2':[-0.17,-0.16,-0.17],'Pt-NO2':[-0.51,-0.59,-0.30],
        'Cu-O2':[-0.42,-0.54,-0.10],'Cu-H2':[-0.03,-0.03,-0.03],'Cu-N2':[-0.26,-0.23,-0.24],
        'Cu-CO2':[-0.36,-0.35,-0.33],'Cu-NO2':[-0.63,-0.82,-0.47],
        'Au-O2':[0,0,0],'Au-H2':[-0.02,-0.02,0],'Au-N2':[-0.09,-0.08,-0.07],
        'Au-CO2':[-0.10,-0.08,-0.08],'Au-NO2':[-0.24,-0.30,-0.14],
        'Pd-CO':[-0.201,-0.202,-0.217],'Pt-CO':[-0.194,-0.235,-0.191],'Cu-CO':[-0.22,-0.22,-0.20],
        'Au-CO':[-0.137,-0.133,-0.076],'Rh-CO':[-0.135,-0.144,-0.175],'Pd-NO':[-0.21,-0.317,-0.223],
        'Pt-NO':[-0.255,-0.311,-0.211],'Cu-NO':[-0.16,-0.16,-0.09],
        'Au-NO':[-0.058,0,0],'Rh-NO':[-0.155,-0.144,-0.192],
        'Cu-H2O-gas':[0.0,0.0,0.0],"Pd-H2O-gas":[0.0,0.0,0.0],"Pt-H2O-gas":[0.0,0.0,0.0],"Au-H2O-gas":[0.0,0.0,0.0],
        'Cu-H2O-liquid':[0.0,0.0,0.0],"Pd-H2O-liquid":[0.0,0.0,0.0],"Pt-H2O-liquid":[0.0,0.0,0.0],"Au-H2O-liquid":[0.0,0.0,0.0]}
     # thetaML
     K={'Cu-H2O-gas':[0.75,0.84,0.60],"Pd-H2O-gas":[1.0,1.0,0.74],"Pt-H2O-gas":[1.0,0.86,0.77],"Au-H2O-gas":[1.0,1.0,0.82],
        'Cu-H2O-liquid':[0.75,0.84,0.60],"Pd-H2O-liquid":[1.0,1.0,0.74],"Pt-H2O-liquid":[1.0,0.86,0.77],"Au-H2O-liquid":[1.0,1.0,0.82]}
     
     entropy_a={'H2':41.362,'N2':82.394,'O2':90.454,'CO':85.142,'NO':93.121,'CO2':76.458}
     entropy_b={'H2':0.201,'N2':0.148,'O2':0.143,'CO':0.147,'NO':0.143,'CO2':0.181}         

     global NPs
     NPs = var_element 
     global latt_para
     latt_para = L[NPs]  #晶格常数
     global A_atom
     A_atom = np.array(G[NPs])  #
     global gamma
     gamma = np.array(H[NPs])  #
     global Gas
     Gas=[var_environment1 ,var_environment2 ,var_environment3 ,var_environment4 ,var_environment5]
     global GasPP
     GasPP = np.array([var_PP1 ,var_PP2 ,var_PP3 ,var_PP4 ,var_PP5 ])
     global temperature
     temperature = float(var_temperature)
     global pressure
     pressure = float(var_pressure)
     global d
     d = float(var_radius)

     global Flag    
     global theta
     global filename_xyz
     global E_ads
     global bond_length
     global w
     global R
     global k_b
     global unit_coversion

     nGas = 0
     thetaML = np.ones(3)
     rPP = np.array([]) #ratio of partial pressure
     S_entropy = np.array([]) #entropy of gas
     revised_gamma = np.zeros(3)
     element = NPs
     T = temperature
     P = pressure
     LL = np.array([])  #[[NPs-Gas[0]],...]

    #填入 LL[], 不同气体组分的rpp[] & S_entropy[] & thetaML[]
     for i in range(5):
         if Gas[i] == 'None': 
             continue
         LL = np.append(LL,[NPs+'-'+Gas[i]])
         nGas = nGas + 1
         rPP = np.append(rPP,[float(GasPP[i])*P/100.0])
         filename_xyz = filename_xyz+Gas[i]+'-'
         #if Gas[i]=='H2O-gas' or Gas[i]=='H2O-liquid':
         #    thetaML=K[NPs+'-'+Gas[i]]
         if Gas[i] == 'H2O-gas':
             S_entropy = np.append(S_entropy,(188.84-41)*unit_coversion)
             thetaML = K[NPs+'-'+Gas[i]]
         elif Gas[i] == 'H2O-liquid':
             S_entropy = np.append(S_entropy,(69.95-41)*unit_coversion)
         else:
             factor_a = entropy_a[Gas[i]] # entropy_factor_a, CO
             factor_b = entropy_b[Gas[i]]   # entropy_factor_b, CO
             S_entropy = np.append(S_entropy,(factor_a * pow(T, factor_b) - R * np.log(rPP[i] / 100000.0)) * unit_coversion)
     if np.sum(rPP) != P: 
         showinfo("Partial pressure","Something wrong with input Partial Pressures") 
         return 

     sj_start=time.time()
     coverage=np.zeros((3,nGas))     
     E_ads = np.array([M[LL[i]] for i in range(nGas)]) # adosoption energy
     E_ads = E_ads.T

    # 载体参数    
     Sub = var_support     
     SG = {'Cu111-ZnO0001':9.14}
     SM = {'Cu111-ZnO0001-H2O-gas':-0.404}
     SN = {'Cu111-ZnO0001':2}
     SQ = {'Cu111-ZnO0001':0.104}
     SP = {'Cu111-ZnO0001':[18,25]}
     SR = {'Cu111-ZnO0001':-0.023}
     SNP = {'Cu111-ZnO0001':'Cu'}
     if Sub != 'None':
         global E_ads_sub
         global Eadh
         global A_sub
         global gamma_c_s
         nGas = 1
         Gas[0] = 'H2O-gas'
         E_ads_sub = SM[Sub+'-'+Gas[0]] 
         Eadh = SQ[Sub]
         A_sub = SG[Sub]
         gamma_c_s = SR[Sub]
         Plane_sub = SP[Sub]
         index_sub = SN[Sub]
         NPs = SNP[Sub]    

    # if nGas==1: 
    #     rPP[0]=P
    #     S_entropy[0]=factor_a*pow(T,factor_b)-R*np.log(rPP[0]/100000.0)*unit_coversion

     #interaction between gases
     w = np.zeros((3,nGas,nGas))
     for n in range (3):
         for i in range(nGas):
             for j in range(nGas):
                 w[n,i,j] = -np.sqrt(N[LL[i]][n]*N[LL[j]][n])
     filename_xyz = filename_xyz + NPs + '_cluster.xyz'  # filename for .xyz file

     #print LL,rPP,E_ads,w     
     # const
    ####################substrate definition#######################
     CN = np.array([4.0, 2.0, 6.0])  # coordinate num 100 110 111
     planes = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1),
               (1, 1, 0), (-1, -1, 0), (-1, 1, 0), (1, -1, 0), (1, 0, 1), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (0, 1, 1),
               (0, -1, -1), (0, -1, 1), (0, 1, -1), (1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1), (-1, -1, 1), (-1, 1, -1),
               (1, -1, -1), (-1, -1, -1)]
     # const_end

    ####################################################################
    # defalt not calculate the surfcn, change to the value of flag to True to open it
    ## if button_geometry["state"]=="normal":
    ##      global Flag
    ##      Flag=
    ## print(button_geometry["state"]

     def make_grid(*args):
         '''
         make grid of list of arguments
          suppose we have 4 arguments, each has [M, N, P, Q] entries
          the q arguments vary the fastest
          for an arbitrary permute [m, n, p, q], the index is
          q + pQ + nPQ + mNPQ = q+Q(p+P(n+mN))
          '''
         num_var = len(args)
         num_permut = int(reduce(lambda count, arg: count * len(arg), args, 1))
         result = np.zeros((num_permut, num_var))
         for i in range(num_permut):
             index = i
             for j in reversed(range(num_var)):
                 sub_index = int(index % len(args[j]))
                 result[i, j] = args[j][sub_index]
                 index /= len(args[j])
         return result

     def gen_fcc(dim, latt_param):
         '''
         generate fcc structure
         dim of length 3
         '''
         dim = [dim] * 3
         x_dim, y_dim, z_dim = dim
         x_rep, y_rep, z_rep = int(x_dim / latt_param), int(y_dim / latt_param), int(z_dim / latt_param)
         prim_block = np.array(
             [[0., 0., 0.], [latt_param / 2, latt_param / 2, 0.0], [latt_param / 2, 0.0, latt_param / 2],
              [0.0, latt_param / 2, latt_param / 2]])
         reps = make_grid(range(x_rep), range(y_rep), range(z_rep))
         bulk_xyz = np.empty((reps.shape[0] * 4, 3), dtype=float)
         for i, rep in enumerate(reps):
             disp = rep * latt_param
             bulk_xyz[i * 4: i * 4 + 4, :] = prim_block + disp
         center = np.sum(bulk_xyz, axis=0) / bulk_xyz.shape[0]
         bulk_xyz -= center
         return bulk_xyz

     def gen_cluster(bulk_xyz, planes, length, d):
         planes = np.array(planes).T
         planes_norm = np.sqrt(np.sum(planes ** 2, axis=0))
         under_plane_mask = np.sum(np.dot(bulk_xyz, planes) / planes_norm > length, axis=1) == 0
         valid_atoms = np.arange(bulk_xyz.shape[0])[under_plane_mask]
         return bulk_xyz[valid_atoms]
     
     def gen_cluster_sub(bulk_xyz, planes, length, d, X = 'coor_inside' ):
         planes = np.array(planes).T
         planes_norm = np.sqrt(np.sum(planes ** 2, axis=0))
         under_plane_mask = np.sum(np.dot(bulk_xyz, planes) / planes_norm > length, axis=1) == 0
         valid_atoms = np.arange(bulk_xyz.shape[0])[under_plane_mask]
         invalid_atoms = np.setdiff1d(np.arange(bulk_xyz.shape[0]), valid_atoms)
         X = X
         if X == 'coor_inside':
            return bulk_xyz[valid_atoms]
         else:
            return bulk_xyz[invalid_atoms]

     def func(TTT):
         TT = np.zeros((nGas))
         for k in range(nGas):
             TT[k] = TTT[k]
         TT_0 = TTT[nGas]
         cal_w = np.mat(w[i]) * np.mat(TT).T
         cal_w = cal_w.getA().flatten() # 让(1,n)大小的矩阵，变成(n,)大小的数组
        # print cal_w
        # if TTT_0==0.0: TTT_0=10**(-6)
         KC = rPP * np.exp(-(E_ads[i] + T*S_entropy - cal_w*CN[i]) / (k_b*T))
        # print cal_w*CN[i]
        # cal_w=np.exp(cal_w*CN[i]/(k_b*T))             
         ffff = TT / TT_0 - KC
         fff_0 = np.sum(TT) + TT_0 - thetaML[i]
         ffff = np.append(ffff, fff_0)
        # print TT,TT_0,KC,ffff
         return ffff  
     # main body    ####################################################################

     #print(Flag)
     for i in range(3):
         r_ads = np.zeros((nGas))
         r_gamma = 0.0       
         while True:
             theta = fsolve(func,np.random.rand(nGas+1))
            # print E_ads[i],w[i],S_entropy,CN[i],k_b,T,rPP,thetaML[i]
            # print theta,func(theta),thetaML[i]
             if abs(np.sum(func(theta))) < 10**(-6) or theta[0] >= thetaML[i]:
                 break 
         for j in range(nGas):
             if theta[j] >= thetaML[i]:
                 theta = np.zeros((nGas))
                 theta[j] = thetaML[i]
         for k in range(nGas):
             coverage[i][k] = theta[k]
        # print coverage[i]
         cal_w = np.mat(w[i]) * np.mat(coverage[i]).T
         cal_w = cal_w.getA().flatten()
         r_ads = (E_ads[i] - CN[i] * cal_w) / A_atom[i]
         r_gamma = np.mat(coverage[i]) * np.mat(r_ads).T
         revised_gamma[i] = gamma[i] + r_gamma

     surface_energies = [revised_gamma[0]] * 6 + [revised_gamma[1]] * 12 + [revised_gamma[2]] * 8
     length = [i * d / np.min(surface_energies) for i in surface_energies]
     length = np.array(length)
     bulk_dim = np.min(length) * 3
     bulk = gen_fcc(bulk_dim, latt_para)
     
     if sum(revised_gamma > 0) == 3:
         coor_valid = gen_cluster(bulk, planes, length, d)
         if Sub!="None":
             theta_sub = np.arange(0, 1, 0.001)
             f_sub = theta_sub / (1 - theta_sub) - P * np.exp(-(E_ads_sub + T * S_entropy) / (k_b * T)) # HOW ABOUT the entropy change
             coverage_sub = theta_sub[np.argmin(np.abs(f_sub))]
             revised_gamma_c_s = gamma_c_s - coverage_sub*E_ads_sub/A_sub
             if revised_gamma_c_s >= 0:
              # print 'revised_gamma_c_s >= 0'
                height = revised_gamma_c_s*length[Plane_sub[0]]/revised_gamma[index_sub]
                coor_support = gen_cluster_sub(coor_valid, [planes[Plane_sub[0]]], height, d)
             else:
              # print 'revised_gamma_c_s < 0'
                height = - revised_gamma_c_s*length[Plane_sub[0]]/revised_gamma[index_sub]
                #print 'length[-1] = %f' %(length[-1])
                #print 'height = %f' %(height)
                coor_support = gen_cluster_sub(coor_valid, [planes[Plane_sub[1]]], height, d, 'coor_outside')
         #print(Flag)
         if Flag:
             n100, n110, n111, e100110, e100111, e110111, conner, nedge, nsurf, ntotal, surfcn = surf_count(coor_valid, bond_length)
             #shape_factor: SF, SF=Asurf/(V^2/3), Asurf is the surface area, V is the volum of the cluster. Volum per atom = a^3/4
             Asurf = n100*A_atom[0] + n110*A_atom[1] + n111*A_atom[2] + e100110*(A_atom[0]+A_atom[1])/2 + e100111*(A_atom[0]+A_atom[2])/2 + e110111*(A_atom[1]+A_atom[2])/2 + conner*(A_atom[0]+A_atom[1]+A_atom[2])/3
             Volum = ntotal*(latt_para**3)/4
             SF = Asurf/(Volum**(2/3))
             #surface_area_fraction: SAF
             Surf100 = (n100 + e100110/2.0 + e100111/2.0 + conner/3.0 )* A_atom[0]
             Surf110 = (n110 + e100110/2.0 + e110111/2.0 + conner/3.0 )* A_atom[1]
             Surf111 = (n111 + e100111/2.0 + e110111/2.0 + conner/3.0 )* A_atom[2]                     
         else:
             n100, n110, n111, nedge, nsurf, ntotal, surfcn, SF, Surf100, Surf110, Surf111 = [1] * 11         
     else:
         coor_valid = np.array(([0, 0, 0], [0, 0, 0]))
         #print (T, P, 'Negtive surface energy')
         n100, n110, n111, nedge, nsurf, ntotal, surfcn, SF, Surf100, Surf110, Surf111 = [0] * 11

     N_atom = coor_valid.shape[0]

     '''with open(filename_record, 'w') as fp_record:
          djf_record1=['T','P','theta_100','theta_110','theta_111','revised_gamma110','revised_gamma111','tevised_gamma100','n100',
		'n110','n111','nedge','nsurf','ntotal','surfcn','SF','Surf100','Surf110','Surf111']
          djf_record2=[T, P, coverage[0], coverage[1], coverage[2], revised_gamma[0], revised_gamma[1], revised_gamma[2],
                       n100, n110, n111,nedge, nsurf, ntotal, surfcn, SF, Surf100, Surf110,Surf111]
          for i,j in zip(djf_record1,djf_record2):
               fp_record.write(i+'      '+str(j)+'\n')
         fp_record.write('T\nP\ntheta_100\ntheta_110\ntheta_111\nrevised_gamma100\nrevised_gamma110\nrevised_gamma111\nn100\nn110\nn111\nnedge\nnsurf\nntotal\nsurfcn\nSF\nSurf100\nSurf110\nSurf111\n')
         fp_record.write('%d\n%.2f\n%.3f\n%.3f\n%.3f\n%.3f\n%.3f\n%.3f\n%d\n%d\n%d\n%d\n%d\n%d\n%.3f\n%.3f\n%.3f\n%.3f\n%.3f\n' % (
         T, P, coverage[0], coverage[1], coverage[2], revised_gamma[0], revised_gamma[1], revised_gamma[2], n100, n110, n111,
         nedge, nsurf, ntotal, surfcn, SF, Surf100, Surf110,Surf111))'''
     if Sub=='None':
         #print (Sub)
         with open(filename_xyz, 'w') as fp_xyz:
             fp_xyz.write('%d\n' %(N_atom))
             fp_xyz.write('cluster_%d_%d.xyz\n' % (T, P))
             elements = [element] * N_atom
             for i in range(N_atom):
                 fp_xyz.write('%s  %.3f  %.3f  %.3f\n' % (elements[i], coor_valid[i][0], coor_valid[i][1], coor_valid[i][2]))             
     else:
         N_support_atom = coor_support.shape[0]
         with open(filename_xyz, 'w') as fp_xyz:
             fp_xyz.write('%d\n' %(N_support_atom+N_atom))
             fp_xyz.write('cluster_%d_%d.xyz\n' % (T, P))
             elements = [element] * N_support_atom
             for i in range(N_support_atom):
                 fp_xyz.write('%s  %.3f  %.3f  %.3f\n' % (elements[i], coor_support[i][0], coor_support[i][1], coor_support[i][2]))
             for i in range(N_atom):
                 fp_xyz.write('%s  %.3f  %.3f  %.3f\n' % ('H', coor_valid[i][0], coor_valid[i][1], coor_valid[i][2]))
     sj_elapsed = round(time.time() - sj_start , 4)

     #print("Time used:",sj_elapsed)
    
     def show_info():
          if sum(revised_gamma > 0) != 3:
               showerror( '   Unsuitable Condition Entered',"Nanoparticle broken \n\nNegative surface energy ")               
          else:
               answer=askyesno('Geometry Message Info Box','Job Completed. Total Cost About: '+str(sj_elapsed)+' Seconds'+'\nShow the result with Visualization Software?')
               if answer==True:
                    graphics_()
     show_info()

#############################################################
def surf_count(coors, distance_threshold):
    # count the low-index surface (111), (100), (110) atom numbers
    # input coors is numpy array of shape (natoms, 3)
    # return (#atoms of 111, #atoms of 100, #atoms of 110)
    # 111 surface the sum of CN of neighbouring atoms is 90
    # 100 surface the sum of CN of neighbouring atoms is 80
    # 110 surface the sum of CN of neighbouring atoms is 70
    natoms = coors.shape[0]
    cn_mat = np.zeros((natoms, 13), dtype=int) # matrix of coordinate atoms
    coor_number = np.zeros(natoms, dtype=float)
    for i, icoor in enumerate(coors):
        dist_list = np.sqrt(np.sum((icoor - coors) ** 2, axis=1)) # list of distance
        coor_number[i] = np.sum(dist_list < distance_threshold) - 1 # list of coordinate number
        icn = np.arange(natoms)[dist_list < distance_threshold] # list of coordinate atoms
        idx = np.argwhere(icn == i) 
        icn = np.delete(icn, idx)  # delete the element which equal to i
        cn_mat[i][0] = icn.shape[0] 
        cn_mat[i, 1: 1 + icn.shape[0]] = icn
    accum_cn_mat = np.zeros(natoms, dtype=float)
    for i, line in enumerate(cn_mat):
        for j in range(1, line[0] + 1):
            accum_cn_mat[i] += cn_mat[cn_mat[i, j], 0]

    n100 = np.sum(coor_number == 8)
    n111 = np.sum(coor_number == 9)
    nsurf = np.sum(coor_number < 10)
    n110 = 0
    e110111 = 0
    e100111 = 0
    e100110 = 0
    conner = 0
    for i, tn in enumerate(coor_number):
        if tn == 7:
            if accum_cn_mat[i] >= 68:
                n110 += 1
            elif accum_cn_mat[i] >= 65:
                e110111 += 1
            else:
                e100111 += 1
        elif tn == 6:
            if accum_cn_mat[i] >= 56:
                e100110 += 1
            else:
                conner += 1
    nedge = nsurf - n100 - n110 - n111
    ntotal = natoms
    surfcn = float(np.sum(cn_mat[:, 0][cn_mat[:, 0] < 10])) / nsurf
    
    return (n100, n110, n111, e100110, e100111, e110111, conner, nedge, nsurf, ntotal, surfcn)

def graphics_():
    global filename_xyz
    #k=os.path.dirname(os.getcwd())+'\\vmd.exe'
    os.startfile(filename_xyz)
    #win32api.ShellExecute(0,'open',k,filename_xyz,'',1)
#######################################################################

def property_():
     global Flag
     Flag=True
     if (var_element  != 'Pd' ) & (var_element  != 'Pt' ) & (var_element  != 'Cu' ) & (var_element  != 'Au' ) & (var_element  != 'Rh' ):
          showinfo('Element','Please choose the given metal as a selection')
          return
     #if (var_support =='Cu111-ZnO0001' and (var_element !='Cu' or var_environment1 !='H2O-gas')) :
         #showinfo('Support','The only supported nanoparticle in this version is: \nElement : Cu(111)\nSupport : ZnO(0001)\nEnvironment : H2O-gas\nIf you want to calculate the structure of supported metal NP, Please choose the parameters metioned above.')
         #return  
     if (var_radius =='') or (float(var_radius )<=15):
          showinfo('Radius','Please enter the radius greater than 15 angstrom')
          return
     if var_pressure =='':
          showinfo('Pressure','Please enter the pressure')
          return
     if var_temperature =='':
          showinfo('Pressure','Please enter the temperature')
          return
     if (var_environment1 !='' and var_environment1 !='None') and (var_PP1 =='' and
                                                                             (float(var_PP1 )<=0 or float(var_PP1 )>100)):
          showinfo('Partial Pressure 1','Please enter the partial pressure 1 (from 0-100) of environment1')
          return
     if (var_environment2 !='' and var_environment2 !='None') and (var_PP2 =='' and
                                                                             (float(var_PP2 )<=0 or float(var_PP2 )>100)):
          showinfo('Partial Pressure 2','Please enter the partial pressure 2 (from 0-100) of environment2')
          return
     if (var_environment3 !='' and var_environment3 !='None') and (var_PP3 =='' and
                                                                             (float(var_PP3 )<=0 or float(var_PP3 )>100)):
          showinfo('Partial Pressure 3','Please enter the partial pressure 3 (from 0-100) of environment3')
          return
     if (var_environment4 !='' and var_environment4 !='None') and (var_PP4 =='' and
                                                                             (float(var_PP4 )<=0 or float(var_PP4 )>100)):
          showinfo('Partial Pressure 4','Please enter the partial pressure 4 (from 0-100) of environment4')
          return
     if (var_environment5 !='' and var_environment5 !='None') and (var_PP5 =='' and
                                                                             (float(var_PP5 )<=0 or float(var_PP5 )>100)):
          showinfo('Partial Pressure 5','Please enter the partial pressure 5 (from 0-100) of environment5')
          return
     L={'Pd':3.891,'Pt':3.924,'Cu':3.615,'Au':4.078,'Rh':3.804}
     G={'Pd':[7.569,10.704,6.555],'Pt':[7.698,10.887,6.667],'Cu':[6.553,9.239,5.658],
                  'Au':[8.316,11.761,7.202],'Rh':[7.237,10.234,6.627]}
     H={'Pd':[0.094,0.099,0.081],'Pt':[0.111,0.117,0.088],'Cu':[0.091,0.095,0.081],
                  'Au':[0.050,0.057,0.040],'Rh':[0.148,0.146,0.125]}
     M={'Pd-O2':[-0.96,-0.95,-0.39],'Pd-H2':[-0.35,-0.43,-0.25],'Pd-N2':[-0.460,-0.654,-0.233],
        'Pd-CO2':[-0.47,-0.54,-0.02], 'Pd-NO2':[-1.14,-1.33,-0.88],'Pd-CO':[-1.564,-1.611,-1.591],
        'Pd-NO':[-1.75,-1.90,-1.79],
        'Pt-O2':[-0.83,-1.19,-0.21],'Pt-H2':[-0.69,-0.62,-0.48],'Pt-N2':[-0.24,-0.42,-0.01],
        'Pt-CO2':[0,0,-0.02],'Pt-NO2':[-1.35,-1.65,-0.88],'Pt-CO':[-1.84,-1.91,-1.48],
        'Pt-NO':[-1.92,-2.07,-1.47],     
        'Cu-O2':[-1.28,-1.02,-0.16],'Cu-H2':[-0.02,0,-0.01],'Cu-N2':[-0.02,0,-0.01],
        'Cu-CO2':[-0.01,-0.01,-0.03],'Cu-NO2':[-1.22,-1.42,-0.89],'Cu-CO':[-0.58,-0.64,-0.47],
        'Cu-NO':[-0.99,-0.95,-0.70],
        'Au-O2':[-0.22,-0.34,0],'Au-H2':[-0.09,-0.04,0],'Au-N2':[-0.08,-0.03,0],
        'Au-CO2':[-0.10,-0.06,-0.01],'Au-NO2':[-0.78,-0.85,-0.50],'Au-CO':[-0.60,-0.62,-0.22],
        'Au-NO':[-0.62,-0.59,-0.31],
        'Rh-CO':[-1.638,-1.724,-1.624],'Rh-NO':[-2.17,-2.26,-2.11],'Rh-O2':[0,0,0],
        'Cu-H2O-gas':[-0.40,-0.47,-0.29],"Pd-H2O-gas":[-0.39,-0.50,-0.39],"Pt-H2O-gas":[-0.47,-0.56,-0.38],"Au-H2O-gas":[-0.3,-0.3,-0.27],
        'Cu-H2O-liquid':[-0.40,-0.47,-0.29],"Pd-H2O-liquid":[-0.39,-0.50,-0.39],"Pt-H2O-liquid":[-0.47,-0.56,-0.38],"Au-H2O-liquid":[-0.3,-0.3,-0.27]}

     N={'Pd-O2':[-0.33,-0.32,-0.21],'Pd-H2':[-0.06,-0.08,-0.06],'Pd-N2':[-0.189,-0.222,-0.121],
        'Pd-CO2':[-0.31,-0.44,-0.19], 'Pd-NO2':[-0.45,-0.53,-0.31],
        'Pt-O2':[-0.36,-0.40,-0.15],'Pt-H2':[-0.01,0,-0.02],'Pt-N2':[-0.18,-0.18,-0.13],
        'Pt-CO2':[-0.17,-0.16,-0.17],'Pt-NO2':[-0.51,-0.59,-0.30],
        'Cu-O2':[-0.42,-0.54,-0.10],'Cu-H2':[-0.03,-0.03,-0.03],'Cu-N2':[-0.26,-0.23,-0.24],
        'Cu-CO2':[-0.36,-0.35,-0.33],'Cu-NO2':[-0.63,-0.82,-0.47],
        'Au-O2':[0,0,0],'Au-H2':[-0.02,-0.02,0],'Au-N2':[-0.09,-0.08,-0.07],
        'Au-CO2':[-0.10,-0.08,-0.08],'Au-NO2':[-0.24,-0.30,-0.14],
        'Pd-CO':[-0.201,-0.202,-0.217],'Pt-CO':[-0.194,-0.235,-0.191],'Cu-CO':[-0.22,-0.22,-0.20],
        'Au-CO':[-0.137,-0.133,-0.076],'Rh-CO':[-0.135,-0.144,-0.175],'Pd-NO':[-0.21,-0.317,-0.223],
        'Pt-NO':[-0.255,-0.311,-0.211],'Cu-NO':[-0.16,-0.16,-0.09],
        'Au-NO':[-0.058,0,0],'Rh-NO':[-0.155,-0.144,-0.192],
        'Cu-H2O-gas':[0.0,0.0,0.0],"Pd-H2O-gas":[0.0,0.0,0.0],"Pt-H2O-gas":[0.0,0.0,0.0],"Au-H2O-gas":[0.0,0.0,0.0],
        'Cu-H2O-liquid':[0.0,0.0,0.0],"Pd-H2O-liquid":[0.0,0.0,0.0],"Pt-H2O-liquid":[0.0,0.0,0.0],"Au-H2O-liquid":[0.0,0.0,0.0]}

     K={'Cu-H2O-gas':[0.75,0.84,0.60],"Pd-H2O-gas":[1.0,1.0,0.74],"Pt-H2O-gas":[1.0,0.86,0.77],"Au-H2O-gas":[1.0,1.0,0.82],
        'Cu-H2O-liquid':[0.75,0.84,0.60],"Pd-H2O-liquid":[1.0,1.0,0.74],"Pt-H2O-liquid":[1.0,0.86,0.77],"Au-H2O-liquid":[1.0,1.0,0.82]}

     entropy_a={'H2':41.362,'N2':82.394,'O2':90.454,'CO':85.142,'NO':93.121,'CO2':76.458}
     entropy_b={'H2':0.201,'N2':0.148,'O2':0.143,'CO':0.147,'NO':0.143,'CO2':0.181}         

     R = 8.314
     unit_coversion = 0.0000103643
     k_b = 0.000086173303
     coverage = np.zeros(3)
     revised_gamma = np.zeros(3)

     global latt_para
     global A_atom
     global gamma
     global NPs
     global nGas     
     global Gas
     global GasPP
     global bond_length
     global T
     global d
     global P
     element = NPs  # element of metal
     NPs=var_element 
     latt_para=L[NPs]
     A_atom=np.array(G[NPs])
     gamma=np.array(H[NPs])     
     d=float(var_radius )
     T=float(var_temperature )
     P=float(var_pressure )
     nGas=0
     Gas=[var_environment1 ,var_environment2 ,var_environment3 ,var_environment4 ,var_environment5 ]
     GasPP=np.array([var_PP1 ,var_PP2 ,var_PP3 ,var_PP4 ,var_PP5 ])

     global LL
     LL=np.array([])
     global thetaML
     thetaML=np.ones(3)     
     #GasPP=GasPP*P/100.0
     global rPP
     rPP=np.array([])
     global S_entropy
     S_entropy=np.array([])
     global filename_xyz
     filename_xyz=''

     for i in range(5):
         if Gas[i]=='None': continue
         LL=np.append(LL,[NPs+'-'+Gas[i]])
         nGas=nGas+1
         rPP=np.append(rPP,[float(GasPP[i])*P/100.0])
         filename_xyz=filename_xyz+Gas[i]+'-'
         if Gas[i]=='H2O-gas' or Gas[i]=='H2O-liquid':
             thetaML=K[NPs+'-'+Gas[i]]
         if Gas[i]=='H2O-gas':
             S_entropy=np.append(S_entropy,(188.84-41)*unit_coversion)
         elif Gas[i]=='H2O-liquid':
             S_entropy=np.append(S_entropy,(69.95-41)*unit_coversion)
         else:
             factor_a =entropy_a[Gas[i]] # entropy_factor_a, CO
             factor_b = entropy_b[Gas[i]]   # entropy_factor_b, CO
             S_entropy = np.append(S_entropy,(factor_a * pow(T, factor_b) - R * np.log(rPP[i] / 100000.0)) * unit_coversion)

     '''if nGas==1: 
         rPP[0]=P
         S_entropy[0]=factor_a*pow(T,factor_b)-R*np.log(rPP[0]/100000.0)*unit_coversion'''
     if np.sum(rPP)!=P: 
         print ("Something wrong with input Partial Pressures")
         
     global E_ads
     global w
     E_ads=np.array([M[LL[i]] for i in range(nGas)])
     E_ads=E_ads.T
     w=np.zeros((3,nGas,nGas))
     for i in range (3):
         for j in range(nGas):
             for n in range(nGas):
                 w[i,j,n]=-np.sqrt(N[LL[j]][i]*N[LL[n]][i])
    # print LL,rPP,E_ads,w     
     # const

     filename_xyz = filename_xyz+NPs+'_cluster.xyz'  # filename for .xyz file
     CN = np.array([4.0, 2.0, 6.0])  # coordinate num 100 110 111
     planes = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1),
               (1, 1, 0), (-1, -1, 0), (-1, 1, 0), (1, -1, 0), (1, 0, 1), (-1, 0, -1), (1, 0, -1), (-1, 0, 1), (0, 1, 1),
               (0, -1, -1), (0, -1, 1), (0, 1, -1), (1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1), (-1, -1, 1), (-1, 1, -1),
               (1, -1, -1), (-1, -1, -1)]
     # const_end

    ####################################################################

     #global Flag #= False # defalt not calculate the surfcn, change to the value of flag to True to open it
    ## if button_geometry["state"]=="normal":
    ##      global Flag
    ##      Flag=
    ## print(button_geometry["state"]

     def make_grid(*args):
         '''
         make grid of list of arguments
          suppose we have 4 arguments, each has [M, N, P, Q] entries
          the q arguments vary the fastest
          for an arbitrary permute [m, n, p, q], the index is
          q + pQ + nPQ + mNPQ = q+Q(p+P(n+mN))
          '''
         num_var = len(args)
         num_permut = int(reduce(lambda count, arg: count * len(arg), args, 1))
         result = np.zeros((num_permut, num_var))
         for i in range(num_permut):
             index = i
             for j in reversed(range(num_var)):
                 sub_index = int(index % len(args[j]))
                 result[i, j] = args[j][sub_index]
                 index /= len(args[j])
         return result

     def gen_fcc(dim, latt_param):
         '''
         generate fcc structure
         dim of length 3
         '''
         dim = [dim] * 3
         x_dim, y_dim, z_dim = dim
         x_rep, y_rep, z_rep = int(x_dim / latt_param), int(y_dim / latt_param), int(z_dim / latt_param)
         prim_block = np.array(
             [[0., 0., 0.], [latt_param / 2, latt_param / 2, 0.0], [latt_param / 2, 0.0, latt_param / 2],
              [0.0, latt_param / 2, latt_param / 2]])
         reps = make_grid(range(x_rep), range(y_rep), range(z_rep))
         bulk_xyz = np.empty((reps.shape[0] * 4, 3), dtype=float)
         for i, rep in enumerate(reps):
             disp = rep * latt_param
             bulk_xyz[i * 4: i * 4 + 4, :] = prim_block + disp
         center = np.sum(bulk_xyz, axis=0) / bulk_xyz.shape[0]
         bulk_xyz -= center
         return bulk_xyz

     def gen_cluster(bulk_xyz, planes, length, d):
         planes = np.array(planes).T
         planes_norm = np.sqrt(np.sum(planes ** 2, axis=0))
         under_plane_mask = np.sum(np.dot(bulk_xyz, planes) / planes_norm > length, axis=1) == 0
         valid_atoms = np.arange(bulk_xyz.shape[0])[under_plane_mask]
         return bulk_xyz[valid_atoms]

     # main body
     #print(Flag)
     sj_start=time.time()
     coverage=np.zeros((3,nGas))
     for i in range(3):
         r_ads=np.zeros((nGas))
         r_gamma=0.0
         def func(TTT):
            
             TT=np.zeros((nGas))
             for k in range(nGas):
                 TT[k]=TTT[k]
             TT_0=TTT[nGas]
             cal_w=np.mat(w[i])*np.mat(TT).T
            # print cal_w
             cal_w=cal_w.getA().flatten()
            # print cal_w

             
            # if TTT_0==0.0: TTT_0=10**(-6)
             KC=rPP*np.exp(-(E_ads[i]+T*S_entropy-cal_w*CN[i])/(k_b*T))
            # print cal_w*CN[i]
            # cal_w=np.exp(cal_w*CN[i]/(k_b*T))
             
             ffff=TT/TT_0-KC
             fff_0=np.sum(TT)+TT_0-thetaML[i]
             ffff=np.append(ffff,fff_0)
            # print TT,TT_0,KC,ffff
             return ffff
         
         while True:
             theta=fsolve(func,np.random.rand(nGas+1))
        #     print E_ads[i],w[i],S_entropy,CN[i],k_b,T,rPP,thetaML[i]
        #     print theta,func(theta),thetaML[i]
             if abs(np.sum(func(theta)))<10**(-6) or theta[0]>=thetaML[i]: break 
         for j in range(nGas):
             if theta[j]>=thetaML[i]:
                 theta=np.zeros((nGas))
                 theta[j]=thetaML[i]
         for k in range(nGas):
             coverage[i][k] = theta[k]
   #      print coverage[i]
         cal_w=np.mat(w[i])*np.mat(coverage[i]).T
         cal_w=cal_w.getA().flatten()
         r_ads=(E_ads[i]-CN[i]*cal_w)/A_atom[i]
         r_gamma=np.mat(coverage[i])*np.mat(r_ads).T
         revised_gamma[i] = gamma[i] + r_gamma

     surface_energies = [revised_gamma[0]] * 6 + [revised_gamma[1]] * 12 + [revised_gamma[2]] * 8
     length = [i*d/np.min(surface_energies) for i in surface_energies]
     length = np.array(length)
     bulk_dim = np.min(length)*3
     bulk = gen_fcc(bulk_dim, latt_para)
     
     if sum(revised_gamma > 0) == 3:
         #print (T, P, revised_gamma)
         coor_valid = gen_cluster(bulk, planes, length, d)
         #print(Flag)
         if Flag:
             n100, n110, n111, e100110, e100111, e110111, conner, nedge, nsurf, ntotal, surfcn = surf_count(coor_valid, bond_length)
             #shape_factor: SF, SF=Asurf/(V^2/3), Asurf is the surface area, V is the volum of the cluster. Volum per atom = a^3/4
             Asurf = n100*A_atom[0] + n110*A_atom[1] + n111*A_atom[2] + e100110*(A_atom[0]+A_atom[1])/2 + e100111*(A_atom[0]+A_atom[2])/2 + e110111*(A_atom[1]+A_atom[2])/2 + conner*(A_atom[0]+A_atom[1]+A_atom[2])/3
             Volum = ntotal*(latt_para**3)/4
             SF = Asurf/(Volum**(2/3))
             #surface_area_fraction: SAF
             Surf100 = (n100 + e100110/2.0 + e100111/2.0 + conner/3.0 )* A_atom[0]
             Surf110 = (n110 + e100110/2.0 + e110111/2.0 + conner/3.0 )* A_atom[1]
             Surf111 = (n111 + e100111/2.0 + e110111/2.0 + conner/3.0 )* A_atom[2]
                     
         else:
             n100, n110, n111, nedge, nsurf, ntotal, surfcn, SF, Surf100, Surf110, Surf111 = [1] * 11         
     else:
         coor_valid = np.array(([0, 0, 0], [0, 0, 0]))
         #print (T, P, 'Negtive surface energy')
         n100, n110, n111, nedge, nsurf, ntotal, surfcn, SF, Surf100, Surf110, Surf111 = [0] * 11

     N_atom = coor_valid.shape[0]

     with open(filename_record, 'w') as fp_record:
          djf_record1=['T','P','theta_100','theta_110','theta_111','revised_gamma110','revised_gamma111','tevised_gamma100','n100',
		'n110','n111','nedge','nsurf','ntotal','surfcn','SF','Surf100','Surf110','Surf111']
          djf_record2=[T, P, coverage[0], coverage[1], coverage[2], revised_gamma[0], revised_gamma[1], revised_gamma[2],
                       n100, n110, n111,nedge, nsurf, ntotal, surfcn, SF, Surf100, Surf110,Surf111]
          for i,j in zip(djf_record1,djf_record2):
               fp_record.write(i+'           '+str(j)+'\n')


    ## with open(filename_record, 'w') as fp_record:
    ##     fp_record.write('T\nP\ntheta_100\ntheta_110\ntheta_111\nrevised_gamma100\nrevised_gamma110\nrevised_gamma111\nn100\nn110\nn111\nnedge\nnsurf\nntotal\nsurfcn\nSF\nSurf100\nSurf110\nSurf111\n')
    ##     fp_record.write('%d\n%.2f\n%.3f\n%.3f\n%.3f\n%.3f\n%.3f\n%.3f\n%d\n%d\n%d\n%d\n%d\n%d\n%.3f\n%.3f\n%.3f\n%.3f\n%.3f\n' % (
    ## T, P, coverage[0], coverage[1], coverage[2], revised_gamma[0], revised_gamma[1], revised_gamma[2], n100, n110, n111,
    ## nedge, nsurf, ntotal, surfcn, SF, Surf100, Surf110,Surf111))

     with open(filename_xyz, 'w') as fp_xyz:
         fp_xyz.write('%d\n' %(N_atom))
         fp_xyz.write('cluster_%d_%d.xyz\n' % (T, P))
         elements = [element] * N_atom
         for i in range(N_atom):
             fp_xyz.write('%s  %.3f  %.3f  %.3f\n' % (elements[i], coor_valid[i][0], coor_valid[i][1], coor_valid[i][2]))
             
     sj_elapsed = round(time.time() - sj_start , 4)
     #print("Time used:",sj_elapsed)
    
     def show_info():
          if sum(revised_gamma > 0) != 3:
               showerror( '   Unsuitable Condition Entered',"Nanoparticle broken \n\nNegative surface energy ")
          else:
               answer=askyesno('Property Message Info Box','Job Completed. Total Cost About: '+str(sj_elapsed)+' Seconds'+
                               '\nShow the record file?')
               if answer==True:
                    global filename_record
                    #k=os.path.dirname(os.getcwd())+'\\vmd.exe'
                    os.startfile(filename_record)
                    #win32api.ShellExecute(0,'open',k,filename_xyz,'',1)
     show_info()    

######################################################
