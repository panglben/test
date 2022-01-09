from tkinter.messagebox import *
from scipy.optimize import fsolve
from functools import reduce
from itertools import permutations
import tkinter
import warnings
import time
import os
import sys
import numpy as np

os.chdir(os.path.dirname(os.path.realpath(sys.argv[0]))) #设置目录为程序所在文件夹
top = tkinter.Tk() 
top.withdraw()

########################################################
NPs = ''  # nanoparticles
structure = ''  # crystal structure
latt_para_a = 0.0 # lattice parameter
latt_para_c = 0.0
A_atom = np.array([]) # area per atom
gamma = np.array([]) # surface energy
nGas = 3 # number of gases
rPP = np.array([]) # partial pressure
S_gas = np.array([[]]) # adsorption entropy (of each gas)
thetaML = np.array([]) # average surface coverage (obtained from MD)
face_num = 0 # number of faces
face_index = np.array([]) # miller index of each face
E_ads = np.array([[]]) # adosorption energy (of each gas on each face)
S_ads = np.array([[]]) # adosorption entropy (of each gas on each face)
w = np.array([[[]]]) # interaction energy between gases (on each face)
P = 0.0 # pressure
T = 0.0 # temperature
d = 0.0 # radius of nanoparticle
#Flag = False
R = 8.314 
unit_coversion = 0.0000103643
k_b = 0.000086173303
bond_length = 3.0
filename_xyz=''
filename_record = ''  # filename for record coverage file
#############################################################

# pause if warning
def handle_warning(message, category, filename, lineno, file=None, line=None):
    q = 'A warning occurred:\n"' + str(message) + '"\nDo you wish to continue?'
    response = askyesno(title = 'Warnings', message = q)

    if response == False:
        raise category(message)


def get_para():
    global NPs, structure, latt_para_a, latt_para_c, A_atom, gamma, nGas, rPP, thetaML, P, T, d, filename_xyz, filename_record
    global face_num, face_index, CN, E_ads, S_ads, S_gas, w
    gas_flag = [1,1,1]
    
    filepath = os.path.join('INP_MSR.txt')

    input = np.loadtxt(filepath, dtype=str, delimiter=': ')

    NPs = input[0][1]
    structure = input[1][1]
    latt_para = input[2][1].split(',')
    latt_para_a = float(latt_para[0])
    try:
        latt_para_c = float(latt_para[1])
    except:
        if structure == 'HCP':
            showerror('Lattice parameter','Two lattice parameters(a and c) are needed for HCP structure')
            return 0
    d = float(input[3][1])
    T = float(input[4][1])
    P = float(input[5][1])
    if d == 0:
        showerror('Radius','Please enter the NP radius')
        return 0
    gas1 = input[6][1].split(',')
    gas2 = input[7][1].split(',')
    gas3 = input[8][1].split(',')
    PP = np.array([float(gas1[0]), float(gas2[0]), float(gas3[0])])
    if sum(PP) != 100:
        showerror('Partial Pressure','Please enter the resonable pratial presure')
        return 0
    S_gas_buff = np.array([float(gas1[1]), float(gas2[1]), float(gas3[1])])
    filename_xyz = filename_xyz + NPs + '_' + structure + '_T_' + str(T) + '_P_' + str(P) + '_cluster.xyz'
    filename_record = filename_record + NPs + '_' + structure + '_T_' + str(T) + '_P_' + str(P) + '_cluster.xyz'
    # 检测气体数
    for i in range(3):
        if float(PP[i]) == 0:
            gas_flag[i] = 0
            nGas -= 1
        else:
            rPP = np.append(rPP, float(PP[i])*P/100.0)
    S_gas = np.array([])
    for i in range(3):
        if gas_flag[i] == 1:
            S_gas = np.append(S_gas, S_gas_buff[i])
    # 记录晶面参数
    face_num = int(input[9][1])
    face_buff = [] #存储晶面名
    face_para = {} #存储晶面参数
    w = np.zeros(face_num)
    for i in range(face_num):
        face = input[10+i][0]
        if face in face_para:
            showerror('index', 'Repeat crystal plane index')
            return 0
        else:
            face_buff.append(face)
            face_para[face] = input[10+i][1]
            
    E_ads = np.zeros((face_num, nGas))
    S_ads = np.zeros((face_num, nGas))
    w = np.zeros((face_num, nGas, nGas))
    thetaML = np.ones((face_num))
    for m in range(face_num):
        face = face_buff[m]
        face_index = np.append(face_index, face) #存入米勒指数
        index = list(face)
        h = float(index[0])
        k = float(index[1])
        l = float(index[2])
        # 根据 Miller index 与 latatice parameter 计算 晶面单位原子面积
        if structure == 'FCC':
            if h/2 != 0 and k/2 != 0 and l/2 != 0:
                A = latt_para_a**2 * np.sqrt(h**2 + k**2 + l**2)/4.0
            else:
                A = latt_para_a**2 * np.sqrt(h**2 + k**2 + l**2)/2.0
        elif structure == 'BCC':
            if (h+k+l)/2 != 0:
                A = latt_para_a**2 * np.sqrt(h**2 + k**2 + l**2)
            else:
                A = latt_para_a**2 * np.sqrt(h**2 + k**2 + l**2)/2.0
        elif structure == 'HCP':
            if (2*h+k)/3 ==0 and l/2 !=0:
                A = np.sqrt(3.0) * latt_para_a**2 * latt_para_c * np.sqrt(0.75*(h**2 + h*k + k**2)/(latt_para_a**2) + (l/latt_para_c)**2)/4.0
            else:
                A = np.sqrt(3.0) * latt_para_a**2 * latt_para_c * np.sqrt(0.75*(h**2 + h*k + k**2)/(latt_para_a**2) + (l/latt_para_c)**2)/2.0
        A_atom = np.append(A_atom, A)
        # face parameter
        parameter = face_para[face].split(',')
        gamma = np.append(gamma, float(parameter[0]))
        E_ads_buff = np.array([parameter[1], parameter[3], parameter[5]])
        E_ads_buff = E_ads_buff.astype(float)
        S_ads_buff = np.array([parameter[2], parameter[4], parameter[6]])
        S_ads_buff = S_ads_buff.astype(float)
        w_buff = np.zeros((3,3))
        w_buff[0] = [parameter[7], parameter[8], parameter[9]]
        w_buff[1] = [parameter[8], parameter[10], parameter[11]]
        w_buff[2] = [parameter[9], parameter[11], parameter[12]]
        w_buff = w_buff.astype(float)
        x = 0
        for i in range(3):
            if gas_flag[i] == 1:
                E_ads[m,x] = E_ads_buff[i]
                S_ads[m,x] = S_ads_buff[i]
                y = 0
                for j in range(3):
                    if gas_flag[j] == 1:
                        w[m,x,y] = w_buff[i,j]
                        y += 1
                x += 1
    return 1
    #print('S_ADS:' , S_ads)
    #print('s_gas:' , S_gas)
    #print('w:' , w)
    #######################################################################################

def get_planes(index):
    planes = []
    index = list(index)
    h = float(index[0])
    k = float(index[1])
    l = float(index[2])
    if structure == 'BCC' or structure == 'FCC':
        for a in range(2):
            h *= -1
            for b in range(2):
                k *= -1
                for c in range(2):
                    l *= -1
                    for item in permutations([h,k,l]):
                         planes.append(item)
    elif structure == 'HCP':
        H = (2 * h - k) / 3
        K = (2 * k - h) / 3
        I = -(H + K)
        L = l
        for a in range(2):
            H *= -1
            for b in range(2):
                K *= -1
                for c in range(2):
                    I *= -1
                    for d in range(2):
                        L *= -1
                        for item in permutations([H,K,I,L]):
                            h = item[0] - item[2]
                            k = item[1] - item[2]
                            l = item[3]
                            planes.append((h,k,l))
    planes = list(set(planes))
    return planes

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

def gen_bcc(dim, latt_param):
    '''
    generate bcc structure
    dim of length 3
    '''
    dim = [dim] * 3
    x_dim, y_dim, z_dim = dim
    x_rep, y_rep, z_rep = int(x_dim / latt_param), int(y_dim / latt_param), int(z_dim / latt_param)
    prim_block = np.array(
        [[0., 0., 0.], [latt_param / 2, latt_param / 2, latt_param / 2]])
    reps = make_grid(range(x_rep), range(y_rep), range(z_rep))
    bulk_xyz = np.empty((reps.shape[0] * 2, 3), dtype=float)
    for i, rep in enumerate(reps):
        disp = rep * latt_param
        bulk_xyz[i * 2: i * 2 + 2, :] = prim_block + disp
    center = np.sum(bulk_xyz, axis=0) / bulk_xyz.shape[0]
    bulk_xyz -= center
    return bulk_xyz

def gen_hcp(dim, latt_param_a, latt_param_c):
    '''
    generate hcp structure
    dim of length 3
    '''
    dim = [dim] * 3
    x_dim, y_dim, z_dim = dim
    x_rep, y_rep, z_rep = int(x_dim / latt_param_a), int(y_dim / latt_param_a), int(z_dim / latt_param_c)
    prim_block = np.array(
        [[latt_param_a/3.0, 2.0*latt_param_a/3.0, latt_param_c/4.0], [2.0*latt_param_a/3.0, latt_param_a/3.0, 3.0*latt_param_c/4.0]])
    reps = make_grid(range(x_rep), range(y_rep), range(z_rep))
    bulk_xyz = np.empty((reps.shape[0] * 2, 3), dtype=float)
    coord_transform = [[1.0, 0.5, 0.0], [0.0, np.sqrt(3)/2.0, 0.0], [0.0, 0.0, 1.0]]
    for i, rep in enumerate(reps):
        dist = [latt_param_a*rep[0], latt_param_a*rep[1], latt_param_c*rep[2]]
        a1a2c = prim_block + dist
        xyz = np.mat(coord_transform) * np.mat(a1a2c).T
        bulk_xyz[i * 2: i * 2 + 2, :] = xyz.T
    center = np.sum(bulk_xyz, axis=0) / bulk_xyz.shape[0]
    bulk_xyz -= center
    return bulk_xyz

def gen_cluster(bulk_xyz, planes, length, d):
    planes = np.array(planes).T
    planes_norm = np.sqrt(np.sum(planes ** 2, axis=0))
    under_plane_mask = np.sum(np.dot(bulk_xyz, planes) / planes_norm > length, axis=1) == 0
    valid_atoms = np.arange(bulk_xyz.shape[0])[under_plane_mask]
    return bulk_xyz[valid_atoms]

def get_surface_energies():
    planes = []
    surface_energies = []
    planes_num = np.zeros(face_num)  
    coverage = np.zeros((face_num, nGas))  # surface coverage (of each gas on each surface)
    revised_gamma = np.zeros((face_num))
    for m in range(face_num):
        # 求解覆盖度方程组
        theta = np.array([]) 
        def func(TTT):
            TT = np.zeros((nGas))
            for k in range(nGas):
                TT[k] = TTT[k]
            TT_0 = TTT[nGas]
            cal_w = np.mat(w[m]) * np.mat(TT).T
            cal_w = cal_w.getA().flatten() # 让(1,n)大小的矩阵，变成(n,)大小的数组
            # if TTT_0==0.0: TTT_0=10**(-6)
            #KC = rPP * np.exp(-(E_ads[m] - T*(S_ads[m] - S_gas) - cal_w*CN[m]) / (k_b*T))
            KC = rPP * np.exp(-(E_ads[m] - T*(S_ads[m] - S_gas) - cal_w) / (k_b*T))
            # cal_w=np.exp(cal_w*CN[i]/(k_b*T))             
            ffff = TT / TT_0 - KC
            fff_0 = np.sum(TT) + TT_0 - thetaML[m]
            ffff = np.append(ffff,fff_0)
            # print TT,TT_0,KC,ffff
            return ffff      
        r_gamma = 0.0   
        revised_gamma[m] = 0.0
        r_ads = np.zeros((nGas))
        while True:
            theta = fsolve(func,np.random.rand(nGas+1))
            if abs(np.sum(func(theta))) < 10**(-6) or theta[0] >= thetaML[m]:
                break 
        # 计算表面能
        for j in range(nGas):
            if theta[j] >= thetaML[m]:
                theta = np.zeros((nGas))
                theta[j] = thetaML[m]
        for k in range(nGas):
            coverage[m][k] = theta[k]
        cal_w = np.mat(w[m]) * np.mat(coverage[m]).T
        cal_w = cal_w.getA().flatten()
        #r_ads = (E_ads[m] - CN[m] * cal_w) / A_atom[m]
        r_ads = (E_ads[m] - cal_w) / A_atom[m]
        r_gamma = np.mat(coverage[m]) * np.mat(r_ads).T
        revised_gamma[m] = float(gamma[m] + r_gamma)
        #print(face_index[m],':',revised_gamma[m])
        # 记录同族晶面
        plane = get_planes(face_index[m])
        planes_num[m] = len(plane)
        planes += plane
        surface_energies += [revised_gamma[m]] * len(plane)
    return (planes, surface_energies, coverage, revised_gamma)

def geometry_():
    flag_para = get_para() # 该值为0代表参数有问题，不执行计算
    sj_start = time.time()
    if flag_para == 0:
        return
    planes, surface_energies, coverage, revised_gamma = get_surface_energies()
    if sum(revised_gamma > 0) != face_num:
        showerror( '   Unsuitable Condition Entered',"Nanoparticle broken \n\nNegative surface energy ")
        return
    length = [e * d / np.min(surface_energies) for e in surface_energies]
    length = np.array(length)
    bulk_dim = np.min(length) * 3
    if structure == 'FCC':
        bulk = gen_fcc(bulk_dim, latt_para_a)
    elif structure == 'BCC':
        bulk = gen_bcc(bulk_dim, latt_para_a)
    elif structure == 'HCP':
        bulk = gen_hcp(bulk_dim, latt_para_a, latt_para_c)
    coor_valid = gen_cluster(bulk, planes, length, d)
    N_atom = coor_valid.shape[0]
    with open(filename_xyz, 'w') as fp_xyz:
        with open('ini.xyz','w') as kmc_ini:
            fp_xyz.write('%d\n' %(N_atom))
            fp_xyz.write('cluster_%d_%d.xyz\n' % (T, P))
            kmc_ini.write('%d\n\n' %(N_atom))
            elements = [NPs] * N_atom
            for i in range(N_atom):
                fp_xyz.write('%s  %.3f  %.3f  %.3f\n' % (elements[i], coor_valid[i][0], coor_valid[i][1], coor_valid[i][2]))
                kmc_ini.write('%s  %.3f  %.3f  %.3f\n' % (elements[i], coor_valid[i][0], coor_valid[i][1], coor_valid[i][2]))
    sj_elapsed = round(time.time() - sj_start , 4)
    showinfo('MSR Message Info Box','Job Completed. Total Cost About: '+str(sj_elapsed)+' Seconds')

if __name__ == '__main__':
    warnings.showwarning = handle_warning
    geometry_()

