# -*- coding: gb2312 -*-
from tkinter import *
from tkinter.messagebox import *
from tkinter.filedialog import *
import tkinter.ttk as ttk
import mosp_MSRM_func as MSRM
import mosp_KMC_func as KMC
import mosp_Input as INPUT
import os
import sys
import json
#import math
#import sys
##import threading

#sys.coinit_flags = 2
os.chdir(sys.path[0]) #设置目录为程序所在文件夹
#产生主窗体window
window=Tk()
window.title('Multiscale Operando Simulation Package')
window.geometry('970x650+230+50')
#window.resizable(0,0)
direc_tory=os.path.abspath('logo.ico')
window.iconbitmap(direc_tory)
#主窗体结束


#主标签设置
tabControl = ttk.Notebook(window) # Create Tab Control
tab0 = ttk.Frame(tabControl) # Add a second tab
tabControl.add(tab0, text='Customization') # Make second tab visible
tab1 = ttk.Frame(tabControl) # Create a tab
tabControl.add(tab1, text='MSRM') # Add the tab
tabControl.pack(expand=1, fill="both") # Pack to make visible
tab2 = ttk.Frame(tabControl) # Add a second tab
tabControl.add(tab2, text='KMC') # Make second tab visible

###################################################################################
#Input 标签框架
Input=LabelFrame(tab0,text='Input',bd=2)
Input.grid(row=0, column=0, sticky=W+N+S+E)
Ba_sic=LabelFrame(Input,bd=0)
Ba_sic.grid(row=0, column=0, sticky=W)
Input_cal =LabelFrame(tab0, text=' Calculation ')
Input_cal.grid(row=0, column=1, padx=10)
Input_para = {"111": "", "101": "", "100": "", "210": "", "211": "", "221": "", "310": "",\
               "311": "", "320": "", "321": "", "322": "", "331": "", "332": ""}
sub_window_flag = 0

#the 1st line
var_input_element=StringVar()
label_input_element=Label(Ba_sic,text='Element:') 
label_input_element.grid(row=0,column=2,padx=15,pady=15,sticky=W)
Input_element=Entry(Ba_sic,textvariable=var_input_element,width=8)
Input_element.grid(row=0,column=3,padx=5,pady=15)
Input_element.focus()

var_input_structure=StringVar()
label_input_structure=Label(Ba_sic,text='Cystal structure:')
label_input_structure.grid(row=0,column=4,padx=15,pady=15,sticky=W)
input_structure_chosen=ttk.Combobox(Ba_sic,width=8,textvariable=var_input_structure)
input_structure_chosen['values']=('FCC','BCC','HCP')
input_structure_chosen.grid(row=0,column=5,padx=5,pady=15)
input_structure_chosen.current(0)

var_input_constant=StringVar()
label_input_constant=Label(Ba_sic,text='Lattice constant(angstrom):')
label_input_constant.grid(row=0,column=6,padx=15,pady=15,sticky=W)
Input_constant=Entry(Ba_sic,textvariable=var_input_constant,width=8)
Input_constant.grid(row=0,column=7,padx=5,pady=15)
Input_para['Lattice-constant'] = var_input_constant.get()

#the 2nd line
var_input_pressure=StringVar()
label_input_pressure=Label(Ba_sic,text='Pressure(Pa):')
label_input_pressure.grid(row=1,column=2,padx=15,pady=15,sticky=W)
Input_pressure=Entry(Ba_sic,textvariable=var_input_pressure,width=8)
Input_pressure.grid(row=1,column=3,padx=5,pady=15)
Input_para['Pressure'] = var_input_pressure.get()

var_input_temperature=StringVar()
label_input_temperature=Label(Ba_sic,text='Temperature(T):')
label_input_temperature.grid(row=1,column=4,padx=15,pady=15,sticky=W)
Input_temperature=Entry(Ba_sic,textvariable=var_input_temperature,width=8)
Input_temperature.grid(row=1,column=5,padx=5,pady=15)
Input_para['Temperature'] = var_input_temperature.get()

var_input_radius=StringVar()
label_input_radius=Label(Ba_sic,text='Radius(angstrom):')
label_input_radius.grid(row=1,column=6,padx=15,pady=15,sticky=W)
Input_radius=Entry(Ba_sic,textvariable=var_input_radius,width=8)
Input_radius.grid(row=1,column=7,padx=5,pady=15)
Input_para['Radius'] = var_input_radius.get()

#the 3rd line
var_input_nfaces=StringVar()
label_input_nfaces=Label(Ba_sic,text='Number of faces:')
label_input_nfaces.grid(row=2,column=2,padx=15,pady=15,sticky=W)
Input_nfaces=Entry(Ba_sic,textvariable=var_input_nfaces,width=8)
Input_nfaces.grid(row=2,column=3,padx=5,pady=15)

#弹出新窗口用于输入界面能信息
def newwindow():
     surface_window = Toplevel(window)
     surface_window.title('Sepecific sufaces and corresponding surface energy')
     surface_window.geometry('400x450')
     
     new_label1 = Label(surface_window, text='surface')
     new_label1.grid(row=0,column=0,padx=10,pady=5)
     new_label2 = Label(surface_window, text='surface energy')
     new_label2.grid(row=0,column=1,padx=20,pady=5) 
     new_label3 = Label(surface_window, text='surface')
     new_label3.grid(row=0,column=2,padx=10,pady=5)
     new_label4 = Label(surface_window, text='surface energy')
     new_label4.grid(row=0,column=3,padx=10,pady=5)

     var_e111 = StringVar()
     new_label_111 = Label(surface_window, text='111')
     new_label_111.grid(row=1,column=0,padx=10,pady=5)
     new_e111 = Entry(surface_window, textvariable=var_e111, width=8)
     new_e111.grid(row=1,column=1,padx=20,pady=5)

     var_e100 = StringVar()
     new_label_100 = Label(surface_window, text='100')
     new_label_100.grid(row=2,column=0,padx=10,pady=5)
     new_e100 = Entry(surface_window, textvariable=var_e100, width=8)
     new_e100.grid(row=2,column=1,padx=20,pady=5)

     var_e101 = StringVar()
     new_label_101 = Label(surface_window, text='101')
     new_label_101.grid(row=3,column=0,padx=10,pady=5)
     new_e101 = Entry(surface_window, textvariable=var_e101, width=8)
     new_e101.grid(row=3,column=1,padx=20,pady=5)   

     var_e210 = StringVar()
     new_label_210 = Label(surface_window, text='210')
     new_label_210.grid(row=4,column=0,padx=10,pady=5)
     new_e210 = Entry(surface_window, textvariable=var_e210, width=8)
     new_e210.grid(row=4,column=1,padx=20,pady=5) 

     var_e211 = StringVar()
     new_label_211 = Label(surface_window, text='211')
     new_label_211.grid(row=5,column=0,padx=10,pady=5)
     new_e211 = Entry(surface_window, textvariable=var_e211, width=8)
     new_e211.grid(row=5,column=1,padx=20,pady=5)      

     var_e221 = StringVar()
     new_label_221 = Label(surface_window, text='221')
     new_label_221.grid(row=6,column=0,padx=10,pady=5)
     new_e221 = Entry(surface_window, textvariable=var_e221, width=8)
     new_e221.grid(row=6,column=1,padx=20,pady=5)

     var_e310 = StringVar()
     new_label_310 = Label(surface_window, text='310')
     new_label_310.grid(row=7,column=0,padx=10,pady=5)
     new_e310 = Entry(surface_window, textvariable=var_e310, width=8)
     new_e310.grid(row=7,column=1,padx=20,pady=5)

     var_e311 = StringVar()
     new_label_311 = Label(surface_window, text='311')
     new_label_311.grid(row=1,column=2,padx=10,pady=5)
     new_e311 = Entry(surface_window, textvariable=var_e311, width=8)
     new_e311.grid(row=1,column=3,padx=10,pady=5)

     var_e320 = StringVar()
     new_label_320 = Label(surface_window, text='320')
     new_label_320.grid(row=2,column=2,padx=10,pady=5)
     new_e320 = Entry(surface_window, textvariable=var_e320, width=8)
     new_e320.grid(row=2,column=3,padx=10,pady=5)   

     var_e321 = StringVar()
     new_label_321 = Label(surface_window, text='321')
     new_label_321.grid(row=3,column=2,padx=10,pady=5)
     new_e321 = Entry(surface_window, textvariable=var_e321, width=8)
     new_e321.grid(row=3,column=3,padx=10,pady=5)  

     var_e322 = StringVar()
     new_label_322 = Label(surface_window, text='322')
     new_label_322.grid(row=4,column=2,padx=10,pady=5)
     new_e322 = Entry(surface_window, textvariable=var_e322, width=8)
     new_e322.grid(row=4,column=3,padx=10,pady=5)      

     var_e331 = StringVar()
     new_label_331 = Label(surface_window, text='331')
     new_label_331.grid(row=5,column=2,padx=10,pady=5)
     new_e331 = Entry(surface_window, textvariable=var_e331, width=8)
     new_e331.grid(row=5,column=3,padx=10,pady=5)     

     var_e332 = StringVar()
     new_label_332 = Label(surface_window, text='332')
     new_label_332.grid(row=6,column=2,padx=10,pady=5)
     new_e332 = Entry(surface_window, textvariable=var_e332, width=8)
     new_e332.grid(row=6,column=3,padx=10,pady=5)    
     #传递参数
     def get_surface_info():
          global sub_window_flag
          Input_para['111'] = var_e111.get()
          Input_para['101'] = var_e101.get()
          Input_para['100'] = var_e100.get()
          Input_para['210'] = var_e210.get()
          Input_para['211'] = var_e211.get()
          Input_para['221'] = var_e221.get()
          Input_para['310'] = var_e310.get()
          Input_para['311'] = var_e311.get()
          Input_para['320'] = var_e320.get()
          Input_para['321'] = var_e321.get()
          Input_para['322'] = var_e322.get()
          Input_para['331'] = var_e331.get()
          Input_para['332'] = var_e332.get()
          sub_window_flag = 1
     #保留子窗口信息
     if(sub_window_flag == 1):
          var_e111.set(Input_para['111'])
          var_e101.set(Input_para['101'])
          var_e100.set(Input_para['100'])
          var_e210.set(Input_para['210'])
          var_e211.set(Input_para['211'])
          var_e221.set(Input_para['221'])
          var_e310.set(Input_para['310'])
          var_e311.set(Input_para['311'])
          var_e320.set(Input_para['320'])
          var_e321.set(Input_para['321'])
          var_e322.set(Input_para['322'])
          var_e331.set(Input_para['331'])
          var_e332.set(Input_para['332'])
     new_button_submit = Button(surface_window, text='submit', bg='light grey', command=get_surface_info)
     new_button_submit.grid(row=8, column=3, sticky=W) 

button_import = Button(Ba_sic, text='more…', bg='light grey', command=newwindow)
button_import.grid(row=2,column=4,padx=5,pady=15,sticky=W)

separ = ttk.Separator(Ba_sic, orient=HORIZONTAL)
separ.grid(row=3, column=0, columnspan=8, sticky='ew')

#the 4th line-gas
var_input_gas1 = StringVar()
label_input_gas1 = Label(Ba_sic,text='Gas1:')
label_input_gas1.grid(row=4,column=2,padx=15,pady=15,sticky=W)
Input_gas1 = Entry(Ba_sic,textvariable=var_input_gas1,width=8)
Input_gas1.grid(row=4,column=3,padx=5,pady=15)

var_input_gas2 = StringVar()
label_input_gas2 = Label(Ba_sic,text='Gas2:')
label_input_gas2.grid(row=4,column=4,padx=15,pady=15,sticky=W)
Input_gas2 = Entry(Ba_sic,textvariable=var_input_gas2,width=8)
Input_gas2.grid(row=4,column=5,padx=5,pady=15)

var_input_gas3 = StringVar()
label_input_gas3 = Label(Ba_sic,text='Gas3:')
label_input_gas3.grid(row=4,column=6,padx=15,pady=15,sticky=W)
Input_gas3 = Entry(Ba_sic,textvariable=var_input_gas3,width=8)
Input_gas3.grid(row=4,column=7,padx=5,pady=15)

#the 5th line-partial pressure
var_input_pp1 = StringVar()
label_input_pp1 = Label(Ba_sic,text='Partial Pressure1(%):')
label_input_pp1.grid(row=5,column=2,padx=15,pady=15,sticky=W)
Input_pp1 = Entry(Ba_sic,textvariable=var_input_pp1,width=8)
Input_pp1.grid(row=5,column=3,padx=5,pady=15)

var_input_pp2 = StringVar()
label_input_pp2 = Label(Ba_sic,text='Partial Pressure2(%):')
label_input_pp2.grid(row=5,column=4,padx=15,pady=15,sticky=W)
Input_pp2 = Entry(Ba_sic,textvariable=var_input_pp2,width=8)
Input_pp2.grid(row=5,column=5,padx=5,pady=15)

var_input_pp3 = StringVar()
label_input_pp3 = Label(Ba_sic,text='Partial Pressure3(%):')
label_input_pp3.grid(row=5,column=6,padx=15,pady=15,sticky=W)
Input_pp3 = Entry(Ba_sic,textvariable=var_input_pp3,width=8)
Input_pp3.grid(row=5,column=7,padx=5,pady=15)

#the 6th line-adsorption energy
var_input_a_energy1 = StringVar()
label_input_a_energy1 = Label(Ba_sic,text='Adsorption Energy1(eV):')
label_input_a_energy1.grid(row=6,column=2,padx=15,pady=15,sticky=W)
Input_a_energy1 = Entry(Ba_sic,textvariable=var_input_a_energy1,width=8)
Input_a_energy1.grid(row=6,column=3,padx=5,pady=15)

var_input_a_energy2 = StringVar()
label_input_a_energy2 = Label(Ba_sic,text='Adsorption Energy2(eV):')
label_input_a_energy2.grid(row=6,column=4,padx=15,pady=15,sticky=W)
Input_a_energy2 = Entry(Ba_sic,textvariable=var_input_a_energy2,width=8)
Input_a_energy2.grid(row=6,column=5,padx=5,pady=15)

var_input_a_energy3 = StringVar()
label_input_a_energy3 = Label(Ba_sic,text='Adsorption Energy3(eV):')
label_input_a_energy3.grid(row=6,column=6,padx=15,pady=15,sticky=W)
Input_a_energy3 = Entry(Ba_sic,textvariable=var_input_a_energy3,width=8)
Input_a_energy3.grid(row=6,column=7,padx=5,pady=15)

#the 7th line-adsorption entropy
var_input_a_entropy1 = StringVar()
label_input_a_entropy1 = Label(Ba_sic,text='Adsorption Entropy1(ev/T):')
label_input_a_entropy1.grid(row=7,column=2,padx=15,pady=15,sticky=W)
Input_a_entropy1 = Entry(Ba_sic,textvariable=var_input_a_entropy1,width=8)
Input_a_entropy1.grid(row=7,column=3,padx=5,pady=15)

var_input_a_entropy2 = StringVar()
label_input_a_entropy2 = Label(Ba_sic,text='Adsorption Entropy2(ev/T):')
label_input_a_entropy2.grid(row=7,column=4,padx=15,pady=15,sticky=W)
Input_a_entropy2 = Entry(Ba_sic,textvariable=var_input_a_entropy2,width=8)
Input_a_entropy2.grid(row=7,column=5,padx=5,pady=15)

var_input_a_entropy3 = StringVar()
label_input_a_entropy3 = Label(Ba_sic,text='Adsorption Entropy3(ev/T):')
label_input_a_entropy3.grid(row=7,column=6,padx=15,pady=15,sticky=W)
Input_a_entropy3 = Entry(Ba_sic,textvariable=var_input_a_entropy3,width=8)
Input_a_entropy3.grid(row=7,column=7,padx=5,pady=15)

#the 8th line-gaseous entropy
var_input_gas_entropy1 = StringVar()
label_input_gas_entropy1 = Label(Ba_sic,text='Gaseous Entropy1(ev/T):')
label_input_gas_entropy1.grid(row=8,column=2,padx=15,pady=15,sticky=W)
Input_gas_entropy1 = Entry(Ba_sic,textvariable=var_input_gas_entropy1,width=8)
Input_gas_entropy1.grid(row=8,column=3,padx=5,pady=15)

var_input_gas_entropy2 = StringVar()
label_input_gas_entropy2 = Label(Ba_sic,text='Gaseous Entropy2(ev/T):')
label_input_gas_entropy2.grid(row=8,column=4,padx=15,pady=15,sticky=W)
Input_gas_entropy2 = Entry(Ba_sic,textvariable=var_input_gas_entropy2,width=8)
Input_gas_entropy2.grid(row=8,column=5,padx=5,pady=15)

var_input_gas_entropy3 = StringVar()
label_input_gas_entropy3 = Label(Ba_sic,text='Gaseous Entropy3(ev/T):')
label_input_gas_entropy3.grid(row=8,column=6,padx=15,pady=15,sticky=W)
Input_gas_entropy3 = Entry(Ba_sic,textvariable=var_input_gas_entropy3,width=8)
Input_gas_entropy3.grid(row=8,column=7,padx=5,pady=15)

#interaction
separ = ttk.Separator(Ba_sic, orient=HORIZONTAL)
separ.grid(row=9, column=0, columnspan=8, sticky='ew')
label_i = Label(Ba_sic,text='Interaction between gases')
label_i.grid(row=10,column = 4, columnspan=2)

var_input_inter11 = StringVar()
label_input_inter11 = Label(Ba_sic,text='gas1-gas1:')
label_input_inter11.grid(row=11,column=2,padx=15,pady=15,sticky=W)
Input_inter11 = Entry(Ba_sic,textvariable=var_input_inter11,width=8)
Input_inter11.grid(row=11,column=3,padx=5,pady=15)

var_input_inter12 = StringVar()
label_input_inter12 = Label(Ba_sic,text='gas1-gas2:')
label_input_inter12.grid(row=11,column=4,padx=15,pady=15,sticky=W)
Input_inter12 = Entry(Ba_sic,textvariable=var_input_inter12,width=8)
Input_inter12.grid(row=11,column=5,padx=5,pady=15)

var_input_inter13 = StringVar()
label_input_inter13 = Label(Ba_sic,text='gas1-gas3:')
label_input_inter13.grid(row=11,column=6,padx=15,pady=15,sticky=W)
Input_inter13 = Entry(Ba_sic,textvariable=var_input_inter13,width=8)
Input_inter13.grid(row=11,column=7,padx=5,pady=15)

var_input_inter22 = StringVar()
label_input_inter22 = Label(Ba_sic,text='gas2-gas2:')
label_input_inter22.grid(row=12,column=2,padx=15,pady=15,sticky=W)
Input_inter22 = Entry(Ba_sic,textvariable=var_input_inter22,width=8)
Input_inter22.grid(row=12,column=3,padx=5,pady=15)

var_input_inter23 = StringVar()
label_input_inter23 = Label(Ba_sic,text='gas2-gas3:')
label_input_inter23.grid(row=12,column=4,padx=15,pady=15,sticky=W)
Input_inter23 = Entry(Ba_sic,textvariable=var_input_inter23,width=8)
Input_inter23.grid(row=12,column=5,padx=5,pady=15)

var_input_inter33 = StringVar()
label_input_inter33 = Label(Ba_sic,text='gas3-gas3:')
label_input_inter33.grid(row=12,column=6,padx=15,pady=15,sticky=W)
Input_inter33 = Entry(Ba_sic,textvariable=var_input_inter33,width=8)
Input_inter33.grid(row=12,column=7,padx=5,pady=15)

def save_input():
     Input_para['Element'] = var_input_element.get()
     Input_para['Crystal-strcture'] = var_input_structure.get()
     Input_para['Number-of-faces'] = var_input_nfaces.get()
     Input_para['Gas1'] = var_input_gas1.get()
     Input_para['Gas2'] = var_input_gas2.get()
     Input_para['Gas3'] = var_input_gas3.get()
     Input_para['PP1'] = var_input_pp1.get()
     Input_para['PP2'] = var_input_pp2.get()
     Input_para['PP3'] = var_input_pp3.get()
     Input_para['Ads-energy1'] = var_input_a_energy1.get()
     Input_para['Ads-energy2'] = var_input_a_energy2.get()
     Input_para['Ads-energy3'] = var_input_a_energy3.get()
     Input_para['Ads-entropy1'] = var_input_a_entropy1.get()
     Input_para['Ads-entropy2'] = var_input_a_entropy2.get()
     Input_para['Ads-entropy3'] = var_input_a_entropy3.get()
     Input_para['gas-entropy1'] = var_input_gas_entropy1.get()
     Input_para['gas-entropy2'] = var_input_gas_entropy2.get()
     Input_para['gas-entropy3'] = var_input_gas_entropy3.get()
     Input_para['inter11'] = var_input_inter11.get()
     Input_para['inter12'] = var_input_inter12.get()
     Input_para['inter13'] = var_input_inter13.get()
     Input_para['inter22'] = var_input_inter22.get()
     Input_para['inter23'] = var_input_inter23.get()
     Input_para['inter33'] = var_input_inter33.get()

button_property=Button(Input_cal,text='Property',width=10,bd=10,command=lambda:[save_input(),INPUT.get_para(**Input_para)])
button_property.grid(row=1,column=0,padx=10,pady=5,sticky=N)
button_geometry=Button(Input_cal,text='Geometry',width=10,bd=10,command=lambda:[save_input(),INPUT.get_para(**Input_para)])
button_geometry.grid(row=0,column=0,padx=10,pady=5,sticky=N)

###################################################################################
#MSRM 标签框架
MS_RM=LabelFrame(tab1,text='MSRM',bd=2)
MS_RM.grid(row=0,column=0, sticky=W+N+S+E)

Me_tal=LabelFrame(MS_RM,bd=0)
Me_tal.grid(row=0, column=0, sticky=W)

Surr_oundings=LabelFrame(MS_RM, padx=20,bd=0)
Surr_oundings.grid(row=0, column=1,  sticky=W)

Cal_culation =LabelFrame(tab1, text=' Calculation ')
Cal_culation.grid(row=1, column=0)

#MSRM-the 1st line
var_support=StringVar()
label_support=Label(Me_tal,text='Support:')
label_support.grid(row=0,column=2,padx=15,pady=15,sticky=W)
suppoprt_chosen=ttk.Combobox(Me_tal,width=12,textvariable=var_support)
suppoprt_chosen['values']=('None','Cu111-ZnO000-1')
suppoprt_chosen.grid(row=0,column=3,padx=25,sticky=W)
suppoprt_chosen.current(0)

var_environment1=StringVar()
label_environment=Label(Surr_oundings,text='Environment1:')
label_environment.grid(row=0,column=4,sticky=W)
environment1_chosen=ttk.Combobox(Surr_oundings,width=8,textvariable=var_environment1)
environment1_chosen['values']=('H2','O2','N2','CO','NO','CO2','H2O-gas','H2O-liquid')
environment1_chosen.grid(row=0,column=5,padx=5,sticky=W)
environment1_chosen.current(0)

Label(Surr_oundings,text='    Partial Pressure 1 (%):').grid(row=0,column=6,padx=20,sticky=W)
var_PP1=StringVar()
PP1=Entry(Surr_oundings,textvariable=var_PP1,width=8)
PP1.grid(row=0,column=7,padx=5,pady=15)
PP1.insert(END,100)

#MSRM-the 2nd line
var_element=StringVar()
label_element=Label(Me_tal,text='Element:')
label_element.grid(row=1,column=2,padx=15,pady=15,sticky=W)
element_chosen=ttk.Combobox(Me_tal,width=12,textvariable=var_element)
element_chosen['values']=('None','Pd','Cu','Au','Pt','Rh','....')
element_chosen.grid(row=1,column=3,padx=25,sticky=W)
element_chosen.current(0)

var_environment2=StringVar()
label_environment=Label(Surr_oundings,text='Environment2:')
label_environment.grid(row=1,column=4,sticky=W)
environment2_chosen=ttk.Combobox(Surr_oundings,width=8,textvariable=var_environment2)
environment2_chosen['values']=('None','H2','O2','N2','CO','NO','CO2','H2O-gas','H2O-liquid')
environment2_chosen.grid(row=1,column=5,padx=5,pady=15)
environment2_chosen.current(0)

Label(Surr_oundings,text='    Partial Pressure 2 (%):').grid(row=1,column=6,padx=20,sticky=W)
var_PP2=StringVar()
PP2=Entry(Surr_oundings,textvariable=var_PP2,width=8)
PP2.grid(row=1,column=7,padx=5,pady=15)

#MSRM-the 3rd line
Label(Me_tal,text='Enter the radius (angstrom):').grid(row=2,column=2,padx=15,sticky=W)
var_radius=StringVar()
radius_=Entry(Me_tal,textvariable=var_radius,width=8)
radius_.grid(row=2,column=3,padx=5,pady=15)

var_environment3=StringVar()
label_environment=Label(Surr_oundings,text='Environment3:')
label_environment.grid(row=2,column=4,sticky=W)
environment3_chosen=ttk.Combobox(Surr_oundings,width=8,textvariable=var_environment3)
environment3_chosen['values']=('None','H2','O2','N2','CO','NO','CO2','H2O-gas','H2O-liquid')
environment3_chosen.grid(row=2,column=5,padx=5,pady=15)
environment3_chosen.current(0)

Label(Surr_oundings,text='    Partial Pressure 3 (%):').grid(row=2,column=6,padx=20,sticky=W)
var_PP3=StringVar()
PP3=Entry(Surr_oundings,textvariable=var_PP3,width=8)
PP3.grid(row=2,column=7,padx=5,pady=15)

#MSRM-the 4th line
Label(Me_tal,text='Enter the pressure (Pa):').grid(row=3,column=2,padx=15,sticky=W)
var_pressure=StringVar()
pressure_=Entry(Me_tal,textvariable=var_pressure,width=8)
pressure_.grid(row=3,column=3,padx=5,pady=15)

var_environment4=StringVar()
label_environment=Label(Surr_oundings,text='Environment4:')
label_environment.grid(row=3,column=4,sticky=W)
environment4_chosen=ttk.Combobox(Surr_oundings,width=8,textvariable=var_environment4)
environment4_chosen['values']=('None','H2','O2','N2','CO','NO','CO2','H2O-gas','H2O-liquid')
environment4_chosen.grid(row=3,column=5,padx=5,pady=15)
environment4_chosen.current(0)

Label(Surr_oundings,text='    Partial Pressure 4 (%): ').grid(row=3,column=6,padx=20,sticky=W)
var_PP4=StringVar()
PP4=Entry(Surr_oundings,textvariable=var_PP4,width=8)
PP4.grid(row=3,column=7,padx=5,pady=15)

#MSRM-the 5th line
Label(Me_tal,text='Enter the temperature (K): ').grid(row=4,column=2,padx=15,pady=16,sticky=W)
var_temperature=StringVar()
temperature_=Entry(Me_tal,textvariable=var_temperature,width=8)
temperature_.grid(row=4,column=3,padx=5,pady=15)

var_environment5=StringVar()
label_environment=Label(Surr_oundings,text='Environment5:')
label_environment.grid(row=4,column=4,sticky=W)
environment5_chosen=ttk.Combobox(Surr_oundings,width=8,textvariable=var_environment5)
environment5_chosen['values']=('None','H2','O2','N2','CO','NO','CO2','H2O-gas','H2O-liquid')
environment5_chosen.grid(row=4,column=5,padx=5,pady=15)
environment5_chosen.current(0)

Label(Surr_oundings,text='    Partial Pressure 5 (%): ').grid(row=4,column=6,padx=20,sticky=W)
var_PP5=StringVar()
PP5=Entry(Surr_oundings,textvariable=var_PP5,width=8)
PP5.grid(row=4,column=7,padx=5,pady=15)

#MSRM-Calculation bottom
def get_NSRM_parameter():
     MSRM.get_support(var_support.get())
     MSRM.get_element(var_element.get())
     MSRM.get_radiu(var_radius.get())
     MSRM.get_pressure(var_pressure.get())
     MSRM.get_temperature(var_temperature.get())
     MSRM.get_environment1(var_environment1.get())
     MSRM.get_environment2(var_environment2.get())
     MSRM.get_environment3(var_environment3.get())
     MSRM.get_environment4(var_environment4.get())
     MSRM.get_environment5(var_environment5.get())
     MSRM.get_PP1(var_PP1.get())
     MSRM.get_PP2(var_PP2.get())
     MSRM.get_PP3(var_PP3.get())
     MSRM.get_PP4(var_PP4.get())
     MSRM.get_PP5(var_PP5.get())

button_property=Button(Cal_culation,text='Property',width=10,bd=10,command=lambda:[get_NSRM_parameter(),MSRM.property_()])
button_property.grid(row=0,column=1,padx=10,pady=5,sticky=N)
button_geometry=Button(Cal_culation,text='Geometry',width=10,bd=10,command=lambda:[get_NSRM_parameter(),MSRM.geometry_1()])
button_geometry.grid(row=0,column=0,padx=10,pady=5,sticky=N)
###################################################################

#KMC标签框架
KMC=LabelFrame(tab2 ,text='KMC',bd=2)
KMC.grid(row=0,column=0,  sticky=W+N+S+E)

KMC_space_group=Label(KMC,bd=0,text='Space Group : ')
KMC_space_group.grid(row=0,column=0,  sticky=W)
var_space_group=StringVar()
space_group_chosen=ttk.Combobox(KMC,width=8,textvariable=var_space_group)
space_group_chosen['values']=('fcc')
space_group_chosen.grid(row=0,column=1,padx=5,sticky=W,pady=10)
space_group_chosen.current(0)

KMC_element1=Label(KMC,bd=0,text='Element 1 : ')
KMC_element1.grid(row=1,column=0,  sticky=W)
var_KMC_element1=StringVar()
KMC_element_chosen1=ttk.Combobox(KMC,width=8,textvariable=var_KMC_element1)
KMC_element_chosen1['values']=('None','Pd','Au','Pt')
KMC_element_chosen1.grid(row=1,column=1,padx=5,sticky=W,pady=10)
KMC_element_chosen1.current(0)

KMC_element2=Label(KMC,bd=0,text='Element 2 : ')
KMC_element2.grid(row=2,column=0,  sticky=W)
var_KMC_element2=StringVar()
KMC_element_chosen2=ttk.Combobox(KMC,width=8,textvariable=var_KMC_element2)
KMC_element_chosen2['values']=('None','Pd','Au','Pt')
KMC_element_chosen2.grid(row=2,column=1,padx=5,sticky=W,pady=10)
KMC_element_chosen2.current(0)

KMC_step=Label(KMC,bd=0,text='Total Steps : ')
KMC_step.grid(row=4,column=0,  sticky=W)
var_KMC_step=StringVar()
E_KMC_step=Entry(KMC,textvariable=var_KMC_step,width=10)
E_KMC_step.grid(row=4,column=1,padx=5,pady=10)

KMC_gap=Label(KMC, pady=10,bd=0,text='Step interval : ')
KMC_gap.grid(row=5,column=0,  sticky=W)
var_KMC_gap=StringVar()
E_KMC_gap=Entry(KMC,textvariable=var_KMC_gap,width=10)
E_KMC_gap.grid(row=5,column=1,padx=5,pady=10)

KMC_temperature=Label(KMC,bd=0,text='Temperature(K) : ')
KMC_temperature.grid(row=6,column=0,  sticky=W)
var_KMC_temperature=StringVar()
E_KMC_temperature=Entry(KMC,textvariable=var_KMC_temperature,width=10)
E_KMC_temperature.grid(row=6,column=1,padx=5,pady=10)

#KMC-Calculation bottom
def get_KMC_parameter():
     KMC.get_sg(var_space_group.get())
     KMC.get_element1(var_KMC_element1.get())
     KMC.get_element2(var_KMC_element2.get())
     KMC.get_step(var_KMC_step.get())
     KMC.get_gap(var_KMC_gap.get())
     KMC.get_temperature(var_KMC_temperature.get())

KMC_Calculation =LabelFrame(tab2, text=' Calculation ')
KMC_Calculation.grid(column=1, row=0,padx=5)
button_calculation=Button(KMC_Calculation,text='Calculation',width=10,bd=10,command=lambda:[get_KMC_parameter,KMC.KMC_calculation_])
button_calculation.grid(row=0,column=1,padx=10,pady=5,sticky=N)

##################################################################
#产生菜单栏
def open_file():
     filepath = askopenfilename()
     with open(filepath, "r") as f:
          para = f.read()
          paradic = json.loads(para)
     var_input_element.set(paradic['Element'])
     var_input_structure.set(paradic['Crystal-strcture'])
     var_input_nfaces.set(paradic['Number-of-faces'])    
     var_input_gas1.set(paradic['Gas1'])   
     var_input_gas2.set(paradic['Gas2'])   
     var_input_gas3.set(paradic['Gas3'])   
     var_input_pp1.set(paradic['PP1'])   
     var_input_pp2.set(paradic['PP2'])   
     var_input_pp3.set(paradic['PP3'])   
     var_input_a_energy1.set(paradic['Ads-energy1'])   
     var_input_a_energy2.set(paradic['Ads-energy2'])    
     var_input_a_energy3.set(paradic['Ads-energy3'])   
     var_input_a_entropy1.set(paradic['Ads-entropy1'])   
     var_input_a_entropy2.set(paradic['Ads-entropy2'])   
     var_input_a_entropy3.set(paradic['Ads-entropy3'])   
     var_input_gas_entropy1.set(paradic['gas-entropy1'])   
     var_input_gas_entropy2.set(paradic['gas-entropy2'])   
     var_input_gas_entropy3.set(paradic['gas-entropy3'])   
     var_input_inter11.set(paradic['inter11'])   
     var_input_inter12.set(paradic['inter12'])   
     var_input_inter13.set(paradic['inter13'])   
     var_input_inter22.set(paradic['inter22'])   
     var_input_inter23.set( paradic['inter23'])   
     var_input_inter33.set(paradic['inter33'])
     global sub_window_flag
     global Input_para
     sub_window_flag = 1
     Input_para = paradic        

def save_file():
     save_window = Toplevel(window)
     save_window.title('Save')
     save_window.geometry('300x150')
     path_ = StringVar()
     name_ = StringVar()
     def select_path():
          dir = askdirectory()
          path_.set(dir)
     Label(save_window, text='Name:').grid(row=0, column=0,pady=5)
     Entry(save_window,textvariable=name_,width=20).grid(row=0,column=1,padx=5,pady=10,sticky='w')          
     Label(save_window, text='Path:').grid(row=1, column=0,pady=5)
     Entry(save_window,textvariable=path_,width=20).grid(row=1,column=1,padx=5,pady=10)
     Button(save_window,text='select', bg='light grey', command=select_path).grid(row=1,column=2,padx=5,pady=10)
     def save():
          save_input()
          filepath = path_.get() + '/' + name_.get() + '.txt'
          json_para = json.dumps(Input_para)
          try:
               with open(filepath,"w") as f:
                    f.write(json_para)
               showinfo("Save","Project has been saved")
          except:
               showinfo("Error","Please check your path")
     Button(save_window, text='submit', bg='light grey', command=save).grid(row=2, column=1)
     save_window.mainloop()

def search_():
     pass

def metal_():
     pass

def surroundings_():
     pass

def about():
    showinfo('About','Author:\nYi Gao,Beien Zhu,Jun Meng,Manyi Duan,Jifeng Du\n\nVersion:\nMOSP test version tv1.0')

menubar=Menu(window)
window.config(menu=menubar)

filemenu=Menu(menubar,tearoff=0)
editmenu=Menu(menubar,tearoff=0)   
functionmenu=Menu(menubar,tearoff=0)
graphicsmenu=Menu(menubar,tearoff=0)
aboutmenu=Menu(menubar,tearoff=0)
helpmenu=Menu(menubar,tearoff=0)
#filemenu.add_command(label='New Project',accelerator='Ctrl + N',command=New_Project)
#filemenu.add_command(label='Open Project',accelerator='Ctrl + O',command=Open_Project)
#filemenu.add_command(label='Save Project',accelerator='Ctrl + S',command=Save_Project)
#filemenu.add_command(label='Save as..',accelerator='Ctrl + shift + S',command=Save_Project)
#filemenu.add_command(label='Recent Files',command=Recent_Files)
#filemenu.add_command(label='Recent Projects',command=Recent_Projects)

menubar.add_cascade(label='File',menu=filemenu)
menubar.add_cascade(label='Edit',menu=editmenu)
menubar.add_cascade(label='Function',menu=functionmenu)
menubar.add_cascade(label='Graphics',menu=graphicsmenu)
menubar.add_cascade(label='Help',menu=helpmenu)

filemenu.add_command(label='New Project',accelerator='Ctrl + N')
filemenu.add_command(label='Open Project',accelerator='Ctrl + O',command=open_file)
filemenu.add_command(label='Save Project',accelerator='Ctrl + S',command=save_file)
#filemenu.add_command(label='Save as..',accelerator='Ctrl + shift + S',command=save_as)
filemenu.add_separator()
filemenu.add_command(label='Recent Files')
filemenu.add_command(label='Recent Projects')

editmenu.add_command(label='Undo',accelerator='Ctrl + Z')
editmenu.add_command(label='Redo',accelerator='Ctrl + Y')
editmenu.add_separator()
editmenu.add_command(label='Cut',accelerator='Ctrl + X')
editmenu.add_command(label='Copy',accelerator='Ctrl + C')
editmenu.add_command(label='Paste',accelerator='Ctrl + V')
editmenu.add_separator()
editmenu.add_command(label='Select All ',accelerator=' Ctrl + A')
editmenu.add_command(label='Search',accelerator='Ctrl + F',command=search_)
#editmenu.add_command(label='Delete')

##buildmenu.add_command(label='P_and_T',command=p_and_t)
propertymenu=Menu(functionmenu,tearoff=0)
propertymenu.add_command(label='Geometry',command=MSRM.geometry_1)
propertymenu.add_separator()
propertymenu.add_command(label='Property',command=MSRM.property_)

custommenu=Menu(functionmenu,tearoff=0)
custommenu.add_command(label='Geometry',command=MSRM.geometry_1)
custommenu.add_separator()
custommenu.add_command(label='Property',command=MSRM.property_)

functionmenu.add_cascade(label='Customization',menu=custommenu)
functionmenu.add_cascade(label='MSRM',menu=propertymenu)

functionmenu.add_separator()
functionmenu.add_command(label='Statistics')

#helpmenu.add_command(label='Author',command=author)
helpmenu.add_command(label='About',command=about)

graphicsmenu.add_command(label='Visualization',command=MSRM.graphics_)

#菜单栏结束

##t2 = threading.Thread(target=window.mainloop())
##t2.start()
##t2.join()
window.mainloop() # Start GUI













