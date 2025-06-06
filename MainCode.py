# -*- coding: utf-8 -*-
"""
@author: GPwgqzhangli
"""

import numpy as np                  
import matplotlib.pyplot as plt
import time

###############################################################################
g=9.81                     ## Acceleration of gravity
Cz=10                      ## Chezy resistance coefficient
Cf=1/Cz**2
nt=1.5                     
fis=1                      
lamda=0.35                 ## alluvial porosity
R=1.65
D=60*(10)**(-3)            ## grain size
beta=0.015*(10)**(-3)      ## wear coefficient
up_rate_0=3*(10)**(-3)     ## uplift rate
up_rate_1=-4*(10)**(-3)    ## subsidence rate 
Lmr=1                      ## macro-roughness layer thickness
IF=0.005                   ## flood intermittency
pl=0.05                    ## characteristically low probability on the hypsometric curve
ph=1-pl                    ## characteristically high probability on the hypsometric curve
peta=(1-pl)/(ph-pl)
T_critical=0.0495          ## critical shear number
Bc_value=20                ## channel width                 
rB=5
Bd_sub_value=rB*Bc_value
Qw_up=1.45                 ## water discharge  
QWW_up=Qw_up*Bc_value
fb=0                       ## no side wall feed
fb_sub=0                   ## no side wall feed
etaa_initial_value=0.5     ## initial alluvial depth
H_down=0                   ## 0,no boundary      
etab_down=0                ## 0,boundary, etab==0 at downstream; -1, etab is free.
HInitial=0                 ## normal flow
WaterSurface_initial=5     ## if HInitial==1
L=8000                     ## total length

'''Slope series'''
slope_up=0.015592          ## the upstream above knickpoint alluvial reach slope
Slope_series=np.array([0.094,0.094/2])
change_value_series=np.array([0,4000])  
sk0=change_value_series[0]  
sk_initial=sk0
us0=change_value_series[1]  
L1=us0
L_fix_sub=L-us0
slope_qac=0.015592
N_feed=1                    ## q_feed=N_feed x q_ac(S=slope_qac)
NN=Slope_series.size        ## the number of the slopes
Ss=0.7
ETA_0=etaa_initial_value+Slope_series[-1]*(L-change_value_series[-1])+Slope_series[0]*(change_value_series[-1])
H_T_AB=5
eta_top_least=ETA_0+H_T_AB  ## the highest position

'''Input Parameters'''
M=100
dx_bar=1/M
x_bar=np.linspace(0,1,M+1)
dtw=0      
Ntoprint=0      
dta=10**(-3)                 ## alluvial
dtb=10**(-3)                 ## bedrock
Nprint=200                   ## inner loop
NNprint=3600       
pprintwant=1200            ## for loops to print
dwant=600              
Want_to_out=600        
show_plot_unit=2             ## time scale in plot, 0=second,1=day,2=year
show_unit=2                  ## time scale in cvs. 0=second,1=day,2=year
'''for longer run'''
run_num=int(0)
last_run_year=96000.0
##############################################################################
'''Normal depth H and Fr series'''
if Slope_series[-1]==0:
    H_normal_series=(Qw_up**2*Cf/np.delete(Slope_series,-1)/g)**(1/3)
else:
    H_normal_series=(Qw_up**2*Cf/Slope_series/g)**(1/3)   
Fr_normal_series=(Qw_up/H_normal_series)/(g*H_normal_series)**(0.5)
###############################################################################
'''Qac series at initial state'''
if Slope_series[-1]==0:
    T_star_up0_series=((Cz)**2*g)**(-1/3)*(Qw_up)**(2/3)*(np.delete(Slope_series,-1))**(2/3)/R/D
else:
    T_star_up0_series=((Cz)**2*g)**(-1/3)*(Qw_up)**(2/3)*(Slope_series)**(2/3)/R/D
T_star_up00_series=fis*T_star_up0_series-T_critical
T_star_up00_series[np.where(T_star_up00_series<0)[0]]=0          
Qac_up_series=4*D*(R*D*g)**(1/2)*(T_star_up00_series)**nt
'''Qac at upstream end'''
Qac_up=Qac_up_series[0] 
###############################################################################
'''sediment feed'''              
T_star_slope_qac=((Cz)**2*g)**(-1/3)*(Qw_up)**(2/3)*(slope_qac)**(2/3)/R/D
Qac_slope_qac=4*D*(R*D*g)**(1/2)*(fis*T_star_slope_qac-T_critical)**nt  
Qs_feed=N_feed*Qac_slope_qac
Qs_feed_str=str(N_feed)+'x'+'$q_{ac}$(S='+str(slope_qac)+')='+str(round(Qs_feed,6))
###############################################################################
TotalTime=dtb*(365.25*24*60*60)*Nprint*NNprint
everyPlotTime=dtb*(365.25*24*60*60)*Nprint*dwant
everySaveTime=dtb*(365.25*24*60*60)*Nprint*Want_to_out
plot_infor=str(round(TotalTime/(365.25*24*60*60),2))+'yr/'+str(round(everyPlotTime/(365.25*24*60*60),4))+'yr'
save_infor=str(round(TotalTime/(365.25*24*60*60),2))+'yr/'+str(everySaveTime/(365.25*24*60*60))+'yr'
print(plot_infor +'\n'+save_infor)
############x_label#########
dtw_d=str(round(dtw,7))+'s('+str(Ntoprint)+')'
dta_d=str(round(dta,7))+'yrs'
dtb_d=str(round(dtb,7))+'yrs'
labxy_xdxdt='M='+str(M)+' dtw='+dtw_d+' dta='+dta_d+' dtb='+dtb_d
#### time part #####
time_s=TotalTime
time_day=time_s/24/60/60
time_year=time_s/365.25/24/60/60
print('time='+str(time_s)+'s='+str(round(time_day,2))+'day='+str(round(time_year,2))+'yr')
time_str='time=('+str(Nprint)+'x'+str(NNprint)+')='+str(round(time_year,2))+'yr'
time_plot_show=[[]]*3
time_plot_show[0]=str(round(time_s,2))+'s'
time_plot_show[1]=str(round(time_day,2))+'day'
time_plot_show[2]=str(round(time_year,2))+'yrs'
time_series_total=dtb*(365.25*24*60*60)*Nprint*(NNprint+1)
tima_series_pieces=dtb*(365.25*24*60*60)*Nprint
tima_want_series_pieces=dtb*(365.25*24*60*60)*Nprint*Want_to_out
time_total=np.append(0,np.arange(0,time_series_total,tima_series_pieces))
timewant_s=np.append(0,np.arange(0,time_series_total,tima_want_series_pieces))
timewant_day=timewant_s/24/60/60
timewant_year=timewant_s/365.25/24/60/60
timewant_show=[[]]*3
timewant_show[0]=timewant_s
timewant_show[1]=timewant_day
timewant_show[2]=timewant_year
###############################################################################
'''Slope in (x,t)'''
Slope_Eta_up=np.zeros([M+1,M+1])  
Slope_Eta_up[0,0]=1
Slope_Eta_up[0,1]=-1
for i in range(1,M,1):
    Slope_Eta_up[i,i-1]=1/2
    Slope_Eta_up[i,i+1]=-1/2
Slope_Eta_up[M,M-1]=1
Slope_Eta_up[M,M]=-1
'''Slope in (x_bar,t_bar)'''
Slope_Eta_bar_up=np.zeros([M+1,M+1])  
Slope_Eta_bar_up[0,0]=0
Slope_Eta_bar_up[0,1]=0
for i in range(1,M+1,1):
    Slope_Eta_bar_up[i,i-1]=1
    Slope_Eta_bar_up[i,i]=-1

au=0.8
Qa_diff_up=np.zeros([M+1,M+2])
Qa_diff_up[M-1,M-1]=-1
Qa_diff_up[M-1,M]=1
for i in range(0,M-1,1):
    Qa_diff_up[i,i]=-au
    Qa_diff_up[i,i+1]=2*au-1
    Qa_diff_up[i,i+2]=1-au
Qa_diff_up[M,M]=-1
Qa_diff_up[M,M+1]=1

eta_diff_up=np.zeros([M+1,M+1])  
eta_diff_up[0,0]=-1
eta_diff_up[0,1]=1
for i in range(1,M,1):
    eta_diff_up[i,i-1]=-1/2
    eta_diff_up[i,i+1]=1/2
eta_diff_up[M,M-1]=-1
eta_diff_up[M,M]=1
###############################################################################
Qw_2_up=np.zeros([M+1,NNprint+1])
H_2_up=np.zeros([M+1,NNprint+1])
U_2_up=np.zeros([M+1,NNprint+1])
Qac_2_up=np.zeros([M+1,NNprint+1])
Qa_2_up=np.zeros([M+1,NNprint+1])
etaa_2_up=np.zeros([M+1,NNprint+1])
etab_2_up=np.zeros([M+1,NNprint+1])
slope_2_up=np.zeros([M+1,NNprint+1])
slopebed_2_up=np.zeros([M+1,NNprint+1])
p_2_up=np.zeros([M+1,NNprint+1])
pa_2_up=np.zeros([M+1,NNprint+1])
Erosion_2_up=np.zeros([M+1,NNprint+1])
Is_2_up=np.zeros([M+1,NNprint+1])
QAA_2_up=np.zeros([M+1,NNprint+1])
QWW_2_up=np.zeros([M+1,NNprint+1])
Bc_2_up=np.zeros([M+1,NNprint+1])
Bd_2_up=np.zeros([M+1,NNprint+1])
Bs_side_2_up=np.zeros([M+1,NNprint+1])
eta_top_2_up=np.zeros([M+1,NNprint+1])
Is_2_up=np.zeros([M+1,NNprint+1])
QFEED_2_up=np.zeros([M+1,NNprint+1])
xx_2_up=np.zeros([M+1,NNprint+1])
Is_dx_up=np.zeros([M+1,NNprint+1])
sk=np.zeros(NNprint+1)   
sk[0]=sk_initial
skdot=np.zeros(NNprint+1)   
switch_points_d=np.zeros(NNprint+1)   
switch_points_d[0]=L1
switch_x_bar_d=np.zeros(NNprint+1)   
switch_x_bar_d[0]=0.5
switch_points_u=np.zeros(NNprint+1)   
switch_points_u[0]=L1
switch_x_bar_u=np.zeros(NNprint+1)   
switch_x_bar_u[0]=0.5
xx_skd_vary_up=np.zeros(M+1)
Qw_skd_vary_up=np.zeros(M+1)
H_skd_vary_up=np.zeros(M+1)
U_skd_vary_up=np.zeros(M+1)
Qac_skd_vary_up=np.zeros(M+1)
Qa_skd_vary_up=np.zeros(M+1)
etaa_skd_vary_up=np.zeros(M+1)
etab_skd_vary_up=np.zeros(M+1)
slope_skd_vary_up=np.zeros(M+1)
slopebed_skd_vary_up=np.zeros(M+1)
slope_bar_skd_vary_up=np.zeros(M+1)
p_skd_vary_up=np.zeros(M+1)
pa_skd_vary_up=np.zeros(M+1)
E_skd_vary_up=np.zeros(M+1)
Bc_skd_vary_up=np.zeros(M+1)
Bd_skd_vary_up=np.zeros(M+1)
etaT_skd_vary_up=np.zeros(M+1)
Bs_side_skd_vary_up=np.zeros(M+1)
Is_skd_vary_up=np.zeros(M+1)
subQAA_skd_vary_up=np.zeros(M+1)
QWW_skd_vary_up=np.zeros(M+1)
QFEED_skd_vary_up=np.zeros(M+1)
etaa_skd_vary_up_P=np.zeros(M+1)
etab_skd_vary_up_P=np.zeros(M+1)
p_skd_vary_up_P=np.zeros(M+1)
pa_skd_vary_up_P=np.zeros(M+1)
###############################################################################
##### Initial Values #######
if (run_num==0):
    xx_skd_initial_up=x_bar*(L-sk_initial)+sk_initial  
    point_sub_index=np.where(xx_skd_initial_up>L1)[0][0]
    point_sub=xx_skd_initial_up[point_sub_index]
    point_uplift_index=np.where(xx_skd_initial_up<=L1)[0][-1]
    point_uplift=xx_skd_initial_up[point_uplift_index]
    arr_uplift=np.arange(0,point_uplift_index+2,1)
    arr_sub=np.arange(point_sub_index,M+1,1)
   
    M_uplift=xx_skd_initial_up[:point_uplift_index+1].size
    M_sub=xx_skd_initial_up[point_sub_index:].size
   
    etab_skd_initial_up=np.zeros(M+1)
    etab_skd_initial_up[point_sub_index:]=(L-xx_skd_initial_up[point_sub_index:])*Slope_series[1]
    etab_skd_initial_up[:point_uplift_index+1]=(L1-xx_skd_initial_up[:point_uplift_index+1])*Slope_series[0]+ etab_skd_initial_up[point_sub_index]
    etab_skd_initial_up[arr_sub]=(L-xx_skd_initial_up[arr_sub])*Slope_series[1]
    etab_skd_initial_up[arr_uplift]=(xx_skd_initial_up[arr_uplift[-1]]-xx_skd_initial_up[arr_uplift])*Slope_series[0]+ etab_skd_initial_up[point_sub_index]
    etaa_skd_initial_up=etaa_initial_value*np.ones(M+1)
   
    Bc_initial_up=Bc_value*np.ones(M+1)
    Bd_initial_up=np.zeros(M+1)
    Bd_initial_up[:point_uplift_index+1]=Bc_value*np.ones(M_uplift)
    Bd_initial_up[point_sub_index:]=Bd_sub_value*np.ones(M_sub)
    up_rate_initial=np.zeros(M+1)
    up_rate_initial[:point_uplift_index+1]=up_rate_0*np.ones(M_uplift)
    up_rate_initial[point_sub_index:]=up_rate_1*np.ones(M_sub)
   
    etaT_skd_initial_up=(etab_skd_initial_up[0]+H_T_AB)*np.ones(M+1)
    H_skd_initial_up=np.zeros(M+1)
    H_skd_initial_up[arr_uplift]=H_normal_series[0]*np.ones(arr_uplift.size)
    H_skd_initial_up[arr_sub]=H_normal_series[1]*np.ones(arr_sub.size)
    Qw_skd_initial_up=Qw_up*np.ones(M+1)
    U_skd_initial_up=Qw_skd_initial_up/H_skd_initial_up
    Qac_skd_initial_up=Qac_up_series[0]*np.ones(M+1)
   
    slope_skd_initial_up=np.dot(Slope_Eta_up,(etaa_skd_initial_up+etab_skd_initial_up))/dx_bar/(L-sk_initial)
    slopebed_skd_initail_up=np.dot(Slope_Eta_up,(etab_skd_initial_up))/dx_bar/(L-sk_initial)
    slope_bar_skd_initial_up=np.dot(Slope_Eta_bar_up,(etaa_skd_initial_up+etab_skd_initial_up))/dx_bar/(L-sk_initial)
    slope_bar_skd_initial_up[0]=(etaa_skd_initial_up[0]+etab_skd_initial_up[0]-(etaa_skd_initial_up[1]+etab_skd_initial_up[1]))/dx_bar/(L-sk_initial)
    p_skd_initial_up=pl+(ph-pl)*etaa_skd_initial_up/Lmr
    pa_skd_initial_up=(p_skd_initial_up-pl)/(1-pl)
    Qa_skd_initial_up=Qs_feed*np.ones(M+1)
    E_skd_initial_up=IF*beta*Qa_skd_initial_up*(1-pa_skd_initial_up)
    Bs_side_skd_initial_up=(etaT_skd_initial_up-etab_skd_initial_up)/Ss
    Is_skd_initial_up=2*fb* Bs_side_skd_initial_up*E_skd_initial_up
    QAA_skd_initial_up=Qs_feed*Bc_initial_up
    QWW_skd_initial_up=Qw_up*Bc_initial_up
############################################################################
'''Plot BedForm'''
xtitle_p='x (m)'
title_p='Rainbow Canyon - Panamint Valley'

fh=plt.figure(figsize=(18,14),dpi=200)
ax=plt.gca()
plt.plot(xx_skd_initial_up,etab_skd_initial_up,'-ko',label='$\eta_{b}$ initial',lw=1.5)
plt.plot(xx_skd_initial_up,etab_skd_initial_up+etaa_skd_initial_up,'--ro',label='$\eta_{b}$+$\eta_{a}$ initial',lw=1.5)
plt.plot(xx_skd_initial_up,etaT_skd_initial_up,':k',label='H+$\eta_{b}$+$\eta_{a}$ initial',lw=1.5)
plt.plot([-1000,xx_skd_initial_up[0]],[(xx_skd_initial_up[0]+1000)*slope_up+etab_skd_initial_up[0],etab_skd_initial_up[0]],'--r',label='$\eta_{b}$ initial',lw=2.5)
plt.xlabel(xtitle_p,fontsize=24)
plt.ylabel('$\eta_{a}$+$\eta_{b}$ (m)',fontsize=24)
plt.title(title_p,fontsize=24)
plt.legend(loc='upper right',fontsize=24)
ax.tick_params(labelsize=24)
plt.xticks(np.append(0,np.append(change_value_series,L)))
plt.grid()
plt.savefig(str(run_num)+'-Bedform-1.png')

'''plot uplift'''
fh=plt.figure(figsize=(18,14),dpi=200)
ax=plt.gca()
ax1=fh.add_subplot(111)
plt.plot(xx_skd_initial_up[arr_uplift],up_rate_initial[arr_uplift]*10**3,'-ok',lw=2.5,label='uplift')
plt.plot(xx_skd_initial_up[arr_sub],up_rate_initial[arr_sub]*10**3,'--ok',lw=2.5,label='subsidence-lift')
ax1.legend(loc='upper right',fontsize=36)
ax1.set_ylabel('up lift rate $v$ (mm/yr)',fontsize=36,fontweight='bold')
plt.title('up lift & bank location',fontsize=36,fontweight='bold')
plt.xlabel('x (m)',fontsize=36,fontweight='bold')
ax2=ax1.twinx()
ax2.plot(xx_skd_initial_up[arr_uplift],Bc_initial_up[arr_uplift]/2,'-or',lw=2.5,label='bed channel width')
ax2.plot(xx_skd_initial_up[arr_uplift],-Bc_initial_up[arr_uplift]/2,'-or',lw=2.5)
ax2.plot(xx_skd_initial_up[arr_sub],Bd_initial_up[arr_sub]/2,'-.',color='orange',lw=2.5,label='alluvium deposit width')
ax2.plot(xx_skd_initial_up[arr_sub],-Bd_initial_up[arr_sub]/2,'-.',color='orange',lw=2.5)
ax2.set_ylabel('channel bank location (m) ',fontsize=36,fontweight='bold')
ax2.legend(loc='lower left',fontsize=36)
plt.xticks(np.append([0,us0],np.append(change_value_series,L)))
ax1.tick_params(axis='both',which='major',labelsize=36,direction='in',pad=10,width=3,length=8)
ax1.grid()
ax2.tick_params(axis='both',which='major',labelsize=36,direction='in',pad=10,width=3,length=8)
ax2.grid()
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)
ax.spines['right'].set_color('k')
ax.spines['top'].set_color('k')
ax.spines['left'].set_color('k')
ax.spines['bottom'].set_color('k')
ax.spines['right'].set_linewidth(3)
ax.spines['top'].set_linewidth(3)
ax.spines['left'].set_linewidth(3)
ax.spines['bottom'].set_linewidth(3)
plt.savefig(str(run_num)+'-Uplift.png')
#################################################################################    
Qw_2_up[:,0]=Qw_skd_initial_up
H_2_up[:,0]=H_skd_initial_up
U_2_up[:,0]=U_skd_initial_up
Qac_2_up[:,0]=Qac_skd_initial_up
Qa_2_up[:,0]=Qa_skd_initial_up
etaa_2_up[:,0]=etaa_skd_initial_up
etab_2_up[:,0]=etab_skd_initial_up
slope_2_up[:,0]=slope_skd_initial_up
p_2_up[:,0]=p_skd_initial_up
pa_2_up[:,0]=pa_skd_initial_up
Bc_2_up[:,0]=Bc_initial_up
Bd_2_up[:,0]=Bc_initial_up
Bs_side_2_up[:,0]=Bs_side_skd_initial_up
eta_top_2_up[:,0]=etaT_skd_initial_up
Erosion_2_up[:,0]=E_skd_initial_up
Is_2_up[:,0]=Is_skd_initial_up
QAA_2_up[:,0]=QAA_skd_initial_up  
QWW_2_up[:,0]=QWW_skd_initial_up
QFEED_2_up[:,0]=Qa_skd_initial_up
xx_2_up[:,0]=xx_skd_initial_up

Qw_skd_vary_up=Qw_skd_initial_up
H_skd_vary_up=H_skd_initial_up
U_skd_vary_up=U_skd_initial_up
etaa_skd_vary_up=etaa_skd_initial_up
etab_skd_vary_up=etab_skd_initial_up
etaT_skd_vary_up=etaT_skd_initial_up
p_skd_vary_up=p_skd_initial_up
pa_skd_vary_up=pa_skd_initial_up
Bc_skd_vary_up=Bc_initial_up
Bd_skd_vary_up=Bd_initial_up
Qa_skd_vary_up=Qa_skd_initial_up
Bs_side_skd_vary_up=Bs_side_skd_initial_up
E_skd_vary_up=E_skd_initial_up
Is_skd_vary_up=Is_skd_initial_up

xx_skd_vary_up=xx_skd_initial_up
sk_vary_up=sk_initial
skdot0=0     
skdot_initial=0
skdot_vary=skdot_initial

slope_skd_vary_up=slope_skd_initial_up
slopebed_skd_vary_up=slope_skd_initial_up
slope_bar_skd_vary_up=slope_bar_skd_initial_up

skdot[0]=skdot_initial  
###############################################################################
start=time.time()
###############################################################################
for kk in range(0,NNprint,1):
    if kk%pprintwant==0:
         print(kk/pprintwant,'Finish!')
    ###############################################################################
    for k in range(0,Nprint,1):
       ###############################################################################
        slope_bar_skd_vary_up=np.dot(Slope_Eta_up,(etaa_skd_vary_up+etab_skd_vary_up))/(dx_bar*(L-sk_vary_up))
        slope_bar_skd_vary_up[-1]=slope_bar_skd_vary_up[-2]
        H_skd_vary_up=(QWW_up**2*Cf/slope_bar_skd_vary_up/g/Bc_skd_vary_up**2)**(1/3)
        U_skd_vary_up=QWW_up/H_skd_vary_up/Bc_skd_vary_up
        Qw_skd_vary_up=U_skd_vary_up*H_skd_vary_up
        Qw_skd_vary_up[0]=Qw_up   
        QWW_skd_vary_up=Qw_skd_vary_up*Bc_skd_vary_up
       ###################################################################################
        T_star_skd_vary_up=H_skd_vary_up*slope_bar_skd_vary_up/R/D
        T_star_skd_vary1_up=fis*T_star_skd_vary_up-T_critical
        T_star_skd_vary1_up[T_star_skd_vary1_up<0]=0
        Qac_skd_vary_up=4*D*(R*D*g)**(1/2)*(T_star_skd_vary1_up)**nt
        Qa_skd_vary_up=Qac_skd_vary_up*pa_skd_vary_up  
        Qa_skd_vary_up[-1]=0     
        ###############################################################################    
        point_sub_index=np.where(x_bar>(L1-sk_vary_up)/(L-sk_vary_up))[0][0]
        point_uplift_index=np.where(x_bar<=(L1-sk_vary_up)/(L-sk_vary_up))[0][-1]
        arr_uplift=np.arange(0,point_uplift_index+2,1)
        arr_sub=np.arange(point_sub_index,M+1,1)
        switch_vary_d=xx_skd_vary_up[arr_uplift[-1]]
        switch_vary_x_bar_d=x_bar[arr_uplift[-1]]
        switch_vary_u=xx_skd_vary_up[arr_uplift[-2]]
        switch_vary_x_bar_u=x_bar[arr_uplift[-2]]
        delta_u=L1-switch_vary_u
        delta_d=switch_vary_d-L1
        if (delta_u<=delta_d):
            arr_uplift=np.arange(0,point_uplift_index+1,1)
            arr_sub=np.arange(point_sub_index-1,M+1,1)
            switch_vary_d=xx_skd_vary_up[arr_sub[1]]
            switch_vary_x_bar_d=x_bar[arr_sub[1]]
            switch_vary_u=xx_skd_vary_up[arr_uplift[-1]]
            switch_vary_x_bar_u=x_bar[arr_uplift[-1]]
        else:
            arr_uplift=np.arange(0,point_uplift_index+2,1)
            arr_sub=np.arange(point_sub_index,M+1,1)
            switch_vary_d=xx_skd_vary_up[arr_sub[0]]
            switch_vary_x_bar_d=x_bar[arr_sub[0]]
            switch_vary_u=xx_skd_vary_up[arr_uplift[-2]]
            switch_vary_x_bar_u=x_bar[arr_uplift[-2]]            
       
        fb_vary=np.zeros(M+1)
        fb_vary[arr_uplift]=fb*np.ones(arr_uplift.size)
        Bd_skd_vary_up[arr_uplift]=Bc_value*np.ones(arr_uplift.size)
        Bd_skd_vary_up[arr_sub]=Bd_sub_value*np.ones(arr_sub.size)
        E_skd_vary_up=IF*beta*Qa_skd_vary_up*(1-pa_skd_vary_up)
        Bs_side_skd_vary_up=(etaT_skd_vary_up-etab_skd_vary_up)/Ss
        Bs_side_skd_vary_up[arr_sub]=np.zeros(arr_sub.size)              
        Part1=-IF*np.dot(Qa_diff_up,np.append(Bc_value,Bc_skd_vary_up)*np.append(Qs_feed,Qa_skd_vary_up))/dx_bar/(L-sk_vary_up)
        Part2=2*fb_vary*Bs_side_skd_vary_up*E_skd_vary_up
        part30=-skdot_vary*(1-x_bar)
        Part3=part30*np.dot(eta_diff_up,etaa_skd_vary_up)/dx_bar/(L-sk_vary_up)
        Part4=1/(1-lamda)/Bd_skd_vary_up/p_skd_vary_up
        F_etaa_up=(Part4*(Part1+Part2)-Part3)*(365.25*24*60*60)*(dta)
        etaa_skd_vary_up=etaa_skd_vary_up+F_etaa_up
        etaa_skd_vary_up[etaa_skd_vary_up<0]=0
        ###################################################################################
        up_rate_initial[arr_uplift]=up_rate_0*np.ones(arr_uplift.size)
        up_rate_initial[arr_sub]=up_rate_1*np.ones(arr_sub.size)
        Part3B=part30*np.dot(eta_diff_up,etab_skd_vary_up)/dx_bar/(L-sk_vary_up)
        F_etab_up=-Part3B*(365.25*24*60*60)*(dtb)+up_rate_initial*(dtb)-E_skd_vary_up*(365.25*24*60*60)*(dtb)*(Bc_skd_vary_up/Bd_skd_vary_up)
        etab_skd_vary_up=etab_skd_vary_up+F_etab_up#        etab_skd_vary_up[-1]=etab_skd_vary_sub[0]  ## at the subsidence point.
        ##################################################################################
        skdot_vary=0  
        sk_vary_up=sk_vary_up+skdot_vary*dtb*(365.25*24*60*60)
        xx_skd_vary_up=x_bar*(L-sk_vary_up)+sk_vary_up
        XX00=xx_skd_vary_up
        YY00=etab_skd_vary_up[0]+H_T_AB
        etaT_skd_vary_up=np.zeros(M+1)
        etaT_skd_vary_up[arr_uplift]=(eta_top_least+up_rate_0*dtb*(365.25*24*60*60))*np.ones(arr_uplift.size)       
        p_skd_vary_up=pl+(ph-pl)*etaa_skd_vary_up/Lmr
        p_skd_vary_up[etaa_skd_vary_up>=(Lmr*peta)]=1
        p_skd_vary_up[etaa_skd_vary_up<0]=pl                          
        pa_skd_vary_up=(p_skd_vary_up-pl)/(1-pl)        
       ###################################################################################
        slope_skd_vary_up=np.dot(Slope_Eta_up,(etaa_skd_vary_up+etab_skd_vary_up))/dx_bar/(L-sk_vary_up)
        slope_skd_vary_up[-1]=slope_skd_vary_up[-2]
        slopebed_skd_vary_up=np.dot(Slope_Eta_up,etab_skd_vary_up)/dx_bar/(L-sk_vary_up)
        slopebed_skd_vary_up[-1]=slopebed_skd_vary_up[-2]        
        Is_skd_vary_up=2*fb*Bs_side_skd_vary_up*E_skd_vary_up
        QAA_skd_vary_up=Qa_skd_vary_up*Bc_skd_vary_up
        QWW_skd_vary_up=Qw_skd_vary_up*Bc_skd_vary_up
        QFEED_vary_up=Qs_feed*Bc_skd_vary_up
        QAA_skd_vary_up=Qa_skd_vary_up*Bc_skd_vary_up
        QWW_skd_vary_up=Qw_skd_vary_up*Bc_skd_vary_up
        QFEED_vary_up=Qs_feed*Bc_skd_vary_up
 ###############################################################################
 ###############################################################################
    H_2_up[:,kk+1]=H_skd_vary_up
    U_2_up[:,kk+1]=U_skd_vary_up
    Qw_2_up[:,kk+1]=Qw_skd_vary_up
    Qac_2_up[:,kk+1]=Qac_skd_vary_up
    Qa_2_up[:,kk+1]=Qa_skd_vary_up
    etaa_2_up[:,kk+1]=etaa_skd_vary_up
    etab_2_up[:,kk+1]=etab_skd_vary_up
    slope_2_up[:,kk+1]=slope_skd_vary_up
    slopebed_2_up[:,kk+1]=slopebed_skd_vary_up
    p_2_up[:,kk+1]=p_skd_vary_up
    pa_2_up[:,kk+1]=pa_skd_vary_up
    Erosion_2_up[:,kk+1]=E_skd_vary_up
    Bs_side_2_up[:,kk+1]=Bs_side_skd_vary_up
    Bc_2_up[:,kk+1]=Bc_skd_vary_up
    Bd_2_up[:,kk+1]=Bd_skd_vary_up
    eta_top_2_up[:,kk+1]=etaT_skd_vary_up
    Is_2_up[:,kk+1]=Is_skd_vary_up
    QAA_2_up[:,kk+1]=QAA_skd_vary_up
    QWW_2_up[:,kk+1]=QWW_skd_vary_up
    QFEED_2_up[:,kk+1]=QFEED_skd_vary_up
    xx_2_up[:,kk+1]=xx_skd_vary_up
    sk[kk+1]=sk_vary_up
    skdot[kk+1]=skdot_vary
    switch_points_d[kk+1]=switch_vary_d
    switch_x_bar_d[kk+1]=switch_vary_x_bar_d
    switch_points_u[kk+1]=switch_vary_u
    switch_x_bar_u[kk+1]=switch_vary_x_bar_u
###############################################################################
elapsed=(time.time()-start)
start=time.time()
EtaaEtab_up=etaa_2_up+etab_2_up
EtaaH_up=etaa_2_up+H_2_up
WaterSurface_up=EtaaEtab_up+H_2_up
Fr_up=np.zeros([M+1,NNprint+1])
H_Half_up=np.zeros([M+1,NNprint+1])
Fr_up=U_2_up/(g*H_2_up)**(0.5)
Erosion_up=np.zeros([M+1,NNprint+1])
Erosion_up=Erosion_2_up
Erosion_m_yr_up=Erosion_up*(365.25*24*60*60)
Is_2_m_yr_up=Is_2_up*(365.25*24*60*60)
elapsed1=(time.time()-start)
print('Time Used =', elapsed1, 's')
###############################################################################        
TStop=NNprint+int(dwant)
dwant=int(dwant)

Want_N=int(NNprint/Want_to_out)
Qw_2=Qw_2_up
U_2=U_2_up
Qac_2=Qac_2_up
Qa_2=Qa_2_up
QAA_2=QAA_2_up
QWW_2=QWW_2_up
QFEED_2=QFEED_2_up
Fr=Fr_up
Erosion=Erosion_up
Erosion_m_yr=Erosion_m_yr_up
H_2=H_2_up
etaa_2=etaa_2_up
etab_2=etab_2_up
slope_2=slope_2_up
slopebed_2=slopebed_2_up
EtaaEtab=EtaaEtab_up
WaterSurface=WaterSurface_up
p_2=p_2_up
pa_2=pa_2_up
eta_top_2=eta_top_2_up
Is_2=Is_2_up
Is_2_m_yr=Is_2_m_yr_up
Bb_2=Bc_2_up
Ba_2=Bd_2_up
Bs_side_2=Bs_side_2_up
xx_2=xx_2_up

list_data=[Qw_2,U_2,Qac_2,Qa_2,Fr,Erosion_m_yr,H_2,etaa_2,etab_2,slope_2,slopebed_2,EtaaEtab,WaterSurface,p_2,pa_2,Is_2_m_yr,xx_2]
name1=(['qw','U','qac','qa','Fr','Erosion_m_yr'])
name2=(['H','etaa','etab','Slope','SlopeBed','EtaaEtab','WaterSurface','p','pa','SideFeed_m2_yr'])
v_N=int(len(name1))
s_N=int(len(name2))  
list_data1=[[]]*v_N
list_data11=[[]]*v_N
list_data2=[[]]*s_N
list_data21=[[]]*s_N
###############################################################################
title_print=[[]]*2
xtitle_print=[[]]*2
title_print[0]='Rainbow Canyon - Panamint Valley'
title_print[1]=title_p
xtitle_print[0]='reach length x (m)'
xtitle_print[1]=xtitle_p
ytitle_print=(['$q_w$ ($m^2$/s)','U (m/s)','$q_{ac}$ ($m^2$/s)','$q_{a}$ ($m^2$/s)','$F_{r}$','Erosion (m/yr)',
'H (m)','Alluvium $\eta_a$ (m)','Bedrock Elevation $\eta_b$ (m)',
'Slope S','Bedslope $S_b$',
'$\eta_b$+$\eta_a$ (m)','$\eta_b$+$\eta_a$+H (m)',
'p','$p_a$',' Side Wall Feed $I_s$ ($m^2$/yr)'])

for i in range(0,v_N,1):
    print(i)   
    fh=plt.figure(figsize=(18,14))
    ax=plt.gca()
    for j in range(0,TStop,dwant):
        plt.plot(list_data[-1][:,j],list_data[i][:,j],lw=2.5)   
    colormap=plt.cm.brg
    colors=[colormap(idx) for idx in np.linspace(1,0,len(ax.lines))]
    for idx,line in enumerate(ax.lines):
        line.set_color(colors[idx])
    plt.plot(list_data[-1][:,0],list_data[i][:,0],'-*g',label='0',lw=3)
    plt.plot(list_data[-1][:,-1],list_data[i][:,-1],'--*k',label=plot_infor,lw=3)   
    plt.xlabel(xtitle_print[0],fontsize=24)
    plt.ylabel(ytitle_print[i],fontsize=24)
    plt.title(title_print[0],fontsize=24)
    ax.tick_params(labelsize=24)
    plt.xticks(np.append(0,np.append(us0,L)))
    plt.legend(loc='upper right',fontsize=24)
    plt.grid()
    plt.savefig(name1[i]+time_plot_show[show_plot_unit]+'.png')
   
for i in range(v_N,v_N+s_N,1):
    print(i)  
    fh=plt.figure(figsize=(18,14))
    ax=plt.gca()
    for j in range(0,TStop,dwant):
        plt.plot(list_data[-1][:,j],list_data[i][:,j],lw=2.5)  
    colormap=plt.cm.brg
    colors=[colormap(idx) for idx in np.linspace(1,0,len(ax.lines))]
    for idx,line in enumerate(ax.lines):
        line.set_color(colors[idx])       
    plt.plot(list_data[-1][:,0],list_data[i][:,0],'-*g',label='0',lw=3)
    plt.plot(list_data[-1][:,-1],list_data[i][:,-1],'--*k',label=plot_infor,lw=3)       
    plt.xlabel(xtitle_print[0],fontsize=24)
    plt.ylabel(ytitle_print[i],fontsize=24)
    plt.title(title_print[0],fontsize=24)
    ax.tick_params(labelsize=24)
    plt.xticks(np.append(0,np.append(us0,L)))
    plt.legend(loc='upper right',fontsize=24)
    plt.grid()
    plt.savefig(name2[i-v_N]+time_plot_show[show_plot_unit]+'.png')
###############################################################################
for i in range(0,v_N,1):
    print(i)
    list_data11[i]=np.zeros([M+2,Want_N+2])
    list_data11[i][0,:]=timewant_show[show_unit]
    list_data11[i][:,0]=np.zeros(M+2)
    list_data11[i][1:(M+2),1:(Want_N+2)]=list_data[i][:,np.arange(0,(NNprint+1),Want_to_out)]
    np.savetxt(name1[i]+time_plot_show[show_plot_unit]+'.csv',list_data11[i],delimiter=',')
   
for i in range(0,s_N,1):
    print(i)
    if i==(s_N-1):
        list_data21[i]=np.zeros([M+2,Want_N+2])
        list_data21[i][0,:]=timewant_show[show_unit]
        list_data21[i][:,0]=np.zeros(M+2)
        list_data21[i][1:(M+2),1:(Want_N+2)]=list_data[i+v_N][:,np.arange(0,(NNprint+1),Want_to_out)]
        np.savetxt(name2[i]+time_plot_show[show_plot_unit]+'.csv',list_data21[i],delimiter=',')
