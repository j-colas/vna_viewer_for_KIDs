#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 20:01:37 2019
Calibration data KID
INPUT
- sweep IQ 
- header 
- IQ stream
OUTPUT 
- streams calibrés (amp/phase)
@author: jules
"""

#%%

#import sys
import os
import scipy
import numpy as np
import scipy.interpolate
import scipy.optimize as opt
import matplotlib.pyplot  as plt
import matplotlib.patches as mpatches

#from scipy import stats
from multiprocessing import Pool

#%% FONCTIONS

def get_info(fname,fsweep):
    file  = open(fname,'r')
    lines = file.readlines()
    
    run_number = lines[4][20:27:]
    fs     =  int(lines[6][19::])
    ADCres = int(lines[7][17::])
    Irange = int(lines[8][11::])
    Qrange = int(lines[9][11::])
    nTons  = int(lines[10][11::])
    file.close()
    #tone   = int(lines[17][18::]) # à adapter pour plusieurs KID.s !
    
    tone = 0
    fileS  = open(fsweep,'r')
    lineS = fileS.readlines()
    tonelineS  = lineS[21]
    tone = np.float32(tonelineS[20::])
    fileS.close()

    return [run_number,fs,ADCres,Irange,Qrange,nTons,tone]

def savetxt(resdir,fname,s):
    file = open(resdir+fname,'w+')
    file.write(s)
    file.close()

def save_parameters(resdir,params):
    [run_number,fs,ADCres,Irange,Qrange,nTons,tone] = params
    
    s = ''
    s += 'run number\t'+run_number+'\n'
    s += 'fs [Hz] \t'+str(fs)+'\n'
    s += 'ADC [bit]\t'+str(ADCres)+'\n'
    s += 'I range [mV]\t'+str(Irange)+'\n'
    s += 'Q range [mV]\t'+str(Qrange)+'\n'
    s += '# reso. \t'+str(nTons)+'\n'
    s += 'tone(s) [Hz]\t'+str(tone)
    
    savetxt(resdir,'parameters_'+run_number+'.txt',s)

def save_S21_model(resdir,params,i):    
    #[fr,Ql,Qc,a,phi,tau,alpha,c,r,r0,R,P] = params
    #fit_fr,fit_Qc,fit_Ql,fit_phi0,fit_delay,fit_a,fit_alpha,Pz,Rz,offset = fit_params
    
    [fr,Qc,Ql,phi0,tau,a,alpha,Pz,Rz,offset,xc,yc,r] = params
    Qi = 1./(1./Ql-1./Qc)
    
    s = ''
    s += 'fr [Hz]\t'+str(fr)+'\n'
    s += 'Ql\t'+str(Ql)+'\n'
    s += 'Qc\t'+str(Qc)+'\n'
    s += 'Qc\t'+str(Qi)+'\n'
    s += 'a\t'+str(a)+'\n'
    s += 'phi\t'+str(phi0)+'\n'
    s += 'tau\t'+str(tau)+'\n'
    s += 'alpha\t'+str(alpha)+'\n'
    s += 'c\t'+str((xc,yc))+'\n'
    s += 'r\t'+str(r)+'\n'
    #s += 'r0\t'+str(r0)+'\n'
    s += 'R\t'+str(Rz)+'\n'
    s += 'P\t'+str(Pz)
    
    savetxt(resdir,'S21_'+run_number+'KID'+str(i+1)+'.txt',s)

    
#%%

def find_phase(x,y,phi=0):
    return np.unwrap(np.arctan2(x,y)-phi)

def find_mag(x,y,phi=0):
    return np.sqrt(x**2+y**2)

def get_tone(fname):
    ftone = 0
    file  = open(fname,'r')
    lines = file.readlines()
    toneline  = lines[21]
    ftone = np.float32(toneline[20::])
    file.close()
    return ftone

def calc_R(xc,yc): 
    return np.sqrt((x-xc)**2 + (y-yc)**2)

def f_2(c):
    Ri = calc_R(*c)
    return Ri - Ri.mean()

def Df_2b(c):
    xc, yc     = c
    df2b_dc    = np.empty((len(c), x.size))

    Ri = calc_R(xc, yc)
    df2b_dc[0] = (xc - x)/Ri                   # dR/dxc
    df2b_dc[1] = (yc - y)/Ri                   # dR/dyc
    df2b_dc       = df2b_dc - df2b_dc.mean(axis=1)[:, np.newaxis]    
    return df2b_dc

def find_circle(xx,yy):
    global x,y
    x = xx
    y = yy
    center_init =  np.mean(x),np.mean(y)
    
    c,ier = opt.leastsq(f_2,center_init, Dfun=Df_2b, col_deriv=True)
    Ri_2b = calc_R(*c)
    r     = Ri_2b.mean()
    del x,y
    return c,r

def func_phi(f,phi0,Ql,fr):
    return phi0 + 2*np.arctan(2*Ql*(1-(f/fr)))

def fun_model(f,Qc,fr,Ql,phi0):
    return 1.-(Ql*np.exp(1j*phi0)/Qc)/(1+2j*Ql*(f/fr-1.))

def fit_entire_model(f_data,z_data,fr,Qc,Ql,phi0,delay,a=1.,alpha=0.,maxiter=0):
    '''
    fits the whole model: a*exp(i*alpha)*exp(-2*pi*i*f*delay) * [ 1 - {Ql/Qc*exp(i*phi0)} / {1+2*i*Ql*(f-fr)/fr} ]
    '''
    def funcsqr(p,x):
        fr,absQc,Ql,phi0,delay,a,alpha = p
        return np.array([np.absolute( ( a*np.exp(np.complex(0,alpha))*np.exp(np.complex(0,-2.*np.pi*delay*x[i])) * ( 1 - (Ql/absQc*np.exp(np.complex(0,phi0)))/(np.complex(1,2*Ql*(x[i]-fr)/fr)) )  ) )**2 for i in range(len(x))])
    def residuals(p,x,y):
        fr,absQc,Ql,phi0,delay,a,alpha = p
        err = [np.absolute( y[i] - ( a*np.exp(np.complex(0,alpha))*np.exp(np.complex(0,-2.*np.pi*delay*x[i])) * ( 1 - (Ql/absQc*np.exp(np.complex(0,phi0)))/(np.complex(1,2*Ql*(x[i]-fr)/fr)) )  ) ) for i in range(len(x))]
        return err
    p0 = [fr,Qc,Ql,phi0,delay,a,alpha]
    (popt, params_cov, infodict, errmsg, ier) = opt.leastsq(residuals,p0,args=(np.array(f_data),np.array(z_data)),full_output=True,maxfev=maxiter)
    len_ydata = len(np.array(f_data))
    if (len_ydata > len(p0)) and params_cov is not None:  #p_final[1] is cov_x data  #this caculation is from scipy curve_fit routine - no idea if this works correctly...
        s_sq = (funcsqr(popt, np.array(f_data))).sum()/(len_ydata-len(p0))
        params_cov = params_cov * s_sq
    else:
        params_cov = np.inf
    return popt, params_cov, infodict, errmsg, ier

def periodic_boundary(x,bound):
        return np.fmod(x,bound)-np.trunc(x/bound)*bound
    
def fit_S21_model(fIQ):
    f,I,Q,tone = fIQ
    
    (xc,yc),r = find_circle(I,Q)
    phi = find_phase(I-xc,Q-yc)

    z_data = I+1j*Q
    df = f-tone
    #tone_phi = phi[np.argmin(abs(df))]
    
    p0 = (phi[np.argmin(np.abs(df))],1e5,tone)
    popt,pcov = scipy.optimize.curve_fit(func_phi,f,phi,p0,bounds=(-np.inf,np.inf),maxfev=2000)
    
    beta = periodic_boundary(popt[0]+np.pi,np.pi)
    R = np.complex(I[np.argmin(np.abs(f-popt[2]))],Q[np.argmin(np.abs(f-popt[2]))])
    P = np.complex((xc+r*np.sin(beta)),(yc+r*np.cos(beta)))
    
    a     = np.sqrt(P.real**2+P.imag**2)
    alpha = np.arctan(P.imag/P.real)
    
    #phase_fit = func_phi(f,*popt)

    A = a*np.exp(1j*alpha)
    A = 1./A 
    
    z  = z_data*A
    Pz = P*A
    Rz = R*A
    cz = np.complex(xc,yc)*A 
    rz = np.abs(r*A)
    
    phi0 = -np.arcsin(cz.imag/rz)
    fr = popt[2]
    Ql = -popt[1]
    delay = 0
    Qc = np.abs(Ql/(2*rz*np.exp(-1j*phi0)))

    p00=(0,Ql,fr)
    I0 = z.real
    Q0 = z.imag
    q1 = find_phase(I0-cz.real,Q0-cz.imag)
    popt3,pcov3 = scipy.optimize.curve_fit(func_phi,f,-q1,p00,bounds=(-np.inf,np.inf))
    #phase_fit = func_phi(f,*popt3)
    
    fi0 = popt3[0]
    Ql = abs(popt3[1])
    fr = popt3[2]
   
    pp = fit_entire_model(f,z,fr,Qc,Ql,phi0,delay,a=1.,alpha=0.)
    popt2, params_cov, infodict, errmsg, ier = pp
        
    [fr,Qc,Ql,phi0,tau,afit,alphafit] = popt2
    return fr,Qc,Ql,phi0,tau,a,alpha,Pz,Rz,fi0,xc,yc,r


def plot_fit_S21(resdir,I,Q,params):  
    
    fr,Qc,Ql,phi0,tau,a,alpha,Pz,Rz,fi0,xc,yc,r = params
    
    f_fit = np.linspace(f.min(),f.max(),2000)
    S21_fit = a*np.exp(1j*alpha)*fun_model(f_fit,Qc,fr,Ql,phi0)
    
    idmin = np.argmin(np.abs(f_fit-fr))
    X = np.complex(S21_fit[idmin].real,S21_fit[idmin].imag)
    Pz = Pz*a*np.exp(1j*alpha)
    #phi_z = np.unwrap(np.angle(z))
    #phi_f = np.unwrap(np.angle(S21_fit))
    
    plt.figure(1,figsize=(10,8),clear=True,edgecolor='k')
    fig = plt.gcf()
    fig.subplots_adjust(wspace=.2,hspace=0)
    plt.subplot(121)
    plt.title('IQ plane')
    plt.grid('on',which='both',linestyle=':')
    plt.plot(I,Q,'k.',label='data')
    fig = plt.gcf()
    ax  = fig.gca()
    art = mpatches.Circle((xc,yc),radius=r,color='g',linestyle="--",fill=False)
    ax.add_patch(art)
    
    plt.plot(S21_fit.real,S21_fit.imag,'r',label='best fit')
    plt.plot(X.real,X.imag,'bo',label='on resonance ')
    plt.plot(Pz.real,Pz.imag,'ro',label='off resonance')
    
    plt.axis('equal')
    plt.xlabel(r'$\Re (S_{21})$')
    plt.ylabel(r'$\Im (S_{21})$')
    plt.plot(xc,yc,'kx')
    plt.legend(fontsize=16,ncol=1,shadow=True)
    plt.text(0.55,0.45,r'$(x_c,y_c)$',horizontalalignment='center',
     verticalalignment='center', transform = ax.transAxes)

    plt.subplot(222)
    plt.plot((f-fr)/1000,find_mag(I,Q),'k.')
    plt.plot((f_fit-fr)/1000,np.abs(S21_fit),'r')
    
    plt.ylabel(r'$\vert S_{21} \vert$')

    plt.subplot(224)
    ax11 = plt.gca()
    
    ax11.plot((f-fr)/1000,find_phase(I,Q),'k.',label="data")
    ax11.plot((f_fit-fr)/1000,find_phase(S21_fit.real,S21_fit.imag),'r',label=r"$\phi$")
    #mz = np.mean(phi_z)
    #mf = np.mean(phi_f)
    
    ax11.set_xlabel(r'$\Delta f [kHz]$')
    ax11.set_ylabel(r'$\phi$ [rd]',color='r')
    ax11.tick_params(axis='y', labelcolor='r')

    #phi0 = -phase_fit[np.argmin(np.abs(f-fr))]
    # -fi0
    
    ax22 = ax11.twinx()
    ax22.plot((f-fr)/1000,find_phase(I-xc,Q-yc),color='k',marker='.',linestyle='',label="data")
    ax22.plot((f_fit-fr)/1000,find_phase(S21_fit.real-xc,S21_fit.imag-yc),'b',label=r'$\Delta\theta$')
    
    #ax22.set_xlabel(r'$\Delta f [Hz]$')
    ax22.set_ylabel(r'$\Delta \theta [rd]$',color='b')
    ax22.tick_params(axis='y', labelcolor='b')

    plt.tight_layout()
    plt.savefig(resdir+'plot_fit_S21_'+run_number+'.png')
    #plt.close(1)
    
def to_phase(resdir,datadir,fnameI,fnameQ,center,offset):
    [xc,yc] = center
    cpt   = 0
    chunk = 10000 
    Nmax  = 125000000 # = 1Go
    Nmax  = Nmax*2
    # taille disque = Nmax*8 #float64
    # taille dique  = Nmax*4 #float32   
    T     = 0
    idf   = 1
    
    I = np.memmap(datadir+fnameI,dtype=int)
    Q = np.memmap(datadir+fnameQ,dtype=int)
    
    out = open(resdir+"phase_"+fnameI[:-4]+'_%03d.bin'% idf,'wb') 
    while(cpt<I.size):
        cpt_new = cpt+chunk
        if cpt_new > I.size:
            cpt_new = I.size
        np.float32(find_phase(I[cpt:cpt_new]-xc, Q[cpt:cpt_new]-yc, phi=offset)).tofile(out)
        cpt = cpt_new
        T = T+chunk
        if T >= Nmax:
            out.close()
            out = open(resdir+"phase_"+fnameI[:-4]+'_%03d.bin'% idf,'wb')
            out.seek(0)
            T   = 0 
            idf = idf+1
    out.close()
    print('to_phase() done') 
    
def to_amp(resdir,datadir,fnameI,fnameQ,a,alpha):
    cpt   = 0
    chunk = 10000 
    Nmax  = 125000000 # = 1Go
    Nmax  = Nmax*2
    # taille disque = Nmax*8 #float64
    # taille dique  = Nmax*4 #float32   
    T     = 0
    idf   = 1
    
    A = a*np.exp(1j*alpha)
    A = 1./A
    
    I = np.memmap(datadir+fnameI,dtype=int)
    Q = np.memmap(datadir+fnameQ,dtype=int)
    
    out = open(resdir+"amp_"+fnameI[:-4]+'_%03d.bin'% idf,'wb') 
    while(cpt<I.size):
        cpt_new = cpt+chunk
        if cpt_new > I.size:
            cpt_new = I.size
        z = I[cpt:cpt_new]+1j*Q[cpt:cpt_new]
        z *= A
        np.float32(np.sqrt(np.abs(z))).tofile(out)
        
        cpt = cpt_new
        T = T+chunk
        if T >= Nmax:
            out.close()
            out = open(resdir+"amp_"+fnameI[:-4]+'_%03d.bin'% idf,'wb')
            out.seek(0)
            T   = 0 
            idf = idf+1
    out.close()
    print('to_mag() done') 
    
def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0)) 
    return (cumsum[N:] - cumsum[:-N]) / float(N)
    
def phase_to_freq(path,params): 
        files = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path,i)) and 'phase_Sion' in i]
        files   = sorted(files,key=lambda x: int(os.path.splitext(x)[0][-3::])) # sort list
        
        fr,Qc,Ql,phi0,tau,a,alpha,Pz,Rz,fi0,xc,yc,r = fit_params[0]
        f_fit = np.linspace(f.min(),f.max(),2000)
        S21_fit = a*np.exp(1j*alpha)*fun_model(f_fit,Qc,fr,Ql,phi0)
        fit_phase = find_phase(S21_fit.real-xc,S21_fit.imag-yc)
        offset = fit_phase[np.argmin(abs(f_fit-fr))]
        
        df = f_fit-fr
        
        print(offset)
        print(files)
        print('num of file = ',len(files))
        
        data = np.memmap(path+files[0],dtype=np.float32,mode='r')[1000]
        mm = np.median(data)
        
        fit_func = scipy.interpolate.interp1d(df,fit_phase,kind="linear",bounds_error=False,fill_value="extrapolate")    
        fit_func2 = scipy.interpolate.interp1d(fit_phase-offset,df,kind="linear",bounds_error=False,fill_value="extrapolate")

        plt.figure(2,figsize=(9,5),clear=True,edgecolor='k')
        plt.subplot(121)
        plt.plot(df,fit_phase,'k.')
        plt.plot(df,fit_func(df),'r')
        plt.plot(df,fit_phase-offset,'b')
        plt.axhline(offset,color='r',ls='--',lw=2)
        plt.axhline(mm,color='b',ls='--',lw=2)
        plt.subplot(122)
        plt.plot(fit_phase-offset,df,'k.')
        plt.plot(fit_phase-offset,fit_func2(fit_phase-offset),'r')
        plt.plot(mm,fit_func2(mm),'bo')
        plt.tight_layout()
        plt.savefig(path+'phase_to_freq_debug.png')
        plt.close(2)
      
        idf   = 1
        for file in files: 
            print("current file = ",file)
            data = np.memmap(path+file,dtype=np.float32,mode='r')
            N     = data.size
            chunk = 10000
            cpt   = 0
            res = np.array([],np.float32)
            while(cpt<N):
                cpt_new = cpt+chunk
                if cpt_new > N:
                    cpt_new = N
                res = np.append(res,np.float32(fit_func2(data[cpt:cpt_new:])))
                cpt = cpt_new
            filename = 'freq_phase_'+file[10:17]+'.%03d.bin'% idf
            out = open(resdir+filename,'wb')
            out.seek(0)
            idf = idf+1
            res.tofile(out)
            out.close()
        print('freq from phase conversion done')

def amp_to_freq(path): 
        files = [i for i in os.listdir(path) if os.path.isfile(os.path.join(path,i)) and 'amp_Sion' in i]
        files   = sorted(files,key=lambda x: int(os.path.splitext(x)[0][-3::])) # sort list
        
        print(files)
        print('num of file = ',len(files))

#%%
#if __name__ == "__main__" : 
#    
#    plt.close('all')
#    
#    num_kid = 1
#    num_cpu = 2
#    
#    datadir = "/home/jules/Documents/GitRep/NEPAL_KID/data_debug/"
#    resdir  = "/home/jules/Documents/GitRep/NEPAL_KID/results/"
#    
#    
##    file_header = "Sion141_115_X_stream_2min_LB44_200mK_1MHz.txt"
##    file_sweep  = "Sion141_114_S.txt"
##    file_I = "Sion141_115_X_stream_2min_LB44_200mK_1MHzI.bin"
##    file_Q = "Sion141_115_X_stream_2min_LB44_200mK_1MHzQ.bin"
#    file_header = "Sion141_040_X_stream_2min_LB32_200mK_1MHz.txt"
#    file_sweep  = "Sion141_039_S.txt"
#    file_I = "Sion141_040_X_stream_2min_LB32_200mK_1MHzI.bin"
#    file_Q = "Sion141_040_X_stream_2min_LB32_200mK_1MHzQ.bin"
# 
#    params = get_info(datadir+file_header,datadir+file_sweep)
#    save_parameters(resdir,params)
#    
#    global run_number, tone
#    run_number = params[0]
#    tone = params[6]
#    
#    data_sweep = np.loadtxt(datadir+file_sweep,comments='%')
#    
#    f = data_sweep[::,0]
#    I = data_sweep[::,1]
#    Q = data_sweep[::,2]
#    fIQ = [(f,I,Q,tone)]
#    
#    for i in range(1,num_kid):
#        I = np.vstack((I,data_sweep[::,2*i+1]))
#        Q = np.vstack((Q,data_sweep[::,2*i+2]))
#        fIQ.append((f,data_sweep[::,2*i+1],data_sweep[::,2*i+2],tone))
#          
#    #Model fit
#    p = Pool(num_cpu)
#    fit_params = p.map(fit_S21_model,fIQ)
#    p.close()
#    p.join()
#
#    for i in range(len(fit_params)):
#        save_S21_model(resdir,fit_params[i],i)
#        plot_fit_S21(resdir,I,Q,fit_params[i])
#
##%% 
#    """
#    ACCELERATION POSSIBLE EN PROFITANT DU POOL() POUR LES FONCTIONS DE CONVERSION
#    """
#    #Gen Phase stream
#    to_phase(resdir,datadir,file_I,file_Q,(fit_params[0][10],fit_params[0][11]),0)
#    #Gen Amp stream
#    to_amp(resdir,datadir,file_I,file_Q,fit_params[0][6],fit_params[0][7])
#    
##%%
#    #Phase to detuning
#    phase_to_freq(resdir,fit_params)
#
#    #Amp to detuning
#    #amp_to_freq(resdir)
#
##%% 
#    signal_phase = np.memmap(resdir+"phase_Sion141_115_X_stream_2min_LB44_200mK_1MHzI_001.bin",dtype=np.float32)[:10000:]
#    signal_amp = np.memmap(resdir+"amp_Sion141_115_X_stream_2min_LB44_200mK_1MHzI_001.bin",dtype=np.float32)[:10000:]
#    signal_freqFromPhase = np.memmap(resdir+'freq_phase_141_115.001.bin',dtype=np.float32)[:10000:]
#    
#    plt.figure()
#    plt.subplot(311)
#    plt.plot(signal_phase,'k')
#    plt.subplot(312)
#    plt.plot(signal_amp,'r')
#    plt.subplot(313)
#    plt.plot(signal_freqFromPhase,'b')
#    