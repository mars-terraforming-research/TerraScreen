import numpy as np
def read_output(filename,tau_target=None):
    '''
    Unpack radiation driver
    '''
    fname='/home/akling/Code/TerraScreen/output/'+filename
    #===Read header===
    Nlev=24
    L_NSPECTI=96
    L_NSPECTV=84
    f = open(fname,'r')
    part_name=f.readline().split(',')[1]
    Ncase=int(f.readline().split(',')[1])
    dt=float(f.readline().split(',')[1])
    Qext670=float(f.readline().split(',')[1])
    alb_sfc=float(f.readline().split(',')[1])
    conrath_nu=float(f.readline().split(',')[1])
    solar_flux=float(f.readline().split(',')[1])
    bwni=np.array([float(x) for x in f.readline().split(',')[1:]])
    bwnv=np.array([float(x) for x in f.readline().split(',')[1:]])
    solar_WN=np.array([float(x) for x in f.readline().split(',')[1:]])
    f.readline() #skip ====
    wl=np.array([float(x) for x in f.readline().split(',')[1:]])
    Qext=np.array([float(x) for x in f.readline().split(',')[1:]])
    Qscat=np.array([float(x) for x in f.readline().split(',')[1:]])
    g=np.array([float(x) for x in f.readline().split(',')[1:]])
    f.readline() #skip ====
    Pref=np.array([float(x) for x in f.readline().split(',')[1:]])
    Tref=np.array([float(x) for x in f.readline().split(',')[1:]])
    [f.readline() for _ in range(2)]#skip 2 line
    #===
    T=np.zeros((Ncase,Nlev))
    OLR_WL,SFC_dIR= [np.zeros((Ncase,L_NSPECTI)) for _ in range(2)]
    ASR_WL,SFC_dVIS=[np.zeros((Ncase,L_NSPECTV)) for _ in range(2)]
    tau,ts,net_top,net_bot,alb,OLR,ASR=[np.zeros((Ncase)) for _ in range(7)]
    for i in range(Ncase):
        vals=np.array([float(x) for x in f.readline().split(',')])
        tau[i]=vals[1]
        ts[i]=vals[2]
        alb[i]=vals[3]
        OLR[i]=vals[4]
        ASR[i]=vals[5]
        net_top[i]=vals[6]
        net_bot[i]=vals[7]
        T[i,:]=vals[8:8+Nlev]
        OLR_WL[i,:]=vals[9+Nlev:9+Nlev+L_NSPECTI]
        ASR_WL[i,:]=vals[9+Nlev+L_NSPECTI:9+Nlev+L_NSPECTI+L_NSPECTV]
        SFC_dIR[i,:]=vals[9+Nlev+L_NSPECTI+L_NSPECTV:9+Nlev+2*L_NSPECTI+L_NSPECTV]
        SFC_dVIS[i,:]=vals[9+Nlev+2*L_NSPECTI+L_NSPECTV:]
    f.close()

    '''
    #Read Extinction, Scattering
    fname='/home/akling/Code/rtdriver_RCE/data/QEXT_96IR_84_VIS_'+filename[7:-4]
    f = open(fname,'r')
    Qext,Qscat,g,wl=[np.zeros((L_NSPECTI+L_NSPECTV)) for _ in range(4)]
    [f.readline() for _ in range(3)]#skip 3 lines
    for i in range(L_NSPECTV):
        vals=np.array([float(x) for x in f.readline().strip().split()])
        Qext[L_NSPECTV-1-i],Qscat[L_NSPECTV-1-i],g[L_NSPECTV-1-i],BL,BR=vals[1],vals[2],vals[4],vals[5],vals[6]
        wl[L_NSPECTV-1-i]=0.5*(BL+BR)
    [f.readline() for _ in range(3)]#skip 3 lines
    N=L_NSPECTV+L_NSPECTI-1
    for i in range(L_NSPECTI):
        vals=np.array([float(x) for x in f.readline().strip().split()])
        Qext[N-i],Qscat[N-i],g[N-i],BL,BR=vals[1],vals[2],vals[4],vals[5],vals[6]
        wl[N-i]=0.5*(BL+BR)
    '''

    #Compute center wavenumber
    res_wn_bwni=bwni[1:]-bwni[0:-1]
    center_wn_bwni=bwni[0:-1]+res_wn_bwni/2
    res_wn_bwnv=bwnv[1:]-bwnv[0:-1]
    center_wn_bwnv=bwnv[0:-1]+res_wn_bwnv/2
    #Initialize model
    class model(object):
        pass
    #Load output in bject
    MOD=model()
    #Store the filename as an attribute
    setattr(MOD,'name',filename)
    setattr(MOD,'Qext670',Qext670)
    setattr(MOD,'alb_sfc',alb_sfc)
    setattr(MOD,'conrath_nu',conrath_nu)
    setattr(MOD,'solar_flux',solar_flux)
    setattr(MOD,'bwni',bwni)
    setattr(MOD,'bwnv',bwnv)
    setattr(MOD,'solar_WN',solar_WN)
    setattr(MOD,'wni',center_wn_bwni)
    setattr(MOD,'wnv',center_wn_bwnv)
    setattr(MOD,'wli',10**4/center_wn_bwni)
    setattr(MOD,'wlv',10**4/center_wn_bwnv)
    setattr(MOD,'Qext',Qext)
    setattr(MOD,'Qscat',Qscat)
    setattr(MOD,'g',g)
    setattr(MOD,'wl',wl)

    #Save variables
    var_list=['tau','ts','alb','OLR','ASR','net_top','net_bot','OLR_WL','ASR_WL','SFC_dIR','SFC_dVIS']
    #Store variables as attributes
    for ivar in var_list:
        if tau_target is None:
            setattr(MOD,ivar,eval(ivar))
        #Extract one specific value
        else:
            itau=np.argmin(np.abs(tau-tau_target))
            setattr(MOD,ivar,eval('%s[%i,...]'%(ivar,itau)))
    return MOD