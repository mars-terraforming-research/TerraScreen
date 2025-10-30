import numpy as np
import matplotlib.pyplot as plt
from pandas import read_csv
import pandas as pd
import scipy.io

'''
This code is used to bin optical cross section (typically provided in um2) to optical efficicencies [adim]
for the 96 IR and 84 VIS bands used in Terrascreen.

Note, one can convert cross secions  to extinction efficiencies with  with Qext =Cext/(pi R^2) and R the radius of particle of equivalent volume.
If Cext is provided in [um2] and R in [m], the conversion is:
Qext =Cext[m]10^-12/(pi R^2)

Required input:
    w:         wavelength in [um]
    qext:      extinction efficiency
    qscat:     scattering efficiency
    g :        G-factor [+/-1]
    fnameOUT:  fullpath of file to write to disk

Less important input:
    V:     volume of particle in m3 if needed to convert from optical cross section to optical efficiencies
    line1 : a description for the particle, added to the file for reference
    line2 : a description for the particle, added to the file for reference
'''

w=np.linspace(0.2,20000) #[in micron]
qext=np.exp(w)  #extinction efficiency
qscat= qext*0.5 #scattering efficiency
g =w*0      # g factor [+-1]

fnameOUT='/home/akling/Code/rtdriver_RCE/data/QEXT_96IR_84_VIS_PARTICULE_NAME'
V=np.pi*((0.08*1.e-6)**2)*9*1.e-6 #Actual volume of particle in [m3]
r=(V*3./(4*np.pi))**(1./3) #m2
line1="Aluminium rod 9um length, 160nm diameter, r_eq = %.3fum \n"%(r*10**6)
line2="Created from file Al_terraforming_nanoparticles_9um.dat\n"

#Selectr ture or false for plotting
plot_wavelength=True

#==================Hard coded values========
#==================Do NOT TOUCH=============
#==========================================
L_NSPECTI=96
L_NSPECTV=84

bwnv=np.array([2222.220, 2294.316, 2366.412, 2438.508,
2510.603, 2582.699, 2654.795, 2726.891,
2798.987, 2871.082, 2943.178, 3015.274,
3087.370, 3165.975, 3244.580, 3323.185,
3401.790, 3480.395, 3559.000, 3637.605,
3716.210, 3794.815, 3873.420, 3952.025,
4030.630, 4142.292, 4253.953, 4365.615,
4477.277, 4588.938, 4700.600, 4812.262,
4923.923, 5035.585, 5147.247, 5258.908,
5370.570, 5560.615, 5750.660, 5940.705,
6130.750, 6320.795, 6510.840, 6700.885,
6890.930, 7080.975, 7271.020, 7461.065,
7651.110, 8055.184, 8459.258, 8863.333,
9267.407, 9671.481,   10075.555, 10479.629,
10883.703, 11287.778, 11691.852, 12095.926,
12500.000, 13541.667, 14583.333, 15625.000,
16666.667, 17708.333, 18750.000, 19791.667,
20833.333, 21875.000, 22916.667, 23958.333,
25000.000, 26388.889, 27777.778, 29166.667,
30555.557, 31944.446, 33333.335, 34722.224,
36111.113, 37500.003, 38888.892, 40277.781,
41666.670])

bwni=np.array([10.000, 34.167, 58.333, 82.500,
        106.667, 130.833, 155.000, 179.167,
        203.333, 227.500, 251.667, 275.833,
        300.000, 320.788, 341.575, 362.363,
        383.150, 403.938, 424.725, 445.513,
        466.301, 487.088, 507.876, 528.663,
        549.451, 556.742, 564.033, 571.324,
        578.615, 585.906, 593.197, 600.488,
        607.779, 615.070, 622.361, 629.652,
        636.943, 642.550, 648.157, 653.764,
        659.370, 664.977, 670.584, 676.191,
        681.798, 687.405, 693.011, 698.618,
        704.225, 710.139, 716.053, 721.967,
        727.881, 733.795, 739.709, 745.624,
        751.538, 757.452, 763.366, 769.280,
        775.194, 791.501, 807.807, 824.114,
        840.421, 856.727, 873.034, 889.341,
        905.647, 921.954, 938.261, 954.567,
        970.874, 996.805, 1022.737, 1048.668,
        1074.600, 1100.531, 1126.463, 1152.394,
        1178.325, 1204.257, 1230.188, 1256.120,
        1282.051, 1360.399, 1438.746, 1517.094,
        1595.441, 1673.789, 1752.137, 1830.484,
        1908.832, 1987.179, 2065.527, 2143.874,
        2222.222])



qextv,qscatv,gv=[np.zeros((L_NSPECTV)) for _ in range(3)]
qexti,qscati,gi=[np.zeros((L_NSPECTI)) for _ in range(3)]

# Writing to a new file or overwriting an existing file

#Convert micron to wavenumber
wn=10**4/w

with open(fnameOUT, "w") as file:
    file.write(line1)
    file.write(line2)
    file.write("Bin  Qext    Qscat     W0      g     L1     L2\n")
    for i in range(0,L_NSPECTV):
        #==========================Binning==============
        ind=np.where((w<= 10**4/bwnv[i]) & (w > 10**4/bwnv[i+1]))
        #Case where no value is returned, use closest index
        if ind[0].size==0:
            midband=0.5*(10**4/bwnv[i]+10**4/bwnv[i+1])
            ind=np.argmin(np.abs(w-midband))
        print('[',bwnv[i],'-',bwnv[i+1],']=',wn[ind].min(),wn[ind].max())
        qextv[i]=np.mean(np.atleast_1d(qext[ind]),axis=0)
        qscatv[i]=np.mean(np.atleast_1d(qscat[ind]),axis=0)
        gv[i]=np.mean(np.atleast_1d(g[ind]),axis=0)



        #==============================================
        file.write('%2d    %.3f   %.3f   %.3f  %.3f  %.2f  %.2f \n'%(i+1,qextv[i],qscatv[i],qscatv[i]/qextv[i],gv[i],10**4/bwnv[i],10**4/bwnv[i+1]))
    file.write("\n")
    file.write("IR: \n")
    file.write("Bin  Qext    Qscat     W0      g        L1       L2\n")
    for i in range(0,L_NSPECTI):
        ind=np.where((w<= 10**4/bwni[i]) & (w > 10**4/bwni[i+1]))
        if ind[0].size==0:
            midband=0.5*(10**4/bwni[i]+10**4/bwni[i+1])
            ind=np.argmin(np.abs(w-midband))
        print('[',bwni[i],'-',bwni[i+1],']=',wn[ind].min(),wn[ind].max())
        qexti[i]=np.mean(np.atleast_1d(qext[ind]),axis=0)
        qscati[i]=np.mean(np.atleast_1d(qscat[ind]),axis=0)
        gi[i]=np.mean(np.atleast_1d(g[ind]),axis=0)

        file.write('%2d    %.3f   %.3f   %.3f  %.3f  %.2f  %.2f \n'%(i+1,qexti[i],qscati[i],qscati[i]/qexti[i],gi[i],10**4/bwni[i],10**4/bwni[i+1]))


print('Processed format ',format)
print('Written ',fnameOUT)


all_wn_bounds=np.append(bwni[0:-1],bwnv)
res_wn=all_wn_bounds[1:]-all_wn_bounds[0:-1]
center_wn=all_wn_bounds[0:-1]+res_wn/2.

allQext_new=np.append(qexti,qextv)
allQscat_new=np.append(qscati,qscatv)
allG_new=np.append(gi,gv)





plt.close('all')
if not plot_wavelength:
    ax=plt.subplot(311)
    ax.semilogx(10**4/w,qext,'.-r',label=line1)
    wavenuber_plot(ax,all_wn_bounds,allQext_new,style='--',col='k',label='Re-binnning 180 bands',sides=True)


    plt.ylabel('Qext')
    plt.xlim([bwni.min(),bwnv.max()])
    plt.grid()
    plt.legend()
    ax=plt.subplot(312)
    ax.semilogx(10**4/w,qscat,'.-r',label='_')
    wavenuber_plot(ax,all_wn_bounds,allQscat_new,style='--',col='k',label='_',sides=True)

    plt.ylabel('Qscat')
    plt.xlim([bwni.min(),bwnv.max()])
    plt.grid()

    ax=plt.subplot(313)

    ax.semilogx(10**4/w,g,'.-r',label='_')
    wavenuber_plot(ax,all_wn_bounds,allG_new,style='--',col='k',label='_',sides=True)
    plt.xlabel('WaveNumber [cm-1]')
    plt.ylabel('G')
    plt.xlim([bwni.min(),bwnv.max()])
    plt.grid()

    #Plot in wavenumbers
else:


    def wavenuber_plot(ax,bounds_wn,var,style='-',col='orange',label='_',sides=True):
        '''
        Plot the spectrum in wavenumber
        '''
        col_default=col
        label_default='_'
        nband=len(var)
        res_wn=bounds_wn[1:]-bounds_wn[0:-1]
        center_wn=bounds_wn[0:-1]+res_wn/2.
        for i in range(0,nband):
            if sides:
                ax.semilogx([center_wn[i]-res_wn[i]/2,center_wn[i]-res_wn[i]/2],[0.,var[i]],style,color=col,lw=0.5)#left
                ax.semilogx([center_wn[i]+res_wn[i]/2,center_wn[i]+res_wn[i]/2],[0.,var[i]],style,color=col,lw=0.5)#right
            if i==nband-1:label_default=label #If this is the last index, add a label
            ax.semilogx([center_wn[i]-res_wn[i]/2,center_wn[i]+res_wn[i]/2],[var[i],var[i]],style,color=col,lw=0.5,label=label_default)



    ax=plt.subplot(311)
    ax.semilogx(w,qext,'.-r',label=line1)
    plt.ylabel('Qext')
    plt.grid()
    plt.xlim([0.1,100])
    plt.legend()

    ax=plt.subplot(312)
    ax.semilogx(w,qscat,'.-r',label='_')
    #ax.semilogx(w,qscat1,'--',color='brown',label='12um')
    #ax.semilogx(w,qscat2,'--',color='grey',label='9um')
    plt.ylabel('Qscat')
    plt.xlim([0.1,100])
    plt.grid()

    ax=plt.subplot(313)

    ax.semilogx(w,g,'.-r',label='_')
    #sax.semilogx(w,g1,'--',color='brown',label='12um')
    #ax.semilogx(w,g2,'--',color='grey',label='9um')
    plt.xlim([0.1,100])
    plt.xlabel('Wavelength [um]')
    plt.ylabel('G')
    plt.grid()

plt.show()














