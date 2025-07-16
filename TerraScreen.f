      program driver

!     V23 RT code 2010

      use grid_h
      use defines_h
      use standard_h
      use fccsave_h
      use radinc_h
      use radcommon_h
      use cldcommon_h
      use constants_h

      implicit none

!  PL and TL are the GCM pressures and temperatures at the layer
!  boundaries and midpoints.

      real*8 PL(L_LEVELS), TL(L_LEVELS)

!  PLEV & TLEV are the pressure and temperatures at the GCM levels
!  while PMID & TMID are the pressure and temperatures at the GCM
!  midpoint.

      real*8 PLEV(L_LEVELS), TLEV(L_LEVELS)
      real*8 TMID(L_LEVELS), PMID(L_LEVELS)

!  DIFFVT is the total diffuse visible flux for a given spectral 
!  interval

      real*8 DIFFVT

!  VISUAL

      real*8 DTAUV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8 TAUV(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
      real*8 TAUCUMV(L_LEVELS,L_NSPECTV,L_NGAUSS)
      real*8 COSBV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8 WBARV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8 taugsurf(L_NSPECTV,L_NGAUSS-1)
      real*8 taucump(L_NSPECTV,L_NGAUSS)
      integer ngwv(L_NSPECTV)
      real*8  :: detau(L_NSPECTV,L_NGAUSS)

!  IR

      real*8 DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      real*8 COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8 taugsurfi(L_NSPECTI,L_NGAUSS-1)
      integer ngwi(L_NSPECTI)
!============================
      real*8 tsat 
!  Water mixing

      real*8 QH2O(L_LEVELS)

      real*8 scaleht
      real*8 SOL(L_NSPECTV)

      integer gcmlayers
      integer NLAYRAD, NLEVRAD, NSPECTI, NSPECTV, NPREF, NTREF
      integer K, L, NW, nn, nnn
      real*8  albi, albv
      real*8  acosz
      real*8  ans
      real*8  fluxid(L_NLAYRAD),fluxvd(L_NLAYRAD)
      real*8  heatingv(L_NLAYRAD), heatingir(L_NLAYRAD)
      real*8  total(L_NLAYRAD)
      real*8  fdmax, fdmin
      real*8  firmax, firmin, fvmax, fvmin, df
      real*4  press1,press2,press3,Scht, glapse, zheight(L_LEVELS) !added Alex
      real*8  tstratd, gtd

      real*8 FMNETI(L_NLAYRAD), FMNETV(L_NLAYRAD)
      real*8 fluxupi(L_NLAYRAD), fluxdni(L_NLAYRAD), NFLUXTOPI
      real*8 fluxupv(L_NLAYRAD), fluxdnv(L_NLAYRAD), NFLUXTOPV
      real*8 fluxdv(L_NLAYRAD), fluxdi(L_NLAYRAD)
      integer :: J, I

!  Clouds

      integer :: nlev
      real*8  :: pcld, tautotcld

!=====Alex
      real*8 OM(L_LEVELS),YM(L_LEVELS),PLOGADJ(L_LEVELS),TETA(L_LEVELS)
      real*8 TECON,PCON
      real*8  total_levs(L_LEVELS)
      real*8 upfluxi_wl(L_NSPECTI,L_NLAYRAD) !added Alex
      real*8 dnfluxi_wl(L_NSPECTI,L_NLAYRAD) !added Alex
      real*8 upfluxv_wl(L_NSPECTV,L_NLAYRAD) !added Alex
      real*8 dnfluxv_wl(L_NSPECTV,L_NLAYRAD) !added Alex
      real*8  :: taucum_WN(L_NSPECTV)  !added Alex
      real*8  ::taugastot
      real*8  :: taucum_vis_sfc  !added Alex
      real*8  :: aa,bb
      integer ::ii !loop
      integer ::i_tau !size of opacity to try

      real*8 global_avg !Added alex
      integer ::it !loop
      integer ::n_iter
      parameter (n_iter=15000) 
      real*8 TL_old(L_LEVELS)
      real*8 net_top
      real*8 net_bot
      real*8 gtd_old
      real*8 dt
      real*8 cs !soil heat capacity 
      real*8, parameter :: KAPA = 0.25684D0
      real*8, parameter :: eps_conv = 0.99 !W/m2 convergence
      parameter (i_tau=14) !14
      real*8 tau_array(i_tau)
      real*8 alb_array(i_tau)
      logical equilibrate
      integer ::nmax 
      nmax =n_iter
      tau_array = (/0.,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.75,1.,1.25,
     * 1.5,1.75,2./)
     
!      tau_array = (/0.,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.,0.0,
!     * 0.0,0.0,0.0/)
!      tau_array = (/0.,0.05/)


!      alb_array = (/0.2194,0.1783,0.1459,0.1199,0.0989,0.0677,0.0468,
!     * 0.0326,0.0139,0.0067,0.0038,0.0026, 0.0021,0.0019/)


!======================================================================C

      NLAYRAD = L_NLAYRAD
      NLEVRAD = L_NLEVRAD
      NSPECTI = L_NSPECTI
      NSPECTV = L_NSPECTV
      NPREF   = L_NPREF
      NTREF   = L_NTREF
      NLEV    = L_LEVELS
      
      
!  ALBI is the IR surface albedo, ALBV the surface albedo in the
!  visible.

      ALBI    = 0.00
      ALBV    = 0.22 !0.22 !inititially 0.24
      equilibrate =.false.
      
!     Pseudo surface heat capacity cs for a soil depth of 1m
!     cs dTsfc/dt = - NET_sfc - F_conv_sfc
!     Exact value for cs is not important when running to equilibrium
!     cs = L [m] x cp_soil [J/(kg.K)]  x rho_regolith [kg/m3]   
      cs =  1.*735.9*1481 ![J/(m2 K)] 
!     Set up spectral intervals in the Solar (VISUAL) and then IR.
!     Read in the k-coefficients. . .




      call radsetup
      
      if (equilibrate) then
        open(69,file='output/output_'//trim(filename_Qext(18:))//'.txt')
      else 
        open(69,file='output/static_'//trim(filename_Qext(18:))//'.txt')
        write(*,*) 'Running in static mode'
         nmax=1
      end if
       

      call ini_optdst(QEXTV,QSCATV,GV,QEXTI,QSCATI,GI,
     *                QXVDST,QXIDST,QSVDST,QSIDST,GVDST,GIDST,
     *                QEXTREFDST)

      call ini_optcld(QEXTVc,QSCATVc,GVc,QEXTIc,QSCATIc,GIc,
     *                QXVCLD,QXICLD,QSVCLD,QSICLD,GVCLD,GICLD,
     *                QEXTREFCLD)
     

      
      
      
      



!==============Setup temperature profile===============
      psf       = 6.52
      gtd     =  249.1   !218.
      ptrop     = 0.1
      CONRNU    = 0.002
      acosz     = 0.5 
      tautotcld = 0.0
      pcld      = 1.
      glapse = 5. ! K/km ! for 250K
!     glapse = 1000.*3.72/736 ! K/km ! for 250K and Wordsworth
      Scht = 9.6 ! km
      dt   = 1./48. ! 1/48 time step in sols



!===============QH2O is the mixing ratio of water=========

      do k=1,L_LEVELS
        QH2O(K) = 1.0D-7
      end do




!      do k=2,L_LEVELS
!      qh2o(k) = 0.0*6.0e-8
!      if(tl(k).gt.167.) then
!      qh2o(k)=0.*6.11*exp(22.5*(1.0-273.16/tl(k)))/pl(k)
!      endif


      do ii=1,i_tau
     
      TAUTOT    = tau_array(ii)


!---------------------------------------------
      

!  Make the optical depth reference pressure equal to the surface 
!  pressure.

      rptau = psf



!     rsdist = 2.428        ! Ls =   0   SCOSZ = 557
!     rsdist = 2.745        ! Ls =  90   SCOSZ = 493
!     rsdist = 2.147        ! Ls = 180   SCOSZ = 630
!     rsdist = 1.927        ! Ls = 270   SCOSZ = 702
!     rsdist = 2.255        !SOLAR FLUX AT MARS:   601.330
                            !used for tests with GCM 1-D model

!  RSDIST is the square of the sun-Mars distance, in AU.

      rsdist = 2.318
      global_avg= .50731 ! 1: incident sunc | 1/4 for global average
                         !.50731 tuned to give a netflux of 0 W/m2 
      gcmlayers = L_LAYERS
!============================
!!  TSTRATD is the initial temperature; change as needed
!     
!      tstratd = 200.0
!   
!      do k=1,L_LEVELS
!        tl(K) = tstratd
!      end do
!


!============================
!  Calculate the sigma values

      sigma(3) = 0.0
      do L=1,L_LAYERS
        K = 2*L+3
        sigma(K) = sigma(K-2)+DSIG(L)
      end do

      do K=4,L_LEVELS-1,2
        sigma(K) = 0.5*(SIGMA(K+1)+SIGMA(K-1))
      end do

!     Fill local arrays with pressures at layer boundaries and
!     mid-points.

      PL(2) = PTROP/2.0
      DO K=3,L_LEVELS
        PL(K) = SIGMA(K)*(PSF-PTROP)+PTROP
      END DO
!  Calculate the OM exponenent      
      DO K = 3, L_LEVELS-1
        OM(K) = (PL(K)/PL(L_LEVELS))**KAPA
         PLOGADJ(K) = Log(PL(K+1))-Log(PL(K))
      END DO
      OM(2) = (PL(2)/PL(L_LEVELS))**KAPA
      OM(L_LEVELS) = 1.0
      
      
      YM(2) = PTROP*100./GRAV
      DO  K = 4, L_LEVELS-1, 2
                 YM(K) = (PL(K+1)-PL(K-1))*100./GRAV
      END DO
      
      
        call cldprofile(psf,ptrop,nlev,sigma,pcld,tautotcld,taurefcld)

      tl(L_LEVELS)=gtd
      zheight(L_LEVELS)=0.

      do k=L_LEVELS,2,-1
!        tl(k)=T_AL_profile(pl(k))
        tl(k)=press_to_temp(pl(k))
!        tl(k)=instable_profile(pl(k))
!==========================================     
      end do

      tl(1) = tl(2)
      tstratd = tl(1)
      
!  GTD is the ground temperature
      gtd     = TL(L_LEVELS)
      
      


!     Fill cumulative dust optical depth arrays (cum. dust optical
!     depth from the top of the atmosphere to the bottom of level K).

!  TAUREF is the dust optical depth at the reference wavelength, in
!  each layer.  

      IF(TAUTOT.LE.0.0) THEN
        do K=1,L_LEVELS+1
          tauref(K) = 0.0
        end do
        TAUCUM(L_LEVELS) = 0.0D0
      ELSE
        CALL dustprofile
      END IF

!========================================================
!========================================================
!========================================================
!Initialize variables for convergence       

      net_bot = 999999.
      net_top = 999999.
      
      it=1
      do while ((it .le. nmax) .and. (abs(net_top)>eps_conv 
     *  .or. abs(net_bot)>eps_conv)) 

      
 !  FILLPT is the interface subroutine that takes the P & T values
!  on the GCM grid, and puts them on the RT vertical grid.  The
!  radiation code uses PMID, TMID, PLEV, and TLEV values.     
        
      call fillpt(pl,psf,ptrop,gtd,tstratd,tl,plev,tlev,pmid,
     *                   tmid)
     
!  Fill the TAUREF array, the dust column density for each GCM sub-layer

      if(TAUTOT.gt.0.0) then
        call filltaucum
      end if


!  Fill special bottom radiation level to zero.

      TAUREF(L_LEVELS+1) = 0.0

!     And now back to the regular code. . .

!     Calculate solar flux at the current mars distance

      ans = 0.0
      if(acosz.lt.1.0e-4) then
        do NW=1,L_NSPECTV
          SOL(NW) = 0.0
        end do

      else
        do NW=1,L_NSPECTV
          SOL(nw) = SOLARF(NW)/RSDIST*global_avg
          ans     = ans+sol(NW)
        end do
      end if

!      write(6,'("SOLAR FLUX AT MARS:  ",f8.3)') ANS 
 
!     Set up, and solve for, the solar (visual) fluxes, if the sun
!     is up

!  TAUREFCLD(K) is the cloud optical depth at the reference wavelength.
!  The value is referenced at level K, and measured from the next 
!  higher level.  This is not a cumulative value, but just the optical
!  depth in each sub-layer.

!  If the sun is up, calculate the solar fluxes, else set them to zero.

      if(acosz.ge.1.0e-4) then

c        call cldprofile(psf,ptrop,nlev,sigma,pcld,tautotcld,taurefcld)

!  Calculate the optical depths in each layer, spectral interval,
!  and Gauss point.
         
            
        call OPTCV(DTAUV,TAUV,TAUCUMV,CO2V,PLEV,PGASREF,
     *             TGASREF,QXVDST,QSVDST,GVDST,WBARV,COSBV,
     *             TAURAY,TAUREF,TMID,PMID,TAUGSURF,QH2O,WREFH2O,
     *             QEXTREFCLD,TAUREFCLD,QXVCLD,QSVCLD,GVCLD)


      !Alex: dtau is a function of N_GAUSS, values are summmed with GWEIGHT into sfluxv
      
     
        call SFLUXV(DTAUV,TAUV,TAUCUMV,ALBV,WBARV,COSBV,
     *              ACOSZ,SOL,GWEIGHT,NFLUXTOPV,FMNETV,
     *              FLUXUPV,FLUXDNV,DIFFVT,FZEROV,taugsurf,
     *              detau,upfluxv_wl,dnfluxv_wl) !dnfluxv_wl upfluxv_wl added alex


  

!------------------- compute visible gas opacties and total--------
!       taucum_WN(:)=0. 
!       taugastot=0. 
!       do j=1, L_NSPECTV !loop over spectral intervals
!         do i=1,L_NGAUSS !loop over gauss point
!          !wavelength dependent
!          taucum_WN(j)=taucum_WN(j)+TAUGSURF(j,i)*GWEIGHT(i)*        
!     *                                    (1.0-FZEROV(j))
!          !grand total integrated over wavelength AND Gauss points
!          taugastot=taugastot+TAUGSURF(j,i)*GWEIGHT(i)*        
!     *                                    (1.0-FZEROV(j))
!
!         end do
!       end do

!---------------------------  

   

      else
        write(*,*)"acosz. is not <.1.0e-4"
        NFLUXTOPV = 0.0
        do L=1,L_NLAYRAD
          FMNETV(L)  = 0.0
          FLUXUPV(L) = 0.0
          FLUXDNV(L) = 0.0
        end do
      end if



!     Set up, and solve for, the Infrared fluxes

!  Calculate the optical depths in each layer, spectral interval,
!  and Gauss point.


!------------- original-------------------
!      call OPTCI(DTAUI,TAUCUMI,CO2I,PLEV,PGASREF,TGASREF,
!     *           QEXTREFDST,QXIDST,QSIDST,GIDST,COSBI,WBARI,TAUREF,
!     *           TMID,PMID,TAUGSURFI,QH2O,WREFH2O,
!     *           QEXTREFCLD,TAUREFCLD,QXICLD,QSICLD,GICLD)

!--------------CIA-----------------------------
      call OPTCI(DTAUI,TAUCUMI,CO2I,PLEV,TLEV,PGASREF,TGASREF, !added TLEV
     *           QEXTREFDST,QXIDST,QSIDST,GIDST,COSBI,WBARI,TAUREF,
     *           TMID,PMID,TAUGSURFI,QH2O,WREFH2O,
     *           QEXTREFCLD,TAUREFCLD,QXICLD,QSICLD,GICLD)

!  Compute the IR fluxes

      call SFLUXI(PLEV,TLEV,DTAUI,TAUCUMI,UBARI,ALBI,DWNI,
     *     COSBI,WBARI,GWEIGHT,NFLUXTOPI,FMNETI,
     *     fluxupi,fluxdni,FZEROI,taugsurfI,upfluxi_wl,dnfluxi_wl) !added upfluxi_wl Alex

!  Fluxes have been computed.  Below is code that outputs the
!  values for graphics/analysis.


      
!==================
!      do k=2,L_LEVELS
!        qh2o(k) = 0.0*6.0e-8
!        if(tl(k).gt.167.) then
!          qh2o(k)=0.*6.11*exp(22.5*(1.0-273.16/tl(k)))/pl(k)
!        endif
!      end do

!===================================
!     Upward and downward flux

!      write(6,500) fluxupi(1),fluxupi(L_nlayrad)
!  500 format(1x"OLR up at the TOA = ",f10.2,/,
!     *       1x"OLR up at the BOA = ",f10.2)

!      write(6,522) fluxdni(1),fluxdni(L_nlayrad)
!  522 format(1x"IRD down at the TOA = ",f10.2,/,
!     *       1x"IRD down at the BOA = ",f10.2)
        
     
!      write(6,510) fluxupv(1),fluxdnv(1),
!     *             fluxupv(1)/fluxdnv(1),
!     *             fluxdnv(1)-fluxupv(1)
!  510 format(1x"Up vis flux at the TOA = ",f10.2,/,
!     *       1x"Dn vis flux at the TOA = ",f10.2,/,
!     *       1x"Albedo = ",f10.2,/,
!     *       1x"Absorbed Solar = ",f10.2)
    
      aa=0.
      do i=1,L_NSPECTI
        aa=aa+upfluxi_wl(i,1)
      enddo
!      write(*,*)'Up IR flux TOA=',fluxupi(1)
!      write(*,*)'Sum up IR  TOA=',aa

      aa=0.
      bb=0.
      do i=1,L_NSPECTV
        aa=aa+dnfluxv_wl(i,1)
        bb=bb+upfluxv_wl(i,1)
!        write(*,*)DWNV(i),"cm-1, ASR fluxdnv(nw)-fluxupv(nw)=",
!     * dnfluxv_wl(i,1)-upfluxv_wl(i,1)
      enddo
!      write(*,*)'ASR TOA=',fluxdnv(1)-fluxupv(1)
!      write(*,*)'Sum ASR TOA=',aa-bb
!      write(*,*)'Sum dn TOA=',aa
!      write(*,*)'Sum up TOA=',bb



!     Upward and downward flux

      firmax = FLUXUPI(1)
      firmin = FLUXUPI(1)
      fvmax  = FLUXUPV(1)
      fvmin  = FLUXUPV(1)
     
      do L=1,L_NLAYRAD
        if(FLUXUPI(L).gt.firmax) firmax = FLUXUPI(L) 
        if(FLUXUPI(L).lt.firmin) firmin = FLUXUPI(L) 
        if(FLUXUPV(L).gt.fvmax)  fvmax  = FLUXUPV(L) 
        if(FLUXUPV(L).lt.fvmin)  fvmin  = FLUXUPV(L) 
      end do

      df = 0.05*(firmax-firmin)
      firmax = firmax+df
      firmin = firmin-df
      
      df = 0.05*(fvmax-fvmin)
      fvmax = fvmax+df
      fvmin = fvmin-df
      


!     Downward fluxes

      firmax = FLUXDNI(1)
      firmin = FLUXDNI(1)
      fvmax  = FLUXDNV(1)
      fvmin  = FLUXDNV(1)
     
      do L=1,L_NLAYRAD
        if(FLUXDNI(L).gt.firmax) firmax = FLUXDNI(L) 
        if(FLUXDNI(L).lt.firmin) firmin = FLUXDNI(L) 
        if(FLUXDNV(L).gt.fvmax)  fvmax  = FLUXDNV(L) 
        if(FLUXDNV(L).lt.fvmin)  fvmin  = FLUXDNV(L) 
      end do

      df = 0.05*(firmax-firmin)
      firmax = firmax+df
      firmin = firmin-df
      
      df = 0.05*(fvmax-fvmin)
      fvmax = fvmax+df
      fvmin = fvmin-df

!     Net fluxes (as well as T-profile)

      firmax = FMNETI(1)
      firmin = FMNETI(1)
      fvmax  = FMNETV(1)
      fvmin  = FMNETV(1)
     
      do L=1,L_NLAYRAD
        if(FMNETI(L).gt.firmax) firmax = FMNETI(L) 
        if(FMNETI(L).lt.firmin) firmin = FMNETI(L) 
        if(FMNETV(L).gt.fvmax)  fvmax  = FMNETV(L) 
        if(FMNETV(L).lt.fvmin)  fvmin  = FMNETV(L) 
      end do

      df = 0.05*(firmax-firmin)
      firmax = firmax+df
      firmin = firmin-df
      
      df = 0.05*(fvmax-fvmin)
      fvmax = fvmax+df
      fvmin = fvmin-df
      
     

      L = L_NLAYRAD
      scaleht = 0.0

!  Flux divergence

      fluxdv(1) = FMNETV(1)-NFLUXTOPV
      fluxdi(1) = FMNETI(1)-NFLUXTOPI

      do L=2,L_NLAYRAD
        fluxdv(L) = FMNETV(L)-FMNETV(L-1)
        fluxdi(L) = FMNETI(L)-FMNETI(L-1)
      end do

!  Heating rates
     
      heatingv(1)  = (FMNETV(1)-NFLUXTOPV)*88775.0*grav/
     *                      (cp*scalep*PLEV(3))
      heatingir(1) = (FMNETI(1)-NFLUXTOPI)*88775.0*grav/
     *                      (cp*scalep*PLEV(3))
      total(1) = heatingv(1) + heatingir(1)

      fdmax = -10000.0
      fdmin =  10000.0
      
!===========comp3.f=====      
!       DO L=2,L_NLAYRAD
!           IRTOTAL(2*L+1) = FMNETI(L) - FMNETI(L-1)
!       END DO

      do L=2,L_NLAYRAD
    
        heatingv(L)   = (FMNETV(L)-FMNETV(L-1))*88775.0*grav/(cp*scalep*
     *                  (PLEV(2*L+1)-PLEV(2*L-1)))
        heatingir(L)  = (FMNETI(L)-FMNETI(L-1))*88775.0*grav/(cp*scalep*
     *                (PLEV(2*L+1)-PLEV(2*L-1)))
        total(L) = heatingv(L) + heatingir(L)
        total_levs(2*L+1)=total(L)
        
        
      
        if(heatingv(L) .GT. fdmax)  fdmax = heatingv(L)
        if(heatingir(L) .GT. fdmax) fdmax = heatingir(L)
        if(heatingv(L) .LT. fdmin)  fdmin = heatingv(L)
        if(heatingir(L) .LT. fdmin) fdmin = heatingir(L)
        
      end do
 

      
      net_top = fluxdnv(1)-fluxupv(1)+fluxdni(1)-fluxupi(1)
      net_bot = fluxdnv(L_nlayrad)-fluxupv(L_nlayrad)+
     *fluxdni(L_nlayrad)-fluxupi(L_nlayrad)



   
      if (equilibrate) then
      DO K=2,L_LEVELS-1,2
         tl(K)=tl(K)+total_levs(K+1)*dt
!  =====Test for saturation=====         
        tsat = 3182.48/(23.3494-LOG(press2))
        if(tl(k).le.tsat) tl(k)=tsat
        if(tl(k).le.155.) tl(k)=155.
!  =====Stabilize heating rates for projection over large timestep=====            
        if (tl(K) .GT. 400.) then
            tl(K)=400.
         end if   
      
      TETA(K) = TL(K)/OM(K)  
      CALL CONVECT(OM,YM,PLOGADJ,TETA,PL,TECON,PCON)   
      TL(K)=TETA(K)*OM(K)
                      
      END DO
      tl(1) = tl(2)  !TOA layer
! ===== Update Sfc temperature =====      
       TL(L_LEVELS)=TL(L_LEVELS)+net_bot*(dt*88775.0)/cs
       gtd=TL(L_LEVELS)
       end if
      
!=======Header==========      
      if ((ii .eq. 1) .and. (it .eq. 1) ) then 
      write(*,*)'    Tau  , Iter, Tsfc [K] ,  NET TOP ,',
     * '  NET BOT'
      write(*,*)'_______________________________________________'
      END IF
!=============debug=========      
      if (floor(float(it+99)/100.) .eq. float(it+99)/100. ) then
      write(*,'(f10.4,a,I05,200(a,f10.4))')tau_array(ii),',',it,',',gtd,
     * ',',net_top,',',net_bot
      end if
      
      it=it+1
      end do !n_iter    
      write(*,*)'_______________________________________________'
      !===================save on model's layers====     

  
      if (ii .eq. 1) then
       write(69,'(a,a)')         'Particle       ,',fullpath_Qext
       write(69,'(a,I04)')       'Ncase          ,',i_tau
       write(69,'(a,f10.6)')     'dt (sols)      ,',dt
       write(69,'(a,f10.6)')     'Qext 0.67um    ,',QEXTV(L_NREFV)
       write(69,'(a,f10.6)')     'Alb sfc        ,',ALBV
       write(69,'(a,f10.6)')     'Conrath nu     ,',CONRNU
       write(69,'(a,f10.6)')     'Sun Flux (W/m2),',fluxdnv(1)
       write(69,'(a,f12.4,200(a,f12.4))')  'BWN IR (cm-1)  ,',
     * WNOI(1)-DWNI(1)/2,(",",WNOI(L)+DWNI(L)/2,L=1,L_NSPECTI)  
       write(69,'(a,f12.4,200(a,f12.4))')  'BWN VIS (cm-1) ,',
     * WNOV(1)-DWNV(1)/2,(",",WNOV(L)+DWNV(L)/2,L=1,L_NSPECTV)  
       write(69,'(a)')'================'     
       write(69,'(a,200(a,f12.6))')'wavelenght [um]',
     * (",",10**4/WNOV(L),L=L_NSPECTV,1,-1),
     * (",",10**4/WNOI(L),L=L_NSPECTI,1,-1)   
       write(69,'(a,200(a,f12.4))')'Qext           ',
     * (",",QEXTV(L),L=L_NSPECTV,1,-1),
     * (",",QEXTI(L),L=L_NSPECTI,1,-1) 
       write(69,'(a,200(a,f12.4))')'Qscat          ',
     * (",",QSCATV(L),L=L_NSPECTV,1,-1),
     * (",",QSCATI(L),L=L_NSPECTI,1,-1) 
       write(69,'(a,200(a,f10.4))')'G factor       ',
     * (",",GV(L),L=L_NSPECTV,1,-1),
     * (",",GI(L),L=L_NSPECTI,1,-1) 
       write(69,'(a)')'================'     
       write(69,'(a,200(a,f10.4))')'P [mbar]       ',
     * (",",plev(2*L),L=2,L_NLAYRAD)
         write(69,'(a,200(a,f10.4))')'Temp initial[K]',
     * (",",tlev(2*L),L=2,L_NLAYRAD)
        write(69,'(a)')'================'
        write(69,'(a,a,a,a,a,a,a,a,400(a,I0.2))')'it   ,',
     *'   tau    ,',
     *'  Tsfc    ,','   alb    ,','   OLR    ,','    ASR   ,',
     *'  NET top ,','  NET bot ',
     *(", T lev  ",L,L=1,L_NLAYRAD),
     *(",  OLR IR",i,i=1,L_NSPECTI),(', ASR VIS',j,j=1,L_NSPECTV),
     *(", SFC DIR",i,i=1,L_NSPECTI),(',SFC DVIS',j,j=1,L_NSPECTV)
     
       end if 

!tau, gts, alb,OLR,ASR,NET_top,NET_down     
!       On all passes, write total heating rates 

      write(69,'(I05,400(a,f10.4))')it,',',
     * TAUTOT,',',gtd,',',fluxupv(1)/fluxdnv(1),',',
     * fluxupi(1),',',fluxdnv(1)-fluxupv(1),',',net_top,',',net_bot,
     * (",",tl(2*L),L=1,L_NLAYRAD),
     * (",",upfluxi_wl(i,1),i=1,L_NSPECTI),
     * (",",dnfluxv_wl(j,1)-upfluxv_wl(j,1),j=1,L_NSPECTV),
     * (",",dnfluxi_wl(j,L_nlayrad),j=1,L_NSPECTI),
     * (",",dnfluxv_wl(j,L_nlayrad),j=1,L_NSPECTV)      

      
     
      end do !itau
      write(*,*)'Completed.'
      
      
      
      CONTAINS
      
       
       FUNCTION T_AL_profile(Pmbar)
       !=return the temperature profile from Ramses's AL case

       real*8, INTENT(IN) :: Pmbar
       !Output:
       real*8 :: T_AL_profile
       real*8 :: P0,T0,p1,p2,P
       P=100.*Pmbar
  
       if (P .ge. 29.64566) then
            p1= 12.36909458
            p2= -0.33054892
            P0=729.398145
            T0=249.1
        else  
            p1=20.93495584
            p2= -9.96216775
            P0=29.64566
            T0=206.09201
        end if 
        T_AL_profile=T0+p1*(log(P)-log(P0))+p2*(log(P)-log(P0))**2   
        if (T_AL_profile .lt. 155.) then 
          T_AL_profile=155.
         end if    


       END FUNCTION T_AL_profile
       
       FUNCTION press_to_temp(Pmbar)
        !Return the altitude [m] as a function of pressure for a
        !quadratic temperature profile::  
        !    T(z) = T0 + gam(z-z0) + a*(z-z0)^2
        ! valid 610 Pa > P >= 1.2415639872674782 Pa
        real*8, INTENT(IN) :: Pmbar
        real*8::press_to_temp
        real*8::P,Z0,P0,T0,g,gam,a,rgas,delta,sq_delta,Temp,Z
        P=Pmbar*100.
        Z0 = 0.
        P0 = 610.
        T0 = 225.9
        gam = -0.00213479
        a = 1.44823D-8
        rgas = 192.
        g = 3.72
        delta = (4*a*T0 - gam**2)
        sq_delta = sqrt(delta)

        Z= (Z0+ sq_delta/(2*a)*             
     * tan(atan(gam/sq_delta)+log(P0/P) * rgas * sq_delta / (2*g))-
     *              gam / (2*a))
        press_to_temp=(T0 + gam*(Z-Z0) + a*(Z-Z0)**2)
 
        END FUNCTION press_to_temp
        
       FUNCTION instable_profile(Pmbar)
!    Return the pressure in Pa in the troposphere as a function of height for a quadratic temperature profile.
    
       real*8, INTENT(IN) :: Pmbar
       real*8::instable_profile
       real*8::P,Z0,P0,T0,g,gam,a,rgas,delta,sq_delta,Temp,Z
       P=Pmbar*100.
       Z0=0
       T0=225.9
       a=1.e-07
       P0=610.
       rgas=192.
       g=3.72
       gam=-0.008
       delta = (4*a*T0 - gam**2)
       sq_delta = sqrt(delta)
       Z= (Z0+ sq_delta/(2*a)*             
     * tan(atan(gam/sq_delta)+log(P0/P) * rgas * sq_delta / (2*g))-
     *              gam / (2*a))
        instable_profile=(T0 + gam*(Z-Z0) + a*(Z-Z0)**2)
      END FUNCTION instable_profile
        
      
      end
