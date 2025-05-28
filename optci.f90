      SUBROUTINE OPTCI(DTAUI,TAUCUMI,CO2I,PLEV,TLEV,PGASREF,TGASREF,  &
                       QREFV,QXIDST,QSIDST,GIDST,COSBI,WBARI,TAUREF,   &
                       TMID,PMID,TAUGSURF,QH2O,WREFH2O,                &
                       Qextrefcld,TAUREFCLD,Qxicld,Qsicld,gicld)

!     GCM v24   2012
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

! THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE INFRARED
! IT CALCUALTES FOR EACH LAYER, FOR EACH SPECRAL INTERVAL IN THE IR
! LAYER: WBAR, DTAU, COSBAR
! LEVEL: TAU
!
! Qrefv is the extinction coefficient at the reference (visible) 
! wavelength - 0.67 microns.
!
! TAUI(L,LW) is the cumulative optical depth at level L (or alternatively
! at the *bottom* of layer L), LW is the spectral wavelength interval.
!
!     TLEV(L) - Temperature at the layer boundary (i.e. level)
!     PLEV(L) - Pressure at the layer boundary (i.e. level)
!     CO2_KI(NT,NP,NW,NG) - IR CO2 k-coefficients 
!                           CO2_K(temp,Pres,Waveln,gauss)
!                           currently: CO2_K(7,11,5,17)
!
!----------------------------------------------------------------------C

      use grid_h
      use constants_h, only: Cmk
      use radinc_h
!      use radcommon_h, only:  kgbar_tab, khbar_tab, kbbar_tab  

      implicit none

      real*8  :: DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8  :: DTAUKI(L_LEVELS+1,L_NSPECTI,L_NGAUSS)
      real*8  :: TAUI(L_NLEVRAD,L_NSPECTI,L_NGAUSS)
      real*8  :: TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      real*8  :: TAUGAS
      real*8  :: PLEV(L_LEVELS)
      real*8  :: TLEV(L_LEVELS)      
      real*8  :: TMID(L_LEVELS), PMID(L_LEVELS), LPMID(L_LEVELS)
      real*8  :: CO2I(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8  :: TGASREF(L_NTREF), PGASREF(L_NPREF)
      real*8  :: COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8  :: WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)

      integer :: L, NW, NG, K, LK
      integer :: MT(L_LEVELS,2), MP(L_LEVELS,2), MW(L_LEVELS)
      integer :: mtt(2), mpt(2), mn

      real*8  :: ANS, TAUREFL
      real*8  :: TAEROS(L_LEVELS,L_NSPECTI)
      real*8  :: DPR(L_LEVELS), U(L_LEVELS), TAUAC

!     For dust
      real*8  :: QXIDST(L_LEVELS+1,L_NSPECTI)
      real*8  :: QSIDST(L_LEVELS+1,L_NSPECTI)
      real*8  :: GIDST(L_LEVELS+1,L_NSPECTI)
      real*8  :: QREFV(L_LEVELS+1)
      real*8  :: TAUREF(L_LEVELS+1)
      real*8  :: TAUREF_save(L_LEVELS+1)

!     For clouds
      real*8  :: Qxicld(L_LEVELS+1,L_NSPECTI)
      real*8  :: Qsicld(L_LEVELS+1,L_NSPECTI)
      real*8  :: gicld(L_LEVELS+1,L_NSPECTI)
      real*8  :: Qextrefcld(L_LEVELS+1)
      real*8  :: TAUREFCLD(L_LEVELS+1)

      real*8  :: TCLOUD(L_LEVELS,L_NSPECTI),TAUREFLCLD

      real*8  :: TAUREFLK(L_LEVELS+1,L_NSPECTI)
      real*8  :: TAUCLDK(L_LEVELS+1,L_NSPECTI)

! fraction of zeros in each spectral interval, as a function of T, P

      real*8  :: dt, tt
      real*8  :: taugsurf(L_NSPECTI,L_NGAUSS-1)

!  Water mixing ratio variables

      real*8  :: QH2O(L_LEVELS), WREFH2O(L_REFH2O), WRATIO(L_LEVELS)

      logical  :: pinter(L_LEVELS), tinter(L_LEVELS)
      real*8  :: kcoe

! Stuff I added for  CIA calcs

      real*8  :: wnoi(L_NSPECTI)
      real*8  :: dwni(L_NSPECTI)
      real*8  :: dtaucia(L_LEVELS+1,L_NSPECTI)
      real*8  :: dtau_co2cia(L_LEVELS+1,L_NSPECTI)
      real*8  :: dtau_h2n2cia(L_LEVELS+1,L_NSPECTI)
      real*8  :: taucia(L_NSPECTI)
      real*8  :: kgbar,kbbar,khbar
      real*8  :: tmp,pmp
      real*8  :: plength,namgCO2,namgH2
      real*8  :: kcoeff(l_levels,l_nspecti,l_ngauss)
      integer :: j

!! Read in the CIA absorption coefficeints
!
!      open(10,file='kgbar.dat',form='formatted')
!      open(11,file='kbbar.dat',form='formatted')
!      open(12,file='khbar.dat',form='formatted')
!      do nw=1,L_NSPECTI
!      read(10,500) (kgbar_tab(nw,j),j=1,5)
!      read(11,500) (kbbar_tab(nw,j),j=1,5)
!      read(12,500) (khbar_tab(nw,j),j=1,5)
!  500 format(1x,5(1pe13.3))
!      end do
!      close(10)
!      close(11)
!      close(12)

!======================================================================C

      do nw=1,l_nspecti
!      write(6,'("WNOI in OPTCI= ",f10.3)') wnoi(nw)
      taucia(nw)=0.
      end do

      do nw=1,L_NSPECTI
        do ng=1,L_NGAUSS
          DTAUKI(L_LEVELS+1,nw,ng)   = 0.0D0
        end do
        TAUREFLK(L_LEVELS+1,nw) = 0.0D0
        TAUCLDK(L_LEVELS+1,nw)  = 0.0D0
      end do

!  Save old tauref values

      do k=1,L_LEVELS+1
        tauref_save(k) = tauref(k)
      end do

!  Determine the total gas opacity throughout the column, for each
!  spectral interval, NW, and each Gauss point, NG.

      DO NG=1,L_NGAUSS-1
        do NW=1,L_NSPECTI
          TAUGSURF(NW,NG) = 0.0D0
        end do
      end do

!      print*,"Cmk=",Cmk
      do K=2,L_LEVELS
        DPR(k) = PLEV(K)-PLEV(K-1)
        U(k)   = Cmk*DPR(k)
        LPMID(K) = log10(pmid(k)) 

        call tpindex(LPMID(K),TMID(K),QH2O(K),PGASREF,TGASREF,WREFH2O, &
                     MTt,MPt,MW(K),pinter(k),tinter(k))

        do mn=1,2
          mt(k,mn) = mtt(mn)
          mp(k,mn) = mpt(mn)
        end do

        TAUREF(K)    = TAUREF(K)    / Qrefv(K)
        TAUREFCLD(K) = TAUREFCLD(K) / Qextrefcld(K)

        DO NW=1,L_NSPECTI

          TAEROS(K,NW) = TAUREF(K)    * Qxidst(K,NW)
          TCLOUD(K,NW) = TAUREFCLD(K) * Qxicld(K,NW)
     
! ************* Begin CIA Absorption Coefficients ***************

        tmp = .5*(tlev(k)+tlev(k-1))
        pmp = .5*(plev(k)+plev(k-1))

!        call kinter(nw,tmp,kgbar_tab,kgbar)
!        call kinter(nw,tmp,kbbar_tab,kbbar)
!        call kinter(nw,tmp,khbar_tab,khbar)
        
        dtau_co2cia(k,nw) = 0.
        dtau_h2n2cia(k,nw) = 0.
        namgCO2=0.
        namgH2=0.
        if(k.gt.3) then
! plength = path length in meters
        plength=50.8*tmp*(dpr(k)/pmp) !Alex: 50.8= R/(MCo2.g)=8.314/(0.044*3.72) with MCo2 in kg/mol
! namgCO2 = number of amagats of CO2 gas. Pressure units are millibars
        !namgCO2=0.00*(pmp/tmp)*(273.15/1013.25)
        !namgH2=0.00*(pmp/tmp)*(273.15/1013.25)

!        namgCO2=0.95*(pmp/tmp)*(273.15/1013.25) !alex modify here
!        namgH2=0.00*(pmp/tmp)*(273.15/1013.25)  !not implemented yet
!        dtau_co2cia(k,nw)=namgCO2*namgCO2*100.*(kgbar+kbbar)*plength 
        !print*,dtau_co2cia(k,nw)

!        print *,'NW',NW,'k',k,'dtauco2',dtau_co2cia(k,nw),       &
!                namgCO2, kbbar, plength, tmp, pmp
        dtau_h2n2cia(k,nw)=namgH2*namgCO2*100.*khbar*plength
!       
    	dtau_h2n2cia(k,nw)=0.0
        end if
        !dtaucia(k,nw) = dtau_co2cia(k,nw)+dtau_h2n2cia(k,nw)
        dtaucia(k,nw) = 0. !Alex dtau_co2cia(k,nw)+dtau_h2n2cia(k,nw)
        taucia(nw) = taucia(nw) + dtaucia(k,nw)




! ************ End compute CIA/DIMER absorption coeffients *****************
          
        END DO
      end do

      do K=2,L_LEVELS
        do nw=1,L_NSPECTI
          do ng=1,L_NGAUSS-1

!           NOW COMPUTE TAUGAS


            call boxinterp(co2i(mt(k,1),mp(k,1),mw(k),nw,ng),          &
                 co2i(mt(k,2),mp(k,1),mw(k),nw,ng),                    &
                 co2i(mt(k,1),mp(k,2),mw(k),nw,ng),                    &
                 co2i(mt(k,2),mp(k,2),mw(k),nw,ng),TGASREF(mt(k,1)),   &
                 TGASREF(mt(k,2)),PGASREF(MP(k,1)),PGASREF(MP(K,2)),   &
                 tmid(k),lpmid(k),tinter(k),pinter(k),ans)

            TAUGAS          = U(k)*10.0D0**ans

            kcoeff(k,nw,ng) = 10.0D0**ans
            TAUGSURF(NW,NG) = TAUGSURF(NW,NG) + TAUGAS
            DTAUKI(K,nw,ng) = TAUGAS+TAEROS(K,NW)+TCLOUD(K,NW)         &
                              +dtaucia(k,nw)
          end do

!  Now fill in the "clear" part of the spectrum (NG = L_NGAUSS)
!  Which holds continuum opacity only

          NG              = L_NGAUSS
          DTAUKI(K,nw,ng) = TAEROS(K,NW)+TCLOUD(K,NW)+dtaucia(k,nw)
        end do
      end do

!! print out the k coefficients
!      do nw=1,5
!!      do k=2,l_levels
!      k=4
!      write(14,488) k,wnoi(nw),(kcoeff(k,nw,ng),ng=1,l_ngauss)
!  488 format(1x,i3,1xf10.1,1x16(1pe13.3))
!      end do

!  Now the full treatment for the layers, where besides the opacity
!  we need to calculate the scattering albedo and asymmetry factors
!  for each layer

      DO NW=1,L_NSPECTI
        DO K=2,L_LEVELS+1
          TAUREFLK(K,NW) = TAUREF(K)    * QSIDST(K,NW)
          TAUCLDK(K,NW)  = TAUREFCLD(K) * QSICLD(K,NW)
        ENDDO
      ENDDO

      DO NW=1,L_NSPECTI

!  First, the special "clear" channel

        NG = L_NGAUSS

        DO L=1,L_NLAYRAD
          K              = 2*L+1
          DTAUI(L,nw,ng) = DTAUKI(K,NW,NG) + DTAUKI(K+1,NW,NG) + 1.d-50
          if(DTAUI(L,NW,NG) .GT. 1.0E-9) then
            WBARI(L,nw,ng) = (TAUREFLK(K,NW)+ TAUREFLK(K+1,NW) +       &
                              TAUCLDK(K,NW) + TAUCLDK(K+1,NW)   ) /    &
                              DTAUI(L,NW,NG)
          else
            WBARI(L,nw,ng) = 0.0D0
            DTAUI(L,NW,NG) = 1.0E-9
          endif

          TAUAC = TAUREFLK(K,NW)+ TAUREFLK(K+1,NW) + TAUCLDK(K,NW) +   &
                  TAUCLDK(K+1,NW)

          if(TAUAC .GT. 0.0) then
            cosbi(L,NW,NG) = ( GIDST(K,NW)   * TAUREFLK(K,NW) +        &
                               GIDST(K+1,NW) * TAUREFLK(K+1,NW) +      &
                               GICLD(K,NW)   * TAUCLDK(K,NW) +         &
                               GICLD(K+1,NW) * TAUCLDK(K+1,NW) ) /     &
                              (TAUREFLK(K,NW)+ TAUREFLK(K+1,NW) +      &
                               TAUCLDK(K,NW) + TAUCLDK(K+1,NW)   )
          else
            cosbi(L,NW,NG) = 0.0D0
          end if
        END DO

!  . . .Now the other Gauss points, if needed.

        DO NG=1,L_NGAUSS-1

          DO L=1,L_NLAYRAD
            K              = 2*L+1
            DTAUI(L,nw,ng) = DTAUKI(K,NW,NG)+DTAUKI(K+1,NW,NG)+1.d-50
            if(DTAUI(L,NW,NG) .GT. 1.0E-9) then
              WBARI(L,nw,ng) = (TAUREFLK(K,NW)+ TAUREFLK(K+1,NW) +     &
                                TAUCLDK(K,NW) + TAUCLDK(K+1,NW)   ) /  &
                                DTAUI(L,NW,NG)
            else
              WBARI(L,nw,ng) = 0.0D0
              DTAUI(L,NW,NG) = 1.0E-9
            endif

            cosbi(L,NW,NG) = cosbi(L,NW,L_NGAUSS)
          END DO
        END DO

      END DO     ! NW spectral loop

!     TOTAL EXTINCTION OPTICAL DEPTHS

      DO NW=1,L_NSPECTI
        NG = L_NGAUSS
        TAUI(1,NW,NG) = 0.0D0
        DO L=1,L_NLAYRAD
          TAUI(L+1,NW,NG) = TAUI(L,NW,NG)+DTAUI(L,NW,NG)
        END DO

        TAUCUMI(1,NW,NG)=0.0D0
        DO K=2,L_LEVELS
          TAUCUMI(K,NW,NG)=TAUCUMI(K-1,NW,NG)+DTAUKI(K,NW,NG)
        END DO

        DO NG=1,L_NGAUSS-1
          TAUI(1,NW,NG)=0.0D0
          DO L=1,L_NLAYRAD
            TAUI(L+1,NW,NG)=TAUI(L,NW,NG)+DTAUI(L,NW,NG)
          END DO

          TAUCUMI(1,NW,NG)=0.0D0
          DO K=2,L_LEVELS
            TAUCUMI(K,NW,NG)=TAUCUMI(K-1,NW,NG)+DTAUKI(K,NW,NG)
          END DO
        END DO
      END DO

!  Restore old tauref values

      do k=1,L_LEVELS+1
        tauref(k) = tauref_save(k)
      end do

      RETURN
      end subroutine optci
      
      subroutine kinter(nw,temp,kbar_tab,kbar)

!     -------------------------------------------------------------
!     kbar_tab is a 2D array whose rows correspond to wavenumber
!     and whose columns correspond to temperature
!     -------------------------------------------------------------

      implicit none

      integer, parameter :: nwave=8
      integer, parameter :: ntemp=5
      real*8  kbar_tab(nwave,ntemp)
      real*8 kbar
      integer nw,nt

      real*8 t_tab(ntemp)
      real*8 temp
      integer kt
      real*8 logk1,logk2,logk

      data t_tab/150.,200.,250.,300.,350./

! limit the bounds of temperature for now. But need to check this.
      if(temp.lt.150.) temp=150.
      if(temp.gt.350.) temp=350.

      do nt=1,ntemp-1
      if(temp.ge.t_tab(nt).and.temp.le.t_tab(nt+1)) go to 100
      end do
  100 kt=nt

      logk1=log(kbar_tab(nw,kt))
      logk2=log(kbar_tab(nw,kt+1))
      logk=logk1+(temp-t_tab(kt))*(logk2-logk1)/(t_tab(kt+1)-t_tab(kt))
      kbar = exp(logk)
      if(kbar_tab(nw,kt).eq.0.and.kbar_tab(nw,kt+1).eq.0.) kbar=0.

      return
      end
      
