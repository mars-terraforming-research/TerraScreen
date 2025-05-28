      SUBROUTINE OPTCV(DTAUV,TAUV,TAUCUMV,CO2V,PLEV,PGASREF,           &
                       TGASREF,QXVDST,QSVDST,GVDST,WBARV,COSBV,        &
                       TAURAY,TAUREF,TMID,PMID,TAUGSURF,QH2O,WREFH2O,  &
                       Qextrefcld,taurefcld,Qxvcld,Qsvcld,gvcld)

!     GCM v24   2012
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

! THIS SUBROUTINE SETS THE OPTICAL CONSTANTS IN THE VISIBLE 
! IT CALCUALTES FOR EACH LAYER, FOR EACH SPECRAL INTERVAL IN THE VISIBLE
! LAYER: WBAR, DTAU, COSBAR
! LEVEL: TAU
!
! TAUV(L,NW,NG) is the cumulative optical depth at the top of radiation code
! layer L. NW is spectral wavelength interval, ng the Gauss point index.
!
!     TLEV(L) - Temperature at the layer boundary
!     PLEV(L) - Pressure at the layer boundary (i.e. level)
!     CO2V(NT,NPS,NW,NG) - Visible CO2 k-coefficients 
!
!----------------------------------------------------------------------C

      use grid_h
      use constants_h, only: Cmk
      use radinc_h

      implicit none

      real*8  :: DTAUV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8  :: DTAUKV(L_LEVELS+1,L_NSPECTV,L_NGAUSS)
      real*8  :: TAUV(L_NLEVRAD,L_NSPECTV,L_NGAUSS)
      real*8  :: TAUCUMV(L_LEVELS,L_NSPECTV,L_NGAUSS)
      real*8  :: PLEV(L_LEVELS)
      real*8  :: TMID(L_LEVELS), PMID(L_LEVELS), LPMID(L_LEVELS)
      real*8  :: CO2V(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTV,L_NGAUSS)
      real*8  :: TGASREF(L_NTREF), PGASREF(L_NPREF)
      real*8  :: COSBV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8  :: WBARV(L_NLAYRAD,L_NSPECTV,L_NGAUSS)
      real*8  :: TAURAY(L_NSPECTV)

!     For dust
      real*8  :: QXVDST(L_LEVELS+1,L_NSPECTV)
      real*8  :: QSVDST(L_LEVELS+1,L_NSPECTV)
      real*8  :: GVDST(L_LEVELS+1,L_NSPECTV)
      real*8  :: Qextref(L_LEVELS+1)
      real*8  :: TAUREF(L_LEVELS+1)
      real*8  :: TAUREF_save(L_LEVELS+1)

!     For clouds
      real*8  :: Qxvcld(L_LEVELS+1,L_NSPECTV)
      real*8  :: Qsvcld(L_LEVELS+1,L_NSPECTV)
      real*8  :: gvcld(L_LEVELS+1,L_NSPECTV)
      real*8  :: Qextrefcld(L_LEVELS+1)
      real*8  :: TAUREFCLD(L_LEVELS+1)

      real*8  :: TCLOUD(L_LEVELS,L_NSPECTV)

      real*8  :: TAUREFLK(L_LEVELS+1,L_NSPECTV)
      real*8  :: TAUCLDK(L_LEVELS+1,L_NSPECTV)

      integer :: L, NW, NG, K, NG1(L_NSPECTV), LK
      integer :: MT(L_LEVELS,2), MP(L_LEVELS,2), MW(L_LEVELS)
      integer :: mtt(2), mpt(2), mn

      real*8  :: ans, TAUGAS
      real*8  :: TRAY(L_LEVELS,L_NSPECTV)
      real*8  :: TAEROS(L_LEVELS,L_NSPECTV)
      real*8  :: DPR(L_LEVELS), U(L_LEVELS)

      real*8  :: taugsurf(L_NSPECTV,L_NGAUSS-1), TRAYAER

!  Reference wavelength is (now) bin #2 - put into qextref
!      real*8 QextREF

!  Water mixing ratio stuff

      real*8  :: QH2O(L_LEVELS), WREFH2O(L_REFH2O), WRATIO(L_LEVELS)
      integer :: nh2o(L_LEVELS)

      logical :: pinter(L_LEVELS), tinter(L_LEVELS)

!======================================================================C

!  Save old tauref values

      do k=1,L_LEVELS+1
        tauref_save(k) = tauref(k)
      end do

!  Determine the total gas opacity throughout the column, for each
!  spectral interval, NW, and each Gauss point, NG.
!  Calculate the continuum opacities, i.e., those that do not depend on
!  NG, the Gauss index.

      DO NG=1,L_NGAUSS-1
        do NW=1,L_NSPECTV
          TAUGSURF(NW,NG) = 0.0D0
        end do
      end do
     
      do K=2,L_LEVELS
        DPR(k)   = PLEV(K)-PLEV(K-1)
        U(k)     = Cmk*DPR(k)
        LPMID(K) = LOG10(PMID(K))
 
        call tpindex(LPMID(K),TMID(K),QH2O(K),PGASREF,TGASREF,WREFH2O, &
                     MTt,MPt,MW(K),pinter(k),tinter(k))

        do mn=1,2
          mt(K,mn) = mtt(mn)
          mp(K,mn) = mpt(mn)
        end do

        Qextref(K) = QXVDST(K,L_NREFV)

        TAUREF(K)    = TAUREF(K) / Qextref(K)
        TAUREFCLD(K) = TAUREFCLD(K) / Qextrefcld(K)

        DO NW=1,L_NSPECTV
          TRAY(K,NW)   = TAURAY(NW)*DPR(K)
          TAEROS(K,NW) = TAUREF(K)    * QXVDST(K,NW)
          TCLOUD(K,NW) = TAUREFCLD(K) * QXVCLD(K,NW)
        END DO
      end do

!  TRAYAER is Tau RAYleigh scattering, plus AERosol opacity

      do K=2,L_LEVELS
 
        do NW=1,L_NSPECTV

          TRAYAER = TRAY(K,NW) + TAEROS(K,NW) + TCLOUD(K,NW)

          do NG=1,L_NGAUSS-1

!           NOW COMPUTE TAUGAS

            call boxinterp(co2v(mt(k,1),mp(k,1),mw(k),nw,ng),          &
                 co2v(mt(k,2),mp(k,1),mw(k),nw,ng),                    &
                 co2v(mt(k,1),mp(k,2),mw(k),nw,ng),                    &
                 co2v(mt(k,2),mp(k,2),mw(k),nw,ng),TGASREF(mt(k,1)),   &
                 TGASREF(mt(k,2)),PGASREF(MP(k,1)),PGASREF(MP(K,2)),   &
                 tmid(k),lpmid(k),tinter(k),pinter(k),ans)

            TAUGAS          = U(k)*10.0D0**ans
            TAUGSURF(NW,NG) = TAUGSURF(NW,NG) + TAUGAS
            DTAUKV(K,nw,ng) = TAUGAS + TRAYAER 
          end do

!  Now fill in the "clear" part of the spectrum (NG = L_NGAUSS)
!  Which holds continuum opacity only

          NG = L_NGAUSS
          DTAUKV(K,nw,ng) = TAEROS(K,NW)+TRAY(K,NW)+TCLOUD(K,NW)

        end do
      end do

!  Now the full treatment for the layers, where besides the opacity
!  we need to calculate the scattering albedo and asymmetry factors
!  for each layer

      DO NW=1,L_NSPECTV
        DO K=2,L_LEVELS
          TAUREFLK(K,NW) = TAUREF(K)    * QSVDST(K,NW)
          TAUCLDK(K,NW)  = TAUREFCLD(K) * QSVCLD(K,NW)
        ENDDO
      ENDDO

      DO NW=1,L_NSPECTV

!  First, the special "clear" channel
 
        NG = L_NGAUSS
        DO L=1,L_LAYERS
          K              = 2*L+1
          DTAUV(L,nw,ng) = DTAUKV(K,NW,NG)+DTAUKV(K+1,NW,NG)
          COSBV(L,NW,NG) = ( GVDST(K,NW)  * TAUREFLK(K,NW) +           &
                             GVDST(K+1,NW)* TAUREFLK(K+1,NW) +         &
                             GVCLD(K,NW)  * TAUCLDK(K,NW)  +           &
                             GVCLD(K+1,NW)* TAUCLDK(K+1,NW) ) /        &
                           ( TRAY(K,NW)     + TRAY(K+1,NW) +           &
                             TAUREFLK(K,NW) + TAUREFLK(K+1,NW) +       &
                             TAUCLDK(K,NW)  + TAUCLDK(K+1,NW) )

          WBARV(L,nw,ng) = ( TAUREFLK(K,NW) + TAUREFLK(K+1,NW) +       &
                             TAUCLDK(K,NW)  + TAUCLDK(K+1,NW)  +       &
                            (TRAY(K,NW)+TRAY(K+1,NW))*0.9999)/         &
                             DTAUV(L,nw,ng)
        END DO

!  Special bottom layer

        L              = L_NLAYRAD
        K              = 2*L+1
        DTAUV(L,nw,ng) = DTAUKV(K,NW,NG)
        COSBV(L,NW,NG) = ( GVDST(K,NW) * TAUREFLK(K,NW) +              &
                           GVCLD(K,NW) * TAUCLDK(K,NW) ) /             &
                         ( TRAY(K,NW)  + TAUREFLK(K,NW) +              &
                           TAUCLDK(K,NW) )

        WBARV(L,nw,ng) = (TAUREFLK(K,NW) + TAUCLDK(K,NW) +             &
                          TRAY(K,NW)*0.9999)/DTAUV(L,nw,ng)

!  . . .Now the other Gauss points, if needed.

        DO NG=1,L_NGAUSS-1
          DO L=1,L_LAYERS
            K              = 2*L+1
            DTAUV(L,nw,ng) = DTAUKV(K,NW,NG)+DTAUKV(K+1,NW,NG)
            COSBV(L,NW,NG) = COSBV(L,NW,L_NGAUSS)
            WBARV(L,nw,ng) = ( TAUREFLK(K,NW) + TAUREFLK(K+1,NW) +     &
                               TAUCLDK(K,NW)  + TAUCLDK(K+1,NW) +      &
                              (TRAY(K,NW)+TRAY(K+1,NW))*0.9999)/       &
                               DTAUV(L,nw,ng)
          END DO

!  Special bottom layer

          L              = L_NLAYRAD
          K              = 2*L+1
          DTAUV(L,nw,ng) = DTAUKV(K,NW,NG)
          COSBV(L,NW,NG) = COSBV(L,NW,L_NGAUSS)
          WBARV(L,nw,ng) = ( TAUREFLK(K,NW) + TAUCLDK(K,NW) +          &
                             TRAY(K,NW)*0.9999 ) / DTAUV(L,nw,ng)
        END DO

      END DO     ! NW spectral loop

!     TOTAL EXTINCTION OPTICAL DEPTHS

      DO NW=1,L_NSPECTV
        NG = L_NGAUSS
        TAUV(1,NW,NG) = 0.0D0
        DO L=1,L_NLAYRAD
          TAUV(L+1,NW,NG) = TAUV(L,NW,NG)+DTAUV(L,NW,NG)
        END DO

        TAUCUMV(1,NW,NG)=0.0D0
        DO K=2,L_LEVELS
          TAUCUMV(K,NW,NG)=TAUCUMV(K-1,NW,NG)+DTAUKV(K,NW,NG)
        END DO
  
        DO NG=1,L_NGAUSS-1
          TAUV(1,NW,NG)=0.0D0
          DO L=1,L_NLAYRAD
            TAUV(L+1,NW,NG)=TAUV(L,NW,NG)+DTAUV(L,NW,NG)
          END DO

          TAUCUMV(1,NW,NG)=0.0D0
          DO K=2,L_LEVELS
            TAUCUMV(K,NW,NG)=TAUCUMV(K-1,NW,NG)+DTAUKV(K,NW,NG)
          END DO
        END DO
      END DO

!  Restore old tauref values

      do k=1,L_LEVELS+1
        tauref(k) = tauref_save(k)
      end do

      RETURN
      end subroutine optcv
