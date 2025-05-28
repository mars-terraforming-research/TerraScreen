      SUBROUTINE SFLUXI(PLEV,TLEV,DTAUI,TAUCUMI,UBARI,RSFI,DWNI,       &
                        COSBI,WBARI,GWEIGHT,NFLUXTOPI,FMNETI,          &
                        fluxupi,fluxdni,FZEROI,TAUGSURF,upfluxi_wl,    &
                        dnfluxi_wl) !added upfluxi_wl, dnfluxi_wl Alex

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

      use grid_h
      use radinc_h
      use radcommon_h, only: planckir, tlimiti

      implicit none

!  PLEV      - Pressures at levels (i.e. layer boundaries and 
!              mid-point)
!  TLEV      - Temperatures at levels
!  DTAUI     - IR opacity of layer L.
!  TAUCUMI   - IR opacity from top of atmosphere to the bottom
!              of level K.
!  RSFI      - Surface albedo in the IR.
!  DWNI      - Width of each IR spectral interval
!  WBARI     - Scattering parameter
!  COSBI     - Scattering asymetry parameter
!  UBARI     - Cosine of the incidence angle
!  GWEIGHT   - Gauss weight
!  NFLUXTOPI - Net flux at the top of the atmosphere 
!  FMNETI    - Net flux at the bottom of layer L
!  FLUXUPI   - Upward flux at layer L
!  FLUXDNI   - Downward flux at layer L
!  FZEROI    - Fraction of zeros (gas opacity k-coefficient off-line
!              computations) in each spectral interval.
!  TAUGSURF  - Gas optical depth (top of atmosphere to the surface).
!
!----------------------------------------------------------------------!

      integer :: NLEVRAD, L, NW, NG, NTS, NTT

      real*8  :: TLEV(L_LEVELS), PLEV(L_LEVELS)
      real*8  :: TAUCUMI(L_LEVELS,L_NSPECTI,L_NGAUSS)
      real*8  :: FMNETI(L_NLAYRAD)
      real*8  :: DWNI(L_NSPECTI)
      real*8  :: DTAUI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8  :: FMUPI(L_NLEVRAD), FMDI(L_NLEVRAD)
      real*8  :: COSBI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8  :: WBARI(L_NLAYRAD,L_NSPECTI,L_NGAUSS)
      real*8  :: GWEIGHT(L_NGAUSS), NFLUXTOPI
      real*8  :: FTOPUP,FLUXUP,FLUXDN,fluxtopi_up,fluxtopi_dn

      real*8  :: UBARI, RSFI, TSURF, BSURF, TTOP, BTOP, TAUTOP
      real*8  :: PLANCK, PLTOP
      real*8  :: fluxupi(L_NLAYRAD), fluxdni(L_NLAYRAD)
      real*8  :: FZEROI(L_NSPECTI)
      real*8  :: taugsurf(L_NSPECTI,L_NGAUSS-1), fzero


      !------------------
      real*8  :: upfluxi_wl(L_NSPECTI,L_NLAYRAD) !added Alex
      real*8  :: dnfluxi_wl(L_NSPECTI,L_NLAYRAD) !added Alex
!======================================================================C
 
      NLEVRAD = L_NLEVRAD
 
!     ZERO THE NET FLUXES
    
      NFLUXTOPI = 0.0
      fluxtopi_up = 0.0
      fluxtopi_dn = 0.0 

      DO L=1,L_NLAYRAD
        FMNETI(L)  = 0.0
        FLUXUPI(L) = 0.0
        FLUXDNI(L) = 0.0
        DO nw=1,L_NSPECTI
         upfluxi_wl(nw,L)=0.
         dnfluxi_wl(nw,L)=0.
        END DO
      END DO
 
!     WE NOW ENTER A MAJOR LOOP OVER SPECTRAL INTERVALS IN THE INFRARED
!     TO CALCULATE THE NET FLUX IN EACH SPECTRAL INTERVAL

      TTOP  = TLEV(2)
      TSURF = TLEV(L_LEVELS)

      NTS   = TSURF*10.0D0-499
      NTT   = TTOP *10.0D0-499

      DO 501 NW=1,L_NSPECTI

!       SURFACE EMISSIONS - INDEPENDENT OF GAUSS POINTS

        BSURF = (1.-RSFI)*PLANCKIR(NW,NTS)
        PLTOP = PLANCKIR(NW,NTT)

!  If FZEROI(NW) = 1, then the k-coefficients are zero - skip to the
!  special Gauss point at the end.
 
        FZERO = FZEROI(NW)
        IF(FZERO.ge.0.99) goto 40


        DO NG=1,L_NGAUSS-1
         
          if(TAUGSURF(NW,NG).lt. TLIMITI) then
            fzero = fzero + (1.0-FZEROI(NW))*GWEIGHT(NG)
            goto 30
          end if

!         SET UP THE UPPER AND LOWER BOUNDARY CONDITIONS ON THE IR
!         CALCULATE THE DOWNWELLING RADIATION AT THE TOP OF THE MODEL
!         OR THE TOP LAYER WILL COOL TO SPACE UNPHYSICALLY
 
          TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
          BTOP   = (1.0-EXP(-TAUTOP/UBARI))*PLTOP
 
!         WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
!         CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
!         WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER 
          
          CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),         &
                      TAUCUMI(1,NW,NG),                                &
                      WBARI(1,NW,NG),COSBI(1,NW,NG),UBARI,RSFI,BTOP,   &
                      BSURF,FTOPUP,FMUPI,FMDI,FLUXUP,FLUXDN)
 
!         NOW CALCULATE THE CUMULATIVE IR NET FLUX

          NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*GWEIGHT(NG)*           &
                                 (1.0-FZEROI(NW))
          fluxtopi_up = fluxtopi_up+FLUXUP*DWNI(NW)*GWEIGHT(NG)*       &
                                 (1.0-FZEROI(NW))
          fluxtopi_dn = fluxtopi_dn+FLUXDN*DWNI(NW)*GWEIGHT(NG)*       &
                                 (1.0-FZEROI(NW))

          DO L=1,L_NLEVRAD-1

!           CORRECT FOR THE WAVENUMBER INTERVALS

            FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*        &
                                    GWEIGHT(NG)*(1.0-FZEROI(NW))
            FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*GWEIGHT(NG)*   &
                                      (1.0-FZEROI(NW))

            FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*GWEIGHT(NG)*    &
                                      (1.0-FZEROI(NW))
                                      
            upfluxi_wl(nw,L) = upfluxi_wl(nw,L)+ FMUPI(L)*DWNI(NW)*GWEIGHT(NG)*&
                                      (1.0-FZEROI(NW))
                                      
            dnfluxi_wl(nw,L) = dnfluxi_wl(nw,L)+ FMDI(L)*DWNI(NW)*GWEIGHT(NG)*&
                                      (1.0-FZEROI(NW))                                      
          END DO

   30     CONTINUE

       END DO       !End NGAUSS LOOP

   40  CONTINUE

!      SPECIAL 17th Gauss point

       NG     = L_NGAUSS

       TAUTOP = DTAUI(1,NW,NG)*PLEV(2)/(PLEV(4)-PLEV(2))
       BTOP   = (1.0-EXP(-TAUTOP/UBARI))*PLTOP

!      WE CAN NOW SOLVE FOR THE COEFFICIENTS OF THE TWO STREAM
!      CALL A SUBROUTINE THAT SOLVES  FOR THE FLUX TERMS
!      WITHIN EACH INTERVAL AT THE MIDPOINT WAVENUMBER 

       CALL GFLUXI(NLEVRAD,TLEV,NW,DWNI(NW),DTAUI(1,NW,NG),            &
                      TAUCUMI(1,NW,NG),                                &
                      WBARI(1,NW,NG),COSBI(1,NW,NG),UBARI,RSFI,BTOP,   &
                      BSURF,FTOPUP,FMUPI,FMDI,FLUXUP,FLUXDN)
 
!      NOW CALCULATE THE CUMULATIVE IR NET FLUX

       NFLUXTOPI = NFLUXTOPI+FTOPUP*DWNI(NW)*FZERO
       fluxtopi_up = fluxtopi_up+FLUXUP*DWNI(NW)*FZERO
       fluxtopi_dn = fluxtopi_dn+FLUXDN*DWNI(NW)*FZERO

       DO L=1,L_NLEVRAD-1

!        CORRECT FOR THE WAVENUMBER INTERVALS

         FMNETI(L)  = FMNETI(L)+(FMUPI(L)-FMDI(L))*DWNI(NW)*FZERO
         FLUXUPI(L) = FLUXUPI(L) + FMUPI(L)*DWNI(NW)*FZERO
         FLUXDNI(L) = FLUXDNI(L) + FMDI(L)*DWNI(NW)*FZERO
         upfluxi_wl(nw,L) = upfluxi_wl(nw,L) + FMUPI(L)*DWNI(NW)*FZERO
         dnfluxi_wl(nw,L) = dnfluxi_wl(nw,L) + FMDI(L)*DWNI(NW)*FZERO
       END DO
  
  501 CONTINUE      !End Spectral Interval LOOP

! *** END OF MAJOR SPECTRAL INTERVAL LOOP IN THE INFRARED****


!      write(53,*) fluxtopi_up, fluxtopi_dn
!      do l=1,l_nlayrad
!        write(53,*) fluxupi(l),fluxdni(l)
!      enddo



      RETURN
      end subroutine sfluxi
