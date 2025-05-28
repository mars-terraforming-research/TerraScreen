      subroutine tpindex(pw,tw,qh2o,pref,tref,WREFH2O,MT,MP,MW,pinter, &
                         tinter)

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center
 
!  PURPOSE
!    Get the TI, UI values for a 2-dimensional interpolation
!    based on the following (The interpolation is done in interpco2):
!    Interpolate the CO2 K-coefficients to the current P,T values.
!    The CO2 coefficients are given on a P,T grid:
!    P = {1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1, 1E+1, 1E+2, 1E+3, 1E+4},
!    T = {50, 100, 150, 200, 250, 300, 350}.
!
!
!  INPUT PARAMETERS
!    PW                 - The (log10) pressure to interpolate to
!    TW                 - The temperature to interpolate to
!    QH2O               - Water abundance
!    WREFH2O            - Array of water mixing ratios of the opacity
!                         grid.
!    
!    
!  OUTPUT PARAMETERS
!    MT                 - Temperature index
!    MP                 - Pressure index
!                         of bounding box
!    MW                 - index of requested water abundance value for
!                         the opacity grid
!----------------------------------------------------------------------C

      use grid_h
      use radinc_h

      implicit none

      integer :: MT(2), MP(2), MW, N
      real*8  :: PW, TW, Qh2o
      real*8  :: WREFH2O(L_REFH2O)
      real*8  :: pref(L_NPREF), tref(L_NTREF) 
      logical :: pinter, tinter

!======================================================================C

!     Get the upper and lower Temperature-grid indicies that bound the
!     requested temperature.  If the requested temperature is outside
!     the T-grid, set up to extrapolate from the appropriate end.

      IF(TW.le.tref(1)) THEN
        MT(1) = 1
        MT(2) = MT(1)
        tinter = .false.
      ELSEIF(TW.GE.tref(L_NTREF)) then
        MT(1) = L_NTREF
        MT(2) = MT(1)
        tinter = .false.
      ELSE
        do n=1,L_NTREF-1
          if(TW.ge.TREF(n) .and. TW.LT.TREF(n+1)) then
            MT(1) = n
            MT(2) = n + 1
            tinter = .true.
            exit
          endif
        end do
      END IF

!     Get the upper and lower Pressure-grid indicies that bound the
!     requested pressure.  If the requested pressure is outside
!     the P-grid, set up to extrapolate from the appropiate end.

      IF(PW.le.Pref(1)) THEN
        MP(1) = 1
        MP(2) = MP(1)
        Pinter = .false.
      ELSEIF(PW.GE.Pref(L_NPREF)) then
        MP(1) = L_NPREF
        MP(2) = MP(1)
        Pinter = .false.
      ELSE
        do n=1,L_NPREF-1 
          if(PW.ge.PREF(n) .and. PW.LT.PREF(n+1)) then
            MP(1) = n
            MP(2) = n + 1
            Pinter = .true.
            exit
          endif
        end do
      END IF

!  Get the indicies for water abundance.  There are 10 sets of 
!  k-coefficients with differing amounts of water vs. CO2.

      IF(Qh2o.ge.WREFH2O(L_REFH2O)) then
        MW = L_REFH2O
      ELSEIF(Qh2o.le.WREFH2O(1)) then
        MW = 1
      ELSE
        DO N=1,L_REFH2O-1
          IF(QH2O.ge.WREFH2O(N) .and. QH2O.lt.WREFH2O(n+1)) then
            MW = N
            EXIT
          END IF
        END DO
      END IF
 
      return
      end subroutine tpindex
