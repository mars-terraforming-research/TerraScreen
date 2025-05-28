      SUBROUTINE FILLTAUCUM

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!
!  PURPOSE
!      FILLTAUCUM fills the TAUCUM array - the dust column density of each
!      GCM level.  NDUSK(K) is the column density of the sublayer whose
!      bottom is level K, and whose top is at K-1.
!
      use grid_h
      use defines_h
      use radinc_h
      use radcommon_h
      use standard_h
      use comp3cmn_h

      implicit none

      external JSRCHGT

      INTEGER Nstar, N

!     implicit none

      integer :: J, I, K, JSRCHGT
      real*8  :: PSTAR, PSTAR1
 
!#=====================================================================

!  Dust optical depth at the bottom of each sub-layer.

      DO N=1,3
        TAUCUM(N) = 0.0
        TAUREF(N) = 0.0
      END DO

      DO K=4,L_LEVELS
        PSTAR     = PSF*SIGMA(K)+PTROP
        PSTAR1    = MAX(PSTAR,PRDST(1))
        NSTAR     = MIN0(JSRCHGT(NPDST-1,PRDST,1,PSTAR1)-1,NPDST-1)
        TAUCUM(K) = TAUDST(NSTAR)+(PSTAR1-PRDST(NSTAR))*
     *              (TAUDST(NSTAR+1) - TAUDST(NSTAR))/
     *              (PRDST(NSTAR+1)-PRDST(NSTAR))
        TAUREF(K) = TAUCUM(K) - TAUCUM(K-1)
      END DO

      TAUREF(L_LEVELS+1) = 0.0

      RETURN
      END
