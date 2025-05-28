      subroutine ini_optcld(QEXTVc,QSCATVc,GrefVc,QEXTIc,QSCATIc,GrefIc,
     *                  Qxvc,Qxic,Qsvc,Qsic,gvc,gic,Qextrefc)

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!  Initialize cloud optical constants, both IR and visible

      use grid_h
      use radinc_h
      use cldcommon_h

      implicit none

!  Arguments
!  ---------

      real*8  QEXTVc(L_NSPECTV)
      real*8  QSCATVc(L_NSPECTV)
      real*8  GrefVc(L_NSPECTV)

      real*8  QEXTIc(L_NSPECTI)
      real*8  QSCATIc(L_NSPECTI)
      real*8  GrefIc(L_NSPECTI)

      real*8  Qxvc(L_LEVELS+1,L_NSPECTV)
      real*8  Qsvc(L_LEVELS+1,L_NSPECTV)
      real*8  gvc(L_LEVELS+1,L_NSPECTV)

      real*8  Qxic(L_LEVELS+1,L_NSPECTI)
      real*8  Qsic(L_LEVELS+1,L_NSPECTI)
      real*8  gic(L_LEVELS+1,L_NSPECTI)

      real*8  Qextrefc(L_LEVELS+1)

!  Local variables
!  ---------------

      integer i,k

! Initialyze various variables
! ----------------------------

      DO K = 1, L_LEVELS+1
!        Qextrefc(K) = QEXTVc(L_NREFV)
        Qextrefc(K) = 1. !QEXTVc(L_NREFV)
      ENDDO

      do i = 1, nlonv
        DO K = 1, L_LEVELS+1
!          Qxvc(K,i) = QEXTVc(i)
!          Qsvc(K,i) = QSCATVc(i)
!          gvc(K,i)  = GrefVc(i)
          Qxvc(K,i) = 1.0 !QEXTVc(i)
          Qsvc(K,i) = 0.99 !QSCATVc(i)
          gvc(K,i)  = 0. !GrefVc(i)
        ENDDO
      enddo
      do i = 1, nloni
        DO K = 1, L_LEVELS+1
!          Qxic(K,i) = QEXTIc(i)
!          Qsic(K,i) = QSCATIc(i)
!          gic(K,i)  = GrefIc(i)
          Qxic(K,i) = 1.0 !QEXTIc(i)
          Qsic(K,i) = 0.99 !QSCATIc(i)
          gic(K,i)  = 0. !GrefIc(i)
        ENDDO
      enddo

      return

      end
