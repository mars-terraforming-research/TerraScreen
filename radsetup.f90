      subroutine radsetup

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!     PURPOSE:
!        Bundle the new radiation code setup subroutines and call
!     this one subroutine from main, where the three included files
!     are also listed.  Quantities are passed between this driver
!     and the radiation code via modules (eg radcommon_h).
!
!----------------------------------------------------------------------C

      use grid_h
      use radinc_h
      use radcommon_h

      implicit none

      integer i,j,l,nw

!======================================================================C

      call setspv(WNOV,DWNV,WAVEV,SOLARF,TAURAY)
      call setspi(WNOI,DWNI,WAVEI)
      call setrad(TGASREF,PGASREF,CO2V,CO2I,QEXTV,QSCATV,WV,GV,       &
                        QEXTI,QSCATI,WI,GI,QEXTVc,QSCATVc,WVc,GVc,     &
                        QEXTIc,QSCATIc,WIc,GIc,FZEROI,FZEROV)

!========================Alex commented========
      ! Read in the CIA absorption coefficeints

!          print *, 'INIT CLOUD CIA COEFFS'
!      open(70,file='./data/kgbar_8band.dat',form='formatted')
!      open(71,file='./data/kbbar_8band.dat',form='formatted')
!      open(72,file='./data/khbar_8band.dat',form='formatted')
!      do nw=1,L_NSPECTI
!      read(70,500) (kgbar_tab(nw,j),j=1,5)
!      read(71,500) (kbbar_tab(nw,j),j=1,5)
!      read(72,500) (khbar_tab(nw,j),j=1,5)
!  500 format(1x,5(1pe13.3))
!      end do
!      close(70)
!      close(71)
!      close(72)



      return
      end subroutine radsetup
