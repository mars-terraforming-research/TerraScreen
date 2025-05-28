      subroutine dustprofile

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

! Bob's updates 9/17/99 
! Reference the dust optical depth to RPTAU, and modify the way
! the dust mixing ratio is calculated to more accurately reflect the
! pressure-optical depth relationship.
!  Called by EMISS and GMP.

      use grid_h
      use defines_h
      use constants_h, only: PI
      use standard_h

      implicit none

      REAL*8  QRDST(L_NPDST)

!     implicit none

      integer :: i, j, n
      real*8  :: refpr, conr, sum, pave, qrdst0, pdif1, pdif2
      real*8  :: deltals, conrtmp, tetalat

!======================================================================C

! Calculate the Reference Pressure Grid (prdst)

      refpr    = (5.0*psf/ptrop)**(1.0/(float(npdst)-1.0))
      prdst(1) = ptrop
   
      do n=2,npdst
        prdst(n) = refpr*prdst(n-1)
      end do

! Calculate the Mixing Ratio at the Reference Pressure Level

!  GCM1.7  6/28/01   spatially varying dust

      conr    = conrnu

      sum = 0.
      do n = 2,npdst
        if (prdst(n).lt.rptau) then
          pave = 0.5*(prdst(n)+prdst(n-1))
          sum  = sum + 
     &           exp(conr*(1.-(rptau/pave)))*
     &              (prdst(n)-prdst(n-1))
        end if
        if (prdst(n).ge.rptau) go to 10 
      end do

10    continue

      pave = 0.5*(rptau+prdst(n-1))
      sum  = sum + exp(conr*(1.-(rptau/pave)))*(rptau-prdst(n-1))

      qrdst0 = tautot/sum

! Now calculate the mixing ratio at all other levels

      do n=1,npdst-1

! Region 1: Mixing ratio changes continuously through the layer

        if (rptau.gt.prdst(n+1)) then
          pave     = 0.5*(prdst(n+1)+prdst(n))
          qrdst(n) = qrdst0*exp(conr*(1.0-(rptau/pave)))
        end if

! Region 2: Reference pressure level within this layer. 

        if (rptau.le.prdst(n+1).and.rptau.ge.prdst(n)) then
          pave     = 0.5*(prdst(n)+rptau)
          pdif1    = rptau-prdst(n)
          pdif2    = prdst(n+1)-rptau
          qrdst(n) = qrdst0*(
     &               exp(conr*(1.0-(rptau/pave)))*pdif1+pdif2) / 
     &               (prdst(n+1)-prdst(n))
        end if

! Region 3: Mixing ratio constant

        if (rptau.lt.prdst(n)) then
          qrdst(n) = qrdst0
        end if

      end do

! Now compute the optical depths (taudst).

      taudst(1) = 0.0

      do n=2,npdst
        taudst(n) = taudst(n-1) +
     *                   qrdst(n-1)*(prdst(n)-prdst(n-1))
      end do

      RETURN
      END
