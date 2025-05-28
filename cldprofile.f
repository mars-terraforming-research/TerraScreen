      SUBROUTINE cldprofile(psf,ptrop,nlev,sigma,pcld,tautotcld,
     *                      taurefcld)

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!
!  PURPOSE
!     CLDPROFILE fills TAUREFCLD(K) with the cloud optical depth (at 
!     the reference wavelength) for each layer.
!
!  PSF       - Surface pressure (mbar)
!  PTROP     - Tropopause pressure (mbar)
!  NLEV      - Number of levels
!  SIGMA     - Sigma values at each level
!  PCLD      - Pressure where cloud is located
!  TAUTOTCLD - Total cloud optical depth of the column
!  TAUREFCLD - Cloud optical depth, in each layer, at the reference
!              wavelength.
!  WTC       - Weight of cloud center
!  WTW       - Weight of cloud wing
!
!----------------------------------------------------------------------

      implicit none

      integer :: nlev
      real*8  :: psf, ptrop, pcld, sigma(nlev), pl(nlev)
      real*8  :: TAUREFCLD(nlev+1), tautotcld
      real*8  :: wtc = 0.5, wtw = 0.25

!     implicit none

      integer :: K
 
!#=====================================================================

!  Pressures at layer boundaries and midpoints.

      pl(1) = 0.0D0
      pl(2) = ptrop/2.0D0
      pl(3) = ptrop

      do K=3,nlev
        pl(k) = sigma(k)*(psf-ptrop)+ptrop
      end do

!  Zero out all layers

      do K=1,nlev+1
        taurefcld(k) = 0.0D0
      end do

!  Find layer where cloud is located

      if(pcld.le.pl(5)) then
!       put cloud center in the top layer below the tropopause
        taurefcld(4) = (wtc+wtw)*tautotcld/2.0D0
        taurefcld(5) = taurefcld(4)
        taurefcld(6) = wtw*tautotcld/2.0D0
        taurefcld(7) = taurefcld(6)

      elseif(pcld.gt.pl(nlev-2)) then

!       put cloud in the bottom layer

        taurefcld(nlev)   = (wtc+wtw)*tautotcld/2.0D0
        taurefcld(nlev-1) = taurefcld(nlev)
        taurefcld(nlev-2) = wtw*tautotcld/2.0D0
        taurefcld(nlev-3) = taurefcld(nlev-2)

      else
       
        do k=5,nlev-4
          
          if(pcld.gt.pl(k) .and. pcld.le.pl(k+2)) then

!  center layer

            taurefcld(k+1) = wtc*tautotcld/2.0D0
            taurefcld(k+2) = taurefcld(k+1)

!  layers adjacent to the cloud's central layer (wings)

!  layer above
            taurefcld(k)   = wtw*tautotcld/2.0
            taurefcld(k-1) = taurefcld(k)
!  layer below
            taurefcld(k+3) = taurefcld(k)
            taurefcld(k+4) = taurefcld(k)
            exit
          end if
        end do
      end if

      return
      end
