      subroutine setrad(TGASREF,PGASREF,CO2V,CO2I,QEXTV,QSCATV,WV,GV, &
                        QEXTI,QSCATI,WI,GI,QEXTVc,QSCATVc,WVc,GVc,     &
                        QEXTIc,QSCATIc,WIc,GIc,FZEROI,FZEROV)

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!     PURPOSE:
!        Set up values used by the radiation code, such as the CO2 gas
!     absorption coefficients.  True constants are defined, and the 
!     time-independent quantities used by the radiation code are 
!     calculated. 
!
!     TGASREF        - Temperatures of the opacity grid
!     PFGASREF       - Pressures (on fine mesh) of the opacity grid
!     CO2V           - Visible CO2 k-coefficients (CO2 gas opacity)
!     CO2I           - IR CO2 k-coefficients (CO2 gas opacity)
!     FZEROV         - Fraction of zeros in the visible (k-coefficient
!                      off-line calculation)
!     FZEROI         - Fraction of zeros in the IR (k-coefficient
!                      off-line calculation)
!
!     AEROSOL RADIATIVE OPTICAL CONSTANTS
!     Values are at the wavelenght interval center
!
!  Dust
!
!     MIE SCATTERING - Size distribution weighted
!     Qextv    - Extinction efficiency - in the visible.
!     Qscatv   - Scattering efficiency - in the visible.
!     WV       - Single scattering albedo - in the visible.
!     GV       - Asymmetry parameter - in the visible.
!
!     Qexti    - Extinction efficiency - in the infrared.
!     Qscati   - Scattering efficiency - in the infrared.
!     WI       - Single scattering albedo - in the infrared.
!     GI       - Asymmetry parameter - in the infrared.
!     
!  Water ice cloud constants
!
!     MIE SCATTERING - Size distribution weighted
!     Qextvc   - Extinction efficiency - in the visible.
!     Qscatvc  - Scattering efficiency - in the visible.
!     WVc      - Single scattering albedo - in the visible.
!     GVc      - Asymmetry parameter - in the visible.
!
!     Qextic   - Extinction efficiency - in the infrared.
!     Qscatic  - Scattering efficiency - in the infrared.
!     WIc      - Single scattering albedo - in the infrared.
!     GIc      - Asymmetry parameter - in the infrared.
!     
!----------------------------------------------------------------------C

      use grid_h
      use radinc_h

      implicit none

      integer :: N, NS, ios, nd

      real*8  :: CO2I(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8  :: CO2V(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTV,L_NGAUSS)
      real*8  :: PFGASREF(L_PINT)
      real*8  :: PGASREF(L_NPREF), TGASREF(L_NTREF)

      real*8  :: qextv(L_NSPECTV), qextvc(L_NSPECTV)
      real*8  :: qscatv(L_NSPECTV), qscatvc(L_NSPECTV)
      real*8  :: wv(L_NSPECTV), wvc(L_NSPECTV)
      real*8  :: gv(L_NSPECTV), gvc(L_NSPECTV)

      real*8  :: qexti(L_NSPECTI), qextic(L_NSPECTI)
      real*8  :: qscati(L_NSPECTI), qscatic(L_NSPECTI)
      real*8  :: wi(L_NSPECTI), wic(L_NSPECTI)
      real*8  :: gi(L_NSPECTI), gic(L_NSPECTI)

!     real*8 kvis, kir
      integer :: nt, np, nw, ng

      real*8  :: fzeroi(L_NSPECTI)
      real*8  :: fzerov(L_NSPECTV)
      

      integer  :: argc, slash_pos


 
!=======================================================================

!!  Read in the dust and water ice optical parameters:  Qext, Qscat,
!!  w0, and g
!
!!  Note that the extension fname_ext is defined in modules.f90

!============Input file=======================
      if (command_argument_count() < 1) then
        write(*,*) "ERROR: Require optical property input file"
        write(*,*) "Usage:  ./driver  optical.dat"
        stop 
      end if
!   Allocate memory
      call get_command_argument(1, length=arglen)
      allocate(character(len=arglen) :: fullpath_Qext)

      call get_command_argument(1, fullpath_Qext)
      write(*,*)'Reading optical properties:',trim(fullpath_Qext)

       open(20,file=fullpath_Qext,   &      
      status='old',iostat=ios)      
      
      
      if(ios.ne.0) then
        write(6,'("setrad.f90:  Could not open optical table file")')
        stop
      end if
      
!--- locate the last path separator (handles / or \) --------------
       slash_pos = max(index(fullpath_Qext,'/', back=.true.),&
        & index(fullpath_Qext,'\', back=.true.))
     
       if (slash_pos <= 0) then             ! no separator found
         dirname_Qext  = '.'
         filename_Qext = fullpath_Qext
       else
         dirname_Qext  = fullpath_Qext(1:slash_pos-1)        ! everything before the last slash
         filename_Qext = fullpath_Qext(slash_pos+1:)         ! everything after
      end if

      
!
!!  Visible

      do n=1,3
        read(20,*) !header
      end do
!
      do n=1,L_NSPECTV 
        read(20,*) nd, Qextv(n), Qscatv(n), wv(n), gv(n)
      end do
!
!!  IR
!
      do n=1,3
        read(20,*) !header
      end do
!
      do n=1,L_NSPECTI 
        read(20,*) nd, Qexti(n), Qscati(n), wi(n), gi(n)
      end do
!
      close(20)
!
!!  Water ice cloud
!
!      open(20,file='data/QEXT_WATER',status='old',iostat=ios)
!      if(ios.ne.0) then
!        write(6,'("setrad.f90:  Could not open dust Qext file")')
!        stop
!      end if
!
!!  Visible
!
!      do n=1,3
!        read(20,*) !header
!      end do
!
!      do n=1,L_NSPECTV 
!        read(20,*) nd, Qextvc(n), Qscatvc(n), wvc(n), gvc(n)
!      end do
!
!!  IR
!
!      do n=1,3
!        read(20,*) !header
!      end do
!
!      do n=1,L_NSPECTI 
!        read(20,*) nd, Qextic(n), Qscatic(n), wic(n), gic(n)
!      end do
!
!      close(20)
!      
!     Set the reference pressure and temperature arrays.  These are
!     the pressures and temperatures at which we have k-coefficients.

      pgasref( 1) = 1.000D-8
      pgasref( 2) = 3.162D-8
      pgasref( 3) = 1.000D-7
      pgasref( 4) = 3.162D-7
      pgasref( 5) = 1.000D-6
      pgasref( 6) = 3.162D-6
      pgasref( 7) = 1.000D-5
      pgasref( 8) = 3.162D-5
      pgasref( 9) = 1.000D-4
      pgasref(10) = 3.162D-4
      pgasref(11) = 1.000D-3
      pgasref(12) = 3.162D-3
      pgasref(13) = 1.000D-2
      pgasref(14) = 3.162D-2
      pgasref(15) = 1.000D-1
      pgasref(16) = 3.162D-1
      pgasref(17) = 1.000D0
      pgasref(18) = 3.162D0
      pgasref(19) = 1.000D+1
      pgasref(20) = 3.162D+1
      pgasref(21) = 1.000D+2
      pgasref(22) = 3.162D+2
      pgasref(23) = 1.000D+3
      pgasref(24) = 3.162D+3
      pgasref(25) = 1.000D+4

!      tgasref( 1) =  50.0D0
!      tgasref( 2) =  75.0D0
!      tgasref( 3) = 100.0D0
!      tgasref( 4) = 125.0D0
!      tgasref( 5) = 150.0D0
!      tgasref( 6) = 175.0D0
!      tgasref( 7) = 200.0D0
!      tgasref( 8) = 225.0D0
!      tgasref( 9) = 250.0D0
!      tgasref(10) = 275.0D0
!      tgasref(11) = 300.0D0
!      tgasref(12) = 325.0D0
!      tgasref(13) = 350.0D0



      tgasref(1)  =  70.0
      tgasref(2)  = 100.0
      tgasref(3)  = 150.0
      tgasref(4)  = 200.0
      tgasref(5)  = 250.0
      tgasref(6)  = 300.0
      tgasref(7)  = 350.0

!      tgasref(8)  = 400.000D0
!      tgasref(9)  = 450.000D0
!      tgasref(10) = 500.000D0
!      tgasref(11) = 550.000D0
!      tgasref(12) = 600.000D0
!      tgasref(13) = 650.000D0
!      tgasref(14) = 700.000D0
!      tgasref(15) = 750.000D0
!      tgasref(16) = 800.000D0

! 
!     Insure that w0 < 1

      DO N=1,L_NSPECTV
        IF(wv(n).ge.0.9999D0) then
          Qscatv(n) = 0.9999D0*Qextv(n)
          wv(n)     = 0.9999D0
        END IF
      END DO

      DO N=1,L_NSPECTI
        IF(wi(n).ge.0.9999D0) then
          Qscati(n) = 0.9999D0*Qexti(n)
          wi(n)     = 0.9999D0
        END IF
      END DO

!  Fill the water ice cloud variables

!     Insure that w0 < 1

      DO N=1,L_NSPECTV
        IF(wvc(n).ge.0.9999D0) then
          Qscatvc(n) = 0.9999D0*Qextvc(n)
          wvc(n)     = 0.9999D0
        END IF
      END DO

      DO N=1,L_NSPECTI
        IF(wic(n).ge.0.9999D0) then
          Qscatic(n) = 0.9999D0*Qextic(n)
          wic(n)     = 0.9999D0
        END IF
      END DO

!     Interpolate CO2 k coefficients to the finer pressure grid.

!      call laginterp(PGASREF,PFGASREF,CO2I,CO2V,FZEROI,FZEROV)

      call initinterp(PGASREF,CO2I,CO2V,FZEROI,FZEROV)

      return
      end subroutine setrad
