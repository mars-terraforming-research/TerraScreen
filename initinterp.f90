      subroutine initinterp(pgref,co2i,co2v,fzeroi,fzerov)

!     GCM v24   2012
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!----------------------------------------------------------------------!

!  Set up for interpolation (linear in log pressure) of the CO2 
!  k-coefficients in the pressure domain.  Subsequent use of these
!  values will use a simple linear interpolation in pressure.

      use grid_h
      use radinc_h

      implicit none

      integer :: n, nt, np, nh, ng, nw, m, i, k, mn
      real*8  :: co2i8(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8  :: co2v8(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTV,L_NGAUSS)
      real*8  :: pgref(L_NPREF)

      real*8  :: co2i(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8  :: co2v(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTV,L_NGAUSS)

      real*8  :: fzeroi(L_NSPECTI)
      real*8  :: fzerov(L_NSPECTV)
 
      real*8  :: x, xi(3), yi(3), ans
      real*8  :: p

!======================================================================!

!  Take log of the reference pressures

      do n=1,L_NPREF
        pgref(n) = LOG10(PGREF(n))
      end do

!     Get CO2 k coefficients


!      open(20,file='./data/CO2H2O_V_15B_350K_v4', status='old',      &
!      open(20,file='./data/CO2H2O_V_84_box_DRY', status='old',      &
!      open(20,file='./data/CO2H2O_V_7_32GP_DRY', status='old',      &
      open(20,file='./data/CO2H2O_V_84_box_DRY', status='old',      &
              form='unformatted')
      read(20) co2v8
      read(20) fzerov
      close(20)


!      open(20,file='./data/CO2H2O_IR_15B_350K_v4', status='old',      &
!      open(20,file='./data/CO2H2O_IR_8_32GP_DRY', status='old',      &
      open(20,file='./data/CO2H2O_IR_96B_box_DRY', status='old',      &
              form='unformatted')


      read(20) co2i8
      read(20) fzeroi
      close(20)

!  Take Log10 of the values - we interpolate the log10 of the values,
!  not the values themselves.   Smallest value is 1.0E-200.

      do nt=1,L_NTREF
        do np=1,L_NPREF
          do nh=1,L_REFH2O
            do ng = 1,L_NGAUSS

              do nw=1,L_NSPECTV
                if(co2v8(nt,np,nh,nw,ng).gt.1.0d-200) then
                  co2v(nt,np,nh,nw,ng) = log10(co2v8(nt,np,nh,nw,ng))
                else
                  co2v(nt,np,nh,nw,ng) = -200.0
                end if
              end do
  
              do nw=1,L_NSPECTI
                if(co2i8(nt,np,nh,nw,ng).gt.1.0d-200) then
                  co2i(nt,np,nh,nw,ng) = log10(co2i8(nt,np,nh,nw,ng))
                else
                  co2i(nt,np,nh,nw,ng) = -200.0
                end if
              end do
      
            end do
          end do
        end do
      end do

! "zero" out 17th (or 33rd) Gauss point

      do nt=1,L_NTREF
        do nh=1,L_REFH2O
          do nw=1,L_NSPECTV
            do np=1,L_NPREF
              co2v(nt,np,nh,nw,L_NGAUSS) = -200.0D0 
            end do
          end do 
          do nw=1,L_NSPECTI
            do np=1,L_NPREF
              co2i(nt,np,nh,nw,L_NGAUSS) = -200.0D0
            end do
          end do 
        end do
      end do

      return
      end subroutine initinterp
