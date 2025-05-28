      subroutine setspi(WNOI,DWNI,WAVEI)

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!     PURPOSE:
!        Set up the spectral intervals in the infrared.  Based on
!     Chris McKay's SETSPI code.
!
!     INPUT PARAMETERS
!     L_NSPECTI  - Number of spectral intervals in the INFRARED
!
!     OUTPUT PARAMETERS
!     WNOI       - Array of wavenumbers at the spectral interval
!                  centers for the infrared.  Array is NSPECTI
!                  elements long.
!     DWNI       - Array of "delta wavenumber", i.e., the width,
!                  in wavenumbers (cm^-1) of each IR spectral
!                  interval.  NSPECTI elements long.
!     WAVEI      - Array (NSPECTI elements long) of the wavelenght
!                  (in microns) at the center of each IR spectral
!                  interval.
!     
!**********************************************************************C

      use grid_h
      use constants_h, only: PI
      use radinc_h
      use radcommon_h, only: planckir

      implicit none

!     BWNI - Bin wavenumber of the edges of the IR spectral bins
!     units are inverse centimeters.  Dimension needs to be changed
!     if the number of IR bins changes.


      REAL*8  :: WNOI(L_NSPECTI), DWNI(L_NSPECTI), WAVEI(L_NSPECTI)

      real*8  :: a, b, ans, y, bpa, bma, T
      real*8  :: wn1, wn2
      integer :: n, nw, nt, m

!  C1 and C2 values from Goody and Yung (2nd edition)  MKS units
!  These values lead to a "sigma" (sigma*T^4) of 5.67032E-8 W m^-2 K^-4

      real*8 :: c1 = 3.741832D-16      ! W m^-2
      real*8 :: c2 = 1.438786D-2       ! m K
      
      real*8 :: x(12) = [ -0.981560634246719D0,  -0.904117256370475D0, &
                          -0.769902674194305D0,  -0.587317954286617D0, &
                          -0.367831498998180D0,  -0.125233408511469D0, &
                           0.125233408511469D0,   0.367831498998180D0, &
                           0.587317954286617D0,   0.769902674194305D0, &
                           0.904117256370475D0,   0.981560634246719D0 ]

      real*8 :: w(12) = [  0.047175336386512D0,   0.106939325995318D0, &
                           0.160078328543346D0,   0.203167426723066D0, &
                           0.233492536538355D0,   0.249147045813403D0, &
                           0.249147045813403D0,   0.233492536538355D0, &
                           0.203167426723066D0,   0.160078328543346D0, &
                           0.106939325995318D0,   0.047175336386512D0 ]

!======================================================================C

!     Bin wavenumber - wavenumber [cm^(-1)] at the edges of the IR
!     spectral bins.
        REAL*8 :: BWNI(L_NSPECTI+1)=[                                &
        10.000D0, 34.167D0, 58.333D0, 82.500D0,                      &
        106.667D0, 130.833D0, 155.000D0, 179.167D0,                  &
        203.333D0, 227.500D0, 251.667D0, 275.833D0,                  &
        300.000D0, 320.788D0, 341.575D0, 362.363D0,                  &
        383.150D0, 403.938D0, 424.725D0, 445.513D0,                  &          
        466.301D0, 487.088D0, 507.876D0, 528.663D0,                  &
        549.451D0, 556.742D0, 564.033D0, 571.324D0,                  &
        578.615D0, 585.906D0, 593.197D0, 600.488D0,                  &
        607.779D0, 615.070D0, 622.361D0, 629.652D0,                  &
        636.943D0, 642.550D0, 648.157D0, 653.764D0,                  &
        659.370D0, 664.977D0, 670.584D0, 676.191D0,                  &
        681.798D0, 687.405D0, 693.011D0, 698.618D0,                  &
        704.225D0, 710.139D0, 716.053D0, 721.967D0,                  &
        727.881D0, 733.795D0, 739.709D0, 745.624D0,                  &
        751.538D0, 757.452D0, 763.366D0, 769.280D0,                  &
        775.194D0, 791.501D0, 807.807D0, 824.114D0,                  &
        840.421D0, 856.727D0, 873.034D0, 889.341D0,                  &
        905.647D0, 921.954D0, 938.261D0, 954.567D0,                  &
        970.874D0, 996.805D0, 1022.737D0, 1048.668D0,                &
        1074.600D0, 1100.531D0, 1126.463D0, 1152.394D0,              &
        1178.325D0, 1204.257D0, 1230.188D0, 1256.120D0,              &
        1282.051D0, 1360.399D0, 1438.746D0, 1517.094D0,              &
        1595.441D0, 1673.789D0, 1752.137D0, 1830.484D0,              &
        1908.832D0, 1987.179D0, 2065.527D0, 2143.874D0,              &
        2222.222D0]

!     Set up mean wavenumbers and wavenumber deltas.  Units of 
!     wavenumbers is cm^(-1); units of wavelengths is microns.

      do M=1,L_NSPECTI
        WNOI(M)  = 0.5*(BWNI(M+1)+BWNI(M))
        DWNI(M)  = BWNI(M+1)-BWNI(M)
        WAVEI(M) = 1.0E+4/WNOI(M)
      end do

!  For each IR wavelength interval, compute the integral of B(T), the
!  Planck function, divided by the wavelength interval, in cm-1.  The
!  integration is in MKS units, the final answer is the same as the
!  original planck.f; W m^-2 wavenumber^-1, where wavenumber is in CM^-1.

      DO NW=1,L_NSPECTI
        a = 1.0D-2/BWNI(NW+1)
        b = 1.0D-2/BWNI(NW)
        bpa = (b+a)/2.0
        bma = (b-a)/2.0
        do nt=500,9000
          T   = dble(NT)/1.0D+1
          ans = 0.0D0
          do m=1,12
            y    = bma*x(m)+bpa
            ans  = ans + w(m)*c1/(y**5*(exp(c2/(y*T))-1.0D0))
          end do
          planckir(NW,nt-499) = ans*bma/(PI*DWNI(NW))
        end do
      END DO

      return
      end subroutine setspi
