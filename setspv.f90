      SUBROUTINE SETSPV(WNOV,DWNV,WAVEV,SOLARF,TAURAY)

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!     PURPOSE:
!        Set up the spectral intervals in the visible (solar).  Based
!     on Chris McKay's SETSPV code.
!
!     INPUT PARAMETERS
!     L_NSPECTV  - Number of spectral intervals in the visible
!
!     OUTPUT PARAMETERS
!     WNOV       - Array of wavenumbers at the spectral interval
!                  center for the visible.  Array is NSPECTV
!                  elements long.
!     DWNV       - Array of "delta wavenumber", i.e., the width,
!                  in wavenumbers (cm^-1) of each visible spectral
!                  interval.  NSPECTV elements long.
!     WAVEV      - Array (NSPECTV elements long) of the wavelenght
!                  (in microns) at the center of each visible spectral
!                  interval.
!     SOLARF     - Array (NSPECTV elements) of solar flux (W/M^2) in
!                  each spectral interval.  Values are for 1 AU, and
!                  are scaled to the Mars distance elsewhere.
!     TAURAY     - Array (NSPECTV elements) of the wavelength independent
!                  part of Rayleigh Scattering.  The pressure dependent 
!                  part is computed elsewhere (OPTCV).
!
!**********************************************************************C

      use grid_h
      use constants_h, only: GRAV, SCALEP
      use radinc_h

      implicit none

!     BWNV - Bin wavenumber of the edges of the visible spectral bins
!     units are inverse centimeters.  Dimension needs to be changed
!     if the number of visible bins changes.

      REAL*8  :: WNOV(L_NSPECTV), DWNV(L_NSPECTV), WAVEV(L_NSPECTV)
      REAL*8  :: SOLARF(L_NSPECTV), TAURAY(L_NSPECTV)

      REAL*8  :: SUM, WL
      INTEGER :: N, M

!     P0      - Rayleigh scattering reference pressure in pascals.
!     GRAV    - Acceleration due to gravity (g) - MKS

      real*8  :: P0 = 9.423D+6

!     Bin wavenumber - wavenumber [cm^(-1)] at the edges of the visible
!     spectral bins.  Go from smaller to larger wavenumbers, the same as
!     in the IR.

!     2222.22D0    ->   4.50 microns
!     3087.37D0    ->   3.24 microns
!     4030.63D0    ->   2.48 microns
!     5370.57D0    ->   1.86 microns
!     7651.11D0    ->   1.31 microns
!     12500.00D0   ->   0.80 microns
!     25000.00D0   ->   0.40 microns            
!     41666.67D0   ->   0.24 microns

!======================= 7 bands================================
!      REAL*8 :: BWNV(L_NSPECTV+1) = [ 2222.22D0, 3087.37D0, 4030.63D0, &
!                                     5370.57D0, 7651.11D0, 12500.00D0, &
!                                     25000.00D0, 41666.67D0 ]


!     Solar flux within each spectral interval, at 1AU (W/M^2)
!     Sum equals 1356 W/m^2 (values from Wehrli, 1985)

!      real*8 :: SOLAR(L_NSPECTV) = [ 12.7, 24.2, 54.6, 145.9, 354.9,   &
!                                     657.5, 106.3 ]

!======================= 84 bands================================
!     real*8 :: SOLAR(L_NSPECTV) = [0.0000, 0.0000, 1.00000, 1.0000, &
!                                1.00000, 1.00000, 1.00000, 3.00000, &
!                                1.00000, 1.00000, 1.00000, 1.00000, &
!                                2.00000, 1.00000, 2.00000, 2.00000, &
!                                1.00000, 3.00000, 2.00000, 2.00000, &
!                                2.00000, 2.00000, 2.00000, 3.00000, &
!                                3.00000, 4.00000, 4.00000, 3.00000, &
!                                5.00000, 4.00000, 5.00000, 5.00000, &
!                                4.00000, 6.00000, 5.00000, 6.00000, &
!                                 10.0000, 10.0000, 11.0000, 13.0000, &
!                                 11.0000, 13.0000, 12.0000, 12.0000, &
!                                 14.0000, 12.0000, 13.0000, 15.0000, &
!                                 29.0000, 29.0000, 30.0000, 30.0000, &
!                                 29.7000, 30.8000, 30.6000, 27.4000, &
!                                 30.7000, 29.2000, 30.7000, 28.7000, &
!                                 75.3000, 71.7000, 71.0000, 68.8000, &
!                                 63.6000, 59.7000, 52.5000, 48.2000, &
!                                 46.4000, 40.7000, 31.4000, 29.0000, &
!                                 22.4100, 21.2100, 16.1800, 15.0300, &
!                                 10.7100, 7.25000, 6.05500, 2.55200, &
!                                 2.32200, 1.68700, 0.56500, 0.45600]


      real*8 :: SOLAR(L_NSPECTV) = [                                  & 
       0.92112,  0.7875 ,  0.88454,  0.99552,  1.1186 ,  1.03806,     &
        1.1408 ,  1.279  ,  1.1272 ,  1.5518 ,  1.3722 ,  1.5036 ,    &
        1.6576 ,  1.35   ,  1.9596 ,  2.1722 ,  1.7782 ,  1.925  ,    &
        2.0798 ,  2.2586 ,  2.4516 ,  1.7524 ,  2.8072 ,  2.95964,    &
        3.64232,  3.89362,  4.03026,  4.28438,  4.3959 ,  4.5845 ,    &
        4.58102,  5.0122 ,  5.275  ,  5.177  ,  5.65   ,  5.7602 ,    &
        9.878  , 10.7884 , 11.3122 , 11.652  , 11.8336 , 12.878  ,    &
       12.1196 , 12.9144 , 12.4592 , 13.2556 , 12.563  , 13.4482 ,    &
       28.9962 , 28.5648 , 29.2776 , 29.7394 , 30.0156 , 28.8552 ,    &
       29.869  , 30.632  , 31.2554 , 29.0598 , 28.6462 , 29.749  ,    &
       75.208  , 71.668  , 70.878  , 68.7545 , 63.534  , 59.614  ,    &
       52.362  , 48.203  , 46.367  , 40.676  , 31.279  , 29.086  ,    &
       22.461  , 20.6932 , 15.82   , 15.1465 , 10.2489 ,  7.4685 ,    &
        6.0628 ,  2.44658,  2.1696 ,  1.55101,  0.54191,  0.44443]



        REAL*8 :: BWNV(L_NSPECTV+1)=[                                &
        2222.220D0, 2294.316D0, 2366.412D0, 2438.508D0,              &
        2510.603D0, 2582.699D0, 2654.795D0, 2726.891D0,              &
        2798.987D0, 2871.082D0, 2943.178D0, 3015.274D0,              &
        3087.370D0, 3165.975D0, 3244.580D0, 3323.185D0,              &
        3401.790D0, 3480.395D0, 3559.000D0, 3637.605D0,              &
        3716.210D0, 3794.815D0, 3873.420D0, 3952.025D0,              &
        4030.630D0, 4142.292D0, 4253.953D0, 4365.615D0,              &
        4477.277D0, 4588.938D0, 4700.600D0, 4812.262D0,              &
        4923.923D0, 5035.585D0, 5147.247D0, 5258.908D0,              &
        5370.570D0, 5560.615D0, 5750.660D0, 5940.705D0,              &
        6130.750D0, 6320.795D0, 6510.840D0, 6700.885D0,              &
        6890.930D0, 7080.975D0, 7271.020D0, 7461.065D0,              &
        7651.110D0, 8055.184D0, 8459.258D0, 8863.333D0,              &
        9267.407D0, 9671.481D0,   10075.555D0, 10479.629D0,          &
        10883.703D0, 11287.778D0, 11691.852D0, 12095.926D0,          &
        12500.000D0, 13541.667D0, 14583.333D0, 15625.000D0,          &
        16666.667D0, 17708.333D0, 18750.000D0, 19791.667D0,          &
        20833.333D0, 21875.000D0, 22916.667D0, 23958.333D0,          &
        25000.000D0, 26388.889D0, 27777.778D0, 29166.667D0,          &
        30555.557D0, 31944.446D0, 33333.335D0, 34722.224D0,          &
        36111.113D0, 37500.003D0, 38888.892D0, 40277.781D0,          &
        41666.670D0]
!======================================================================C

!     Set up mean wavenumbers and wavenumber deltas.  Units of 
!     wavenumbers is cm^(-1); units of wavelengths is microns.

      do M=1,L_NSPECTV
        WNOV(M)  = 0.5*(BWNV(M+1)+BWNV(M))
        DWNV(M)  = BWNV(M+1)-BWNV(M)
        WAVEV(M) = 1.0E+4/WNOV(M)
      end do

!     Sum the solar flux, and write out the result.  

      sum = 0.0
      do N=1,L_NSPECTV
        SOLARF(N) = SOLAR(N)
        sum       = sum+SOLARF(N)
      end do
      write(6,'("Solar flux at 1AU = ",f7.2," W/M^2")') sum

!     Set up the wavelength independent part of Rayleigh Scattering.
!     The pressure dependent part will be computed elsewhere (OPTCV).
!     WAVEV is in microns.  There is no Rayleigh scattering in the IR.

      do N=1,L_NSPECTV
        WL        = WAVEV(N)
        TAURAY(N) = (8.7/grav)*(1.527*(1.0+0.013/wl**2)/wl**4)*        &
                     scalep/P0
      end do

      RETURN
      end subroutine setspv
