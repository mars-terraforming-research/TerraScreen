!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            modules.f90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     RT based on GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            version_h   
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module version_h

      character(len=9), parameter :: version = "V23 r0029"

      end module version_h
      module grid_h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                             grid_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     Number of atmospheric layers
!      integer, PARAMETER :: L_LAYERS  = 29
      integer, PARAMETER :: L_LAYERS  = 24

!     Number of atmospheric levels:   2 * L_LAYERS + 3
      integer, PARAMETER :: L_LEVELS  = 2*L_LAYERS+3

!     Number of dust particle sizes
      integer, parameter :: NDP = 2

      end module grid_h
      module defines_h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            defines_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use grid_h
      implicit none

!     L_LAYERS + 1
      integer, PARAMETER :: L_LAYERP1 = L_LAYERS+1

!     L_LAYERS - 1
      integer, PARAMETER :: L_LAYERM1 = L_LAYERS-1

!     L_LEVELS - 1
      integer, PARAMETER :: L_LEVELM1 = L_LEVELS-1

!     L_LEVELS - 2
      integer, PARAMETER :: L_LEVELM2 = L_LEVELS-2

!     L_LEVELS - 3
      integer, PARAMETER :: L_LEVELM3 = L_LEVELS-3

!     L_LEVELS - 4
      integer, PARAMETER :: L_LEVELM4 = L_LEVELS-4

!     Number of fine-mesh layers used in calculating dust profiles
      integer, PARAMETER :: L_NPDST   = 100

      end module defines_h
      module constants_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            constants_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

!     The mathematical constant PI.
      real*8, parameter :: PI = 3.1415926535897932D0

!     Acceleration due to gravity (mks)
      real*8, parameter :: GRAV = 3.72D0

!     Gas constant for mars.
      real*8, parameter :: RGAS = 1.8902D+2

!     Stefan-Boltzmann constant
      real*8, parameter :: STBO = 5.67051D-8

!     Heat capacity (or specific heat) of CO2 gas.
!     ( units of joules per ( kg * degrees kelvin))
      real*8, parameter :: Cp = 7.3594D+2

!     Factor to convert pressures from millibars to Pascals
      real*8, parameter :: SCALEP = 1.00D+2

!     A radiation code conversion factor.
      real*8, parameter :: Cmk = 3.51D+22 

      end module constants_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            radinc_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module radinc_h

      use grid_h, only: L_LAYERS
      implicit none

!======================================================================C
!
!     RADINC.H    RADiation INCludes
!
!     Includes for the radiation code; RADIATION LAYERS, LEVELS,
!     number of spectral intervals. . .
!
!     GCM2.0  Feb 2003
! 
!======================================================================C

!     RADIATION parameters

!     In radiation code, layer 1 corresponds to the stratosphere.  Level
!     1 is the top of the stratosphere.  The dummy layer is at the same
!     temperature as the (vertically isothermal) stratosphere, and
!     any time it is explicitly needed, the appropriate quantities will
!     be dealt with (aka "top". . .)

!     L_NLEVRAD corresponds to the surface - i.e., the GCM Level that
!     is at the surface.  PLEV(L_NLEVRAD) = P(J,I)+PTROP, 
!     PLEV(2) = PTROP, PLEV(1) = ptrop

!     L_NLAYRAD is the number of radiation code layers
!     L_NLEVRAD is the number of radiation code levels.  Level N is the
!               top of layer N. 
!
!     L_NSPECTI is the number of IR spectral intervals
!     L_NSPECTV is the number of visible (or Solar) spectral intervals
!     L_NGAUSS  is the number of Gauss points for K-coefficients
!               GAUSS POINT 9 (aka the last one) is the special case
!     L_NWNGI   is L_NSPECTI*L_NGAUSS;  the total number of "intervals"
!               in the IR
!     L_NWNGV   is L_NSPECTV*L_NGAUSS;  the total number of "intervals"
!               in the VISIBLE
!
!     L_NPREF   is the number of reference pressures that the 
!               k-coefficients are calculated on
!     L_PINT    is the number of Lagrange interpolated reference
!               pressures for the CO2 k-coefficients.
!     L_NTREF   is the number of refernce temperatures for the
!               k-coefficients
!     L_TAUMAX  is the largest optical depth - larger ones are set
!               to this value.
!
!     MAXEXP    The largest value (TAU/MU) used in exp(-TAU/MU), where
!               TAU is the optical depth, and MU is the cosine of the
!               zenith angle.
!
!     L_REFH2O  The number of different water-mixing ratio values for
!               the k-coefficients that are now CO2+H2O. 
!
!     L_NREFI   The spectral interval number of the IR reference
!               wavelength (i.e. the 9 micron band) (8-12 microns)
!
!     L_NREFV   The spectral interval number of the visible reference
!               wavelength (i.e. the 0.67 micron band) 
!
!----------------------------------------------------------------------C

      integer, parameter :: L_NLAYRAD  = L_LAYERS+1
      integer, parameter :: L_NLEVRAD  = L_LAYERS+2
      
      integer, parameter :: L_NSPECTI =  96
      integer, parameter :: L_NSPECTV =  84
      integer, parameter :: L_NGAUSS  = 17

      integer, parameter :: L_NPREF   = 25
      integer, parameter :: L_NTREF   = 7
      integer, parameter :: L_TAUMAX  = 35

      real*8, parameter  :: MAXEXP    = 35.0D0

      integer, parameter :: L_PINT    = 51

      integer, parameter :: L_REFH2O  = 10

      integer, parameter :: L_NREFV   = 63 !original 6
      integer, parameter :: L_NREFI   = 72 !original 7

      character(len=300), allocatable :: filename_Qext,dirname_Qext
      character(len=:), allocatable :: fullpath_Qext
      integer :: arglen

      end module radinc_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            radcommon_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      module radcommon_h

      use grid_h, only: L_LEVELS
      use radinc_h
      implicit none

!----------------------------------------------------------------------C
!
!                             radcommon.h
!                         FORTRAN PARAMETERS
!                          GCM2.0  Feb 2003
!
!----------------------------------------------------------------------C
!
!     WNOI       - Array of wavenumbers at the spectral interval
!                  centers for the infrared.  Array is NSPECTI
!                  elements long.
!     DWNI       - Array of "delta wavenumber", i.e., the width,
!                  in wavenumbers (cm^-1) of each IR spectral
!                  interval.  NSPECTI elements long.
!     WAVEI      - Array (NSPECTI elements long) of the wavelenght
!                  (in microns) at the center of each IR spectral
!                  interval.
!     WNOV       - Array of wavenumbers at the spectral interval
!                  center for the VISIBLE.  Array is NSPECTV
!                  elements long.
!     DWNV       - Array of "delta wavenumber", i.e., the width,
!                  in wavenumbers (cm^-1) of each VISIBLE spectral
!                  interval.  NSPECTV elements long.
!     WAVEV      - Array (NSPECTV elements long) of the wavelenght
!                  (in microns) at the center of each VISIBLE spectral
!                  interval.
!     SOLARF     - Array (NSPECTV elements) of solar flux (W/M^2) in
!                  each spectral interval.  Values are for 1 AU, and
!                  are scaled to the Mars distance elsewhere.
!     TAURAY     - Array (NSPECTV elements) of the pressure-independent
!                  part of Rayleigh scattering optical depth.
!     PTOP       - Pressure at the top of the radiation code coordinate;
!                  = smallest k-coefficient pressure (1.0E-6 mbar)
!     FZEROI     - Fraction of zeros in the IR CO2 k-coefficients, for
!                  each temperature, pressure, and spectral interval
!     FZEROV     - Fraction of zeros in the VISIBLE CO2 k-coefficients, for
!                  each temperature, pressure, and spectral interval
!
!     AEROSOL RADIATIVE OPTICAL CONSTANTS
!     Values are at the wavelenght interval center
!
!     MIE SCATTERING - Size distribution weighted
!     Qextv    - Extinction efficiency - in the visible.
!     QextREF  - Reference visible wavelength (.67 micron band)
!     Qscatv   - Scattering efficiency - in the visible.
!     WV       - Single scattering albedo - in the visible.
!     GV       - Asymmetry parameter - in the visible.
!
!     Qexti    - Extinction efficiency - in the infrared.
!     Qscati   - Scattering efficiency - in the infrared.
!     WI       - Single scattering albedo - in the infrared.
!     GI       - Asymmetry parameter - in the infrared.
!     
!  Water ice clouds
!
!     MIE SCATTERING - Size distribution weighted
!     Qextvc   - Extinction efficiency - in the visible.
!     QextREFc - Reference visible wavelength (.67 micron band)
!     Qscatvc  - Scattering efficiency - in the visible.
!     WVc      - Single scattering albedo - in the visible.
!     GVc      - Asymmetry parameter - in the visible.
!
!     Qextic   - Extinction efficiency - in the infrared.
!     Qscatic  - Scattering efficiency - in the infrared.
!     WIc      - Single scattering albedo - in the infrared.
!     GIc      - Asymmetry parameter - in the infrared.
!     

      REAL*8 WNOI(L_NSPECTI), DWNI(L_NSPECTI), WAVEI(L_NSPECTI)
      REAL*8 WNOV(L_NSPECTV), DWNV(L_NSPECTV), WAVEV(L_NSPECTV)
      REAL*8 SOLARF(L_NSPECTV), TAURAY(L_NSPECTV)

!      real*8 CO2I(L_NTREF,L_PINT,L_REFH2O,L_NSPECTI,L_NGAUSS)
!      real*8 CO2V(L_NTREF,L_PINT,L_REFH2O,L_NSPECTV,L_NGAUSS)
      real*8 CO2I(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTI,L_NGAUSS)
      real*8 CO2V(L_NTREF,L_NPREF,L_REFH2O,L_NSPECTV,L_NGAUSS)
      real*8 FZEROI(L_NSPECTI)
      real*8 FZEROV(L_NSPECTV)
      real*8 PGASREF(L_NPREF), TGASREF(L_NTREF)

      real*8 qextv(L_NSPECTV), qscatv(L_NSPECTV), wv(L_NSPECTV)
      real*8 gv(L_NSPECTV)
      real*8 QextREF

      real*8 qexti(L_NSPECTI), qscati(L_NSPECTI), wi(L_NSPECTI)
      real*8 gi(L_NSPECTI)

      real*8 qextvc(L_NSPECTV), qscatvc(L_NSPECTV), wvc(L_NSPECTV)
      real*8 gvc(L_NSPECTV)
      real*8 QextREFc

      real*8 qextic(L_NSPECTI), qscatic(L_NSPECTI), wic(L_NSPECTI)
      real*8 gic(L_NSPECTI)

      real*8 planckir(L_NSPECTI,8501)

      real*8 TAUREF(L_LEVELS+1)
      real*8 PFGASREF(L_PINT)

!     Hemispheric mean in the IR
      real*8, parameter :: UBARI = 0.5D0

!===Alex commented CIA===
!       real*8  :: kgbar_tab(L_NSPECTI,5)
!       real*8  :: kbbar_tab(L_NSPECTI,5)
!       real*8  :: khbar_tab(L_NSPECTI,5)

!     These are for the Gauss-split 0.95 case

      real*8, parameter :: GWEIGHT(L_NGAUSS) = [                       &
                      4.8083554740D-02, 1.0563099137D-01,              &
                      1.4901065679D-01, 1.7227479710D-01,              &
                      1.7227479710D-01, 1.4901065679D-01,              &
                      1.0563099137D-01, 4.8083554740D-02,              &
                      2.5307134073D-03, 5.5595258613D-03,              &
                      7.8426661469D-03, 9.0670945845D-03,              &
                      9.0670945845D-03, 7.8426661469D-03,              &
                      5.5595258613D-03, 2.5307134073D-03,  0.0D0 ]  


!      real*8, parameter :: GWEIGHT(33) = [                             &
!                      1.2897418221D-02, 2.9570423871D-02,              &
!                      4.5200293049D-02, 5.9198761346D-02,              & 
!                      7.1058094688D-02, 8.0349346713D-02,              &
!                      8.6736622146D-02, 8.9989039966D-02,              &
!                      8.9989039966D-02, 8.6736622146D-02,              &
!                      8.0349346713D-02, 7.1058094688D-02,              &
!                      5.9198761346D-02, 4.5200293049D-02,              &
!                      2.9570423871D-02, 1.2897418221D-02,              &
!                      6.7881148529D-04, 1.5563380985D-03,              &
!                      2.3789627921D-03, 3.1157242814D-03,              &
!                      3.7398997204D-03, 4.2289129849D-03,              &
!                      4.5650853761D-03, 4.7362652614D-03,              & 
!                      4.7362652614D-03, 4.5650853761D-03,              & 
!                      4.2289129849D-03, 3.7398997204D-03,              &
!                      3.1157242814D-03, 2.3789627921D-03,              &
!                      1.5563380985D-03, 6.7881148529D-04,  0.0D0 ]  




!     These are for the CO2+H2O k-coefficients

      real*8, parameter :: WREFCO2(L_REFH2O) = [                       &
                     9.999999D-1, 9.99999D-1, 9.9999D-1, 9.999D-1,     &
                     9.99D-1, 9.9D-1, 9.0D-1, 8.0D-1, 7.0D-1, 6.0D-1 ]

      real*8, parameter :: WREFH2O(L_REFH2O) = [                       &
                     1.0D-7, 1.0D-6, 1.0D-5, 1.0D-4, 1.0D-3, 1.0D-2,   &
                     1.0D-1, 2.0D-1, 3.0D-1, 4.0D-1                  ]

!     If the CO2 optical depth (top to the surface) is less than
!     this value, we place that Gauss-point into the "zeros"
!     channel.  NRC parameter.

!  TLIMITS - TLIMIT for solar part of the spectrum
!  TLIMITI - TLIMIT for the IR

      real*8, parameter :: TLIMITS = 1.0D-32
      real*8, parameter :: TLIMITI = 5.0D-32
      


      end module radcommon_h
      module fccsave_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                              fccsave_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  Cloud Radiative properties

      use grid_h,   only: L_LEVELS, L_LAYERS
      use radinc_h, only: L_NSPECTV, L_NSPECTI, L_NGAUSS
      implicit none

      real*8  TAUREFCLD(L_LEVELS+1)
      real*8  QEXTREFCLD(L_LEVELS+1)

      real*8  QXVCLD(L_LEVELS+1,L_NSPECTV)
      real*8  QSVCLD(L_LEVELS+1,L_NSPECTV)
      real*8  GVCLD(L_LEVELS+1,L_NSPECTV)

      real*8  QXICLD(L_LEVELS+1,L_NSPECTI)
      real*8  QSICLD(L_LEVELS+1,L_NSPECTI)
      real*8  GICLD(L_LEVELS+1,L_NSPECTI)

!  Dust Radiative properties

      real*8  QEXTREFDST(L_LEVELS+1)

      real*8  QXVDST(L_LEVELS+1,L_NSPECTV)
      real*8  QSVDST(L_LEVELS+1,L_NSPECTV)
      real*8  GVDST(L_LEVELS+1,L_NSPECTV)

      real*8  QXIDST(L_LEVELS+1,L_NSPECTI)
      real*8  QSIDST(L_LEVELS+1,L_NSPECTI)
      real*8  GIDST(L_LEVELS+1,L_NSPECTI)

      end module fccsave_h
      module standard_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            standard_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use grid_h
      use defines_h

      implicit none

      real*8  :: TAUCUM(L_LEVELS)
      real*8  :: SIGMA(L_LEVELS)
      real*8  :: PRDST(L_NPDST), TAUDST(L_NPDST)

!     Spatially varying dust - TAUSURF(J,I) is the (visible 0.67 micron)
!     dust opacity at the surface for each grid point
!     TAUTOT added 6/28/01  GCM1.7
!     TAUTOT is the dust opacity at the RPTAU reference pressure
!     for each grid point.

      REAL*8  :: TAUTOT, RPTAU

!     variables that were in bkdata.f
!
!   ALICEN  :  SURFACE ALBEDO FOR CO2 ICE ON GROUND - Northern cap
!   ALICES  :  SURFACE ALBEDO FOR CO2 ICE ON GROUND - Southern cap
!
!   CP      :  HEAT CAPACITY (OR SPECIFIC HEAT) OF CO2 GAS.
!                  ( UNITS OF JOULES PER ( KG * DEGREES KELVIN))
!
!   DECMAX  :  MAXIMUM SOLAR DECLINATION.
!   DSIG    :  COORDINATE FOR ATM. LAYERS.
!   DTM     :  UNROUNDED DEFAULT NUMBER OF MINUTES PER TIME STEP.
!
!   ECCN    :  ORBITAL ECCENTRICITY FOR MARS.
!   EGOGND  :  EMISSIVITY OF BARE GROUND OUTSIDE 15 MICRON BAND WIDTH.
!   EG15GND :  EMISSIVITY OF BARE GROUND INSIDE 15 MICRON BAND WIDTH.
!   EG15CO2 :  EMISSIVITY OF GROUND WITH CO2 ICE INSIDE 15 MICRON BAND
!                   WIDTH.
!   EPS     :  SHEAR < EPS  MEANS  NEAR ZERO SHEAR.
!
!   GRAV    :  ACCELERATION DUE TO GRAVITY (MKS)
!
!   IM      :  DEFAULT LONGITUDINAL DIMENSION OF GRID
!   ISIZE   :  X DIMENSION OF GRID ARRAYS
!
!   JM      :  DEFAULT LATITUDINAL DIMENSION OF GRID
!   JSIZE   :  Y DIMENSION OF GRID ARRAYS
!
!   KAPA    :  A THERMODYNAMIC CONSTANT.
!
!   LAYERS  :  Z DIMENSION OF GRID ARRAYS
!
!   NCYCLE  :  DEFAULT NUMBER OF TIME STEPS PER TIME-STEPPING CYCLE.
!                   (USED IN STEP.)
!   NC3     :  THE FULL COMP3 ROUTINE IS EXECUTED ONLY ONCE PER NC3
!                   TIME STEPS.
!   NLAY    :  DEFAULT NUMBER OF ATMOSPHERIC LAYERS
!   NSIZE2D :  SIZE OF 2-D GRID ARRAYS =ISIZE*JSIZE.
!   NTAPE   :  Number of tape volume (ex. fort.11_032   NTAPE = 32)
!              GCM1.7 update
!
!   PI      :  THE MATHEMATICAL CONSTANT PI.
!   PTROP   :  DEFAULT PRESSURE AT THE TROPOPAUSE.
!
!   RAD     :  RADIUS OF MARS IN KM. CONVERTED TO METERS IN
!                   SUBROUTINE INPUT.
!   RGAS    :  GAS CONSTANT FOR MARS.
!   ROT     :  ANGULAR (ROTATIONAL) VELOCITY OF MARS IN RADIANS
!                   PER SECOND.
!   ROTPER  :  ROTATIONAL PERIOD IN MARTIAN HOURS.
!   RPTAU   :  Reference Pressure optical depth;  6.1 mbar for now
!
!   SCALEP  :  MULTIPLY BY 100 TO CONVERT PRESSURE FROM MILLIBARS
!                   TO PASCALS.
!   STBO    :  STEFAN-BOLTZMAN CONSTANT.
!
!   TAUCRT  :  ICE-CLOUD OPTICAL DEPTH
!   TAUH    :  DEFAULT NUMBER OF HOURS BETWEEN WRITING OUTPUT TO HISTORY
!                     TAPE.
!   TAUID   :  DEFAULT STARTING DAY FOR A GCM RUN.
!   TAUIH   :  DEFAULT STARTING HOUR FOR A GCM RUN.
!   TAUTOT  :  TOTAL GLOBAL OPTICAL DEPTH
!   TDINCR  :  MINIMUM AND INCREMENT VALUES
!   TDMIN   :  FOR SETS OF DUST OPT. DEPTHS.
!   THD     :  THERMAL DIFFUSIVITY
!   TIINCR  :  MINIMUM AND INCREMENT VALUES
!   TIMIN   :  FOR SET OF ICE OPT. DEPTHS.
!
!   end original bkdata.f comments

!     THE C ARRAY IS A SMALLER RECORD-TYPE STRUCTURE USED TO GROUP
!     TOGETHER A NUMBER OF GLOBALLY USED CONSTANTS.

      real*8  :: RSDIST

!      real*8, parameter :: dsig(L_LAYERS) = [                          &
!                          3.2718D-07, 5.3944D-07, 8.8940D-07,          &
!                          1.4663D-06, 2.4175D-06, 3.9859D-06,          &
!                          6.5718D-06, 1.0835D-05, 1.7863D-05,          &
!                          2.9452D-05, 4.8556D-05, 8.0060D-05,          &
!                          1.3199D-04, 2.1763D-04, 3.5879D-04,          &
!                          5.9155D-04, 9.7530D-04, 1.6080D-03,          &
!                          2.6512D-03, 4.3714D-03, 7.2060D-03,          &
!                          1.1882D-02, 1.9590D-02, 3.2298D-02,          &
!                          5.3245D-02, 8.7800D-02, 1.4475D-01,          &
!                          2.3865D-01, 3.9347017648D-01 ]

      real*8, parameter :: DSIG(L_LAYERS) = [                          &
                  0.0001237D0, 0.0002186D0, 0.0003877D0, 0.0006877D0,  &
                  0.0012199D0, 0.0021622D0, 0.0038572D0, 0.0057287D0,  &
                  0.0094287D0, 0.0155443D0, 0.0256262D0, 0.0422609D0,  &
                  0.0696418D0, 0.1148966D0, 0.1530633D0, 0.1572269D0,  &
                  0.1352191D0, 0.1212412D0, 0.0718975D0, 0.0435678D0,  &
                  0.014D0,     0.008D0,     0.003D0,     0.001D0      ]

      real*8  :: MWCO2 = 4.41D-2
      real*8  :: MWH2O = 1.80153D-2
      real*8  :: MWRATIO

      integer :: NPDST     = L_NPDST
      integer :: LAYERS    = L_LAYERS

      integer :: NLAY    = L_LAYERS
      real*8  :: PSF
      real*8  :: PTROP

!  Standard value = 0.03
!     DATA  CONRNU / 0.03  /   ! Standard value  ~25km half-height
!     DATA  CONRNU / 0.003 /   ! ~50 km half-height
!     DATA  CONRNU / 0.5   /   ! ~10 km half-height

      real*8  :: CONRNU

      end module standard_h
      module comp3cmn_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            comp3cmn_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use grid_h, only: L_LAYERS, L_LEVELS
      implicit none

      real*8  :: PL(L_LEVELS)
      real*8  :: TL(L_LEVELS)

      end module comp3cmn_h
      module cldcommon_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!                            cldcommon_h
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use grid_h, only: L_LAYERS
      use radinc_h, only: L_NSPECTI, L_NSPECTV
      implicit none

      integer, parameter :: nlonv   = L_NSPECTV
      integer, parameter :: nloni   = L_NSPECTI

      end module cldcommon_h
