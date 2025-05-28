      subroutine funcd(solar,downir,rhouch,rhoucht,scond,stemp,sthick,
     *                 tg,f,df)

C  March 2002   c-grid

      implicit none

C  TG is the "X" that we are solving for in the grand scheme of things.
C  It is, of course, the ground temperature variable.

      real*8 tg, f, df

      real*8 solar, scond, stemp, sthick, stbo
      real*8 downir, rhoucht, rhouch

C======================================================================C

      stbo   = 5.6696E-8

C     F is the function value, df the derivitave

c     print*,"funcd:" ,solar,downir
      f  = solar + downir + rhoucht - rhouch*TG +
     *                      2.0*scond*(stemp-TG)/sthick - stbo*tg**4
      df = -rhouch -2.0*scond/sthick - 4.0*stbo*TG**3

      return
      end
