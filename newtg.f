      subroutine newtg(solar,downir,rhouch,ta,scond,stemp,
     *                 sthick,tg)

C  March 2002   c-grid

C  Use modified Newton-Raphson method for finding the ground
C  temperature.

      implicit none

      integer J, N
      real*8  df, dx, dxold, f, fh, fl, temp, TH, TL
      real*8  T1, T2, xacc
      real*8  solar, scond, stemp, sthick, tg, ta
      real*8  downir, rhouch, rhoucht

C  MAXIT is the maximum number of iterations allowed for convergence.

      integer MAXIT
      parameter(MAXIT = 30)

C  T1 is the initial low temperature
C  T2 is the initial high temperature
C  XACC is the accuracy (in Kelvins)

      parameter(T1   = 50.0)
      parameter(T2   = 350.0)
      parameter(XACC = 0.001)

C======================================================================C

c      print*,' IN NEWTG: Solar=',solar
      rhoucht = rhouch * ta

      call funcd(solar,downir,rhouch,rhoucht,scond,stemp,sthick,T1,FL,
     *           df)
      call funcd(solar,downir,rhouch,rhoucht,scond,stemp,sthick,T2,FH,
     *           df)

      if(FL.eq.0.0) then
        TG = T1
        goto 100
      elseif(FH.eq.0.0) then
        TG = T2
        goto 100
      elseif(FL.LT.0.0) then
        TL = T1
        TH = T2
      else
        TL = T2
        TH = T1
      end if

      TG    = 0.5*(T1+T2)
      dxold =  abs(T2-T1)
      dx    = dxold

      call funcd(solar,downir,rhouch,rhoucht,scond,stemp,sthick,Tg,f,df)
      
      do J=1,MAXIT
c      print*,j,tg, f

        if(((TG-TH)*df-f)*((TG-TL)*df-f) .ge. 0.0   .or.
     *       abs(2.0*f) .gt. abs(dxold*df)               ) then

          dxold = dx
          dx    =  0.5*(TH-TL)
          TG    = TL + dx
          if(TL.eq.TG) then
            goto 100
          end if
        else
          
          dxold = dx
          dx    = f/df
          temp  = TG
          TG    = TG - dx
          if(TEMP.eq.TG) then
            goto 100
          end if
        end if

        if(abs(dx).lt.xacc) then
           goto 100
        end if

        call funcd(solar,downir,rhouch,rhoucht,scond,stemp,sthick,TG,f,
     *             df)

        if(F.lt.0.0) then
          TL = TG
        else
          TH = TG
        end if

      end do

C  If we reach this statement, we've done MAXIT number of iterations
C  without convergence.

      write(6,'("------- SUBROUTINE NEWTG")')
      write(6,'("Maximum number of iterations ",i3," exceeded.")')
      stop

  100 continue
      return
      end
