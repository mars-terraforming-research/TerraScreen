      subroutine ORBIT(snom)
C
      real*8 TEMP(669),snom(669)

      data ec/0.093379/,pi/3.14159/
C
      SECC = SQRT(1. - EC**2)
      YEAR = 669.
      ACOS = PI/2. - ATAN(EC/SECC)
      AA1  = PI - ACOS
      AA2  = PI + ACOS
      ANOM = 0.0
      TDAY = 0.0
      SNOM(669) = 250.
      DO 8 I = 1,668
      DEL = 2.*PI/YEAR
      ANOM = ANOM + DEL
 5006 AA = SIN(ANOM)*SECC/(1. + EC*COS(ANOM))
      BB = ATAN(AA/SQRT(1.-AA**2))
      IF ( ANOM.GE.AA1.AND.ANOM.LT.AA2 ) BB = PI - BB
      IF ( ANOM.GE.AA2 ) BB = 2.0*PI + BB
      TNEW = (YEAR/2./PI)*(BB - EC*AA)
      TEST = TNEW - TDAY
      IF ( ABS(1.-TEST).LT.0.10 ) GO TO 5007
      IF ( 1.-TEST ) 5008,5007,5009
 5008 DEL = DEL/3.
      ANOM = ANOM - DEL
      GO TO 5006
 5009 ANOM = ANOM + DEL/2.
      GO TO 5006
 5007 SNOM(I) = ANOM
      SNOM(I) = 180.*SNOM(I)/PI
      SNOM(I) = SNOM(I) + 250.
      IF ( SNOM(I).GT.360. ) SNOM(I) = SNOM(I) - 360.
      TDAY = TDAY + 1.
    8 CONTINUE
C
C===> LET DAY 1 CORRESPOND TO LS=0
C
      DO 11 I=1,669
      J=I+185
      IF ( J.GT.669 ) J = J-669
      TEMP(I) = SNOM(J)
   11 CONTINUE
C
      DO 12 I = 1,669
      SNOM(I) = TEMP(I)
c     print*,i,snom(i)
   12 CONTINUE
C
      return
      END
