      subroutine boxinterp(ft1p1,ft2p1,ft1p2,ft2p2,tref1,tref2,pref1,  &
                 pref2,tmid,pmid,tinter,pinter,ans)


!     GCM v24   2012
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

!----------------------------------------------------------------------!
!
!     GCM2.2   $Revision: 6 $
!     Last change date:
!     $Date: 2013-03-22 14:43:12 -0700 (Fri, 22 Mar 2013) $
!
!----------------------------------------------------------------------!

!                   T2 .FT2P1                    .FT2P2
!                                                 
!                                                 
!                   T1 .FT1P1                    .FT1P2
!                      P1                        P2

!----------------------------------------------------------------------!

      implicit none

      real*8 :: ft1p1, ft2p1, ft1p2, ft2p2, tref1, tref2
      real*8 :: pref1, pref2, tmid, pmid, ans1, ans2, ans
      logical :: tinter, pinter

!======================================================================C

      if(.not.tinter .and. .not.pinter) then
        ans = ft1p1
      elseif(.not.tinter .and. pinter) then
        ans = ft1p1 + (ft1p2 - ft1p1)*(pmid - pref1)/(pref2 - pref1)
      elseif(tinter .and. .not.pinter) then
        ans = ft1p1 + (ft2p1 - ft1p1)*(tmid - tref1)/(tref2 - tref1)
      elseif(tinter .and. pinter) then
        ans1 = ft1p1 + (ft2p1 - ft1p1)*(tmid - tref1)/(tref2 - tref1)
        ans2 = ft1p2 + (ft2p2 - ft1p2)*(tmid - tref1)/(tref2 - tref1)
        ans  = ans1 + (ans2 - ans1)*(pmid - pref1)/(pref2 - pref1)
      endif 
  
      return
      end
