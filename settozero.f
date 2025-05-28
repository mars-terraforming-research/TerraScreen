       subroutine settozero(nsize,array)

!     GCM v23   2010
!     Ames Mars GCM group
!     Jeffery Hollingsworth, PI
!     NASA Ames Research Center

       implicit none

       integer nsize
       real*8 array(nsize)

       integer i

       do i=1,nsize
         array(i) = 0.0D0
       enddo

       return

       end
       subroutine settozero4(nsize,array)

       implicit none

       integer nsize
       real*4 array(nsize)

       integer i

       do i=1,nsize
         array(i) = 0.0
       enddo

       return

       end
