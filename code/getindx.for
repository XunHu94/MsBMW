      subroutine getindx(nn,mid,size,loc,index,inflag)
c-----------------------------------------------------------------------
c
c     Gets the coordinate index location of a point within a grid
c     ***********************************************************
c
c
c nn       number of "nodes" or "cells" in this coordinate direction
c mid     origin at the center of the first cell
c size     size of the cells
c loc     location of the point being considered
c index   output index within [1,n]
c inflag  true if the location is actually in the grid (false otherwise
c         e.g., if the location is outside then index will be set to
c         nearest boundary
c
c
c
c-----------------------------------------------------------------------
	use varmod
	implicit none
      integer::nn,index
      real::mid,size,loc
      logical::inflag
c
c Compute the index of "loc":
c
      index = int( (loc-mid)/size + 1.5 )
      
c
c Check to see if in or out:
c
      if(index.lt.1) then
            index  = 1
            inflag = .false.
      else if(index.gt.nn) then
            index  = nn
            inflag = .false.
      else
            inflag = .true.
      end if
c
c Return to calling program:
c
      return
      end
