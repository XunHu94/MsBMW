 !***********************************************************************************
    SUBROUTINE upscale(ix,iy,iz)
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    integer :: upix,upiy,upiz,Fineix,Fineiy,Fineiz
    INTEGER, INTENT(IN) :: ix,iy,iz

    upix = (ix-1)/ncoarse+1
    upiy = (iy-1)/ncoarse+1
    upiz = (iz-1)/ncoarse+1
    
    do Fineix = max((upix-1)*ncoarse+1,1),min(upix*ncoarse,nx)
        do Fineiy = max((upiy-1)*ncoarse+1,1),min(upiy*ncoarse,ny)
            do Fineiz = max((upiz-1)*ncoarse+1,1),min(upiz*ncoarse,nz)
                
                if(simim(Fineix,Fineiy,Fineiz)>0) then
                    simim(ix,iy,iz) = simim(Fineix,Fineiy,Fineiz)
                    DenSim(ix,iy,iz) = DenSim(Fineix,Fineiy,Fineiz)
                    VpSim(ix,iy,iz) = VpSim(Fineix,Fineiy,Fineiz)                     
                    SnySim(ix,iy,iz) = SnySim(Fineix,Fineiy,Fineiz)
                    go to 56
                end if
            end do
        end do
    end do 
    
56  continue
    
    END SUBROUTINE upscale
    
    !***********************************************************************************