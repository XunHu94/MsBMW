 !***********************************************************************************
    subroutine Convolution(num1,xtmp,ytmp,ztmp,den,speed,rarry) 
    !c----------------------------------------------------------------------
    !c
    !c----------------------------------------------------------------------
    use varmod
    implicit none
    integer::num1
    integer::i,j,k,xtmp,ytmp,ztmp
    real::dt,f0,sinang
    real,dimension(0:nz-1)::reflect
    real,dimension(MAXZ)::den,speed,speeds,EI,kk
    real,dimension(0:nz+NumSampledWalve-1)::sny
    real,dimension(num1)::rarry

    sny(0:nz+NumSampledWalve-1)=0.0
    reflect(0)=0.0
    EI(1)=den(1)*speed(1)
    do k=1,nz-1
        EI(k+1)=den(k+1)*speed(k+1)
        reflect(k)=(EI(k+1)-EI(k))/(EI(k+1)+EI(k))
    end do
    
    if(NumSampledWalve.gt.nz)then
        do i=0,nz+NumSampledWalve-1
            if(i.lt.nz)then
                sny(i)=0.0
                do j=0,i
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            else if(i.ge.nz.and.i.lt.NumSampledWalve)then
                sny(i)=0.0
                do j=i-nz+1,i
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            else if(i.ge.NumSampledWalve)then
                sny(i)=0.0
                do j=i-nz+1,NumSampledWalve
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            end if
        end do
    else
        do i=0,nz+NumSampledWalve-1
            if(i.lt.NumSampledWalve)then
                sny(i)=0.0
                do j=0,i
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            else if(i.ge.NumSampledWalve.and.i.lt.nz-1)then
                sny(i)=0.0
                do j=0,NumSampledWalve
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            else if(i.ge.nz-1)then
                sny(i)=0.0
                do j=i-nz+1,NumSampledWalve 
                    sny(i)=sny(i)+wavelet(j)*reflect(i-j)
                end do
            end if
        end do
    end if

    do i=1,num1
        rarry(i)=sny(i+HalfNumSampledWalve+ztmp)
    end do
    sny(0:num1+NumSampledWalve-1)=0.0
    end subroutine Convolution
    !***********************************************************************************