subroutine generateCdfs(samplesSize,a,b,mx,sx,q,x0,t)
    implicit none
    integer samplesSize
    real samples(samplesSize), samplesInvCDF(samplesSize)
    real theoreticalCDF(samplesSize), theoreticalPDF(samplesSize)
    integer i
    real u
    real x0,q
    real a,b,mx,sx, inv_cdf,cdf, pdf
    integer t


    do i=1,samplesSize
        call random_number(u)
        samples(i) = u
    end do

    call sortSamples(samples,samplesSize)

    do i=1,samplesSize
        samplesInvCDF(i) = inv_cdf(a,b,mx,sx,samples(i),x0,t)
    end do

    do i=1,samplesSize
        theoreticalCDF(i) = cdf(a,b,mx,sx,samplesInvCDF(i),t)
    end do

    open(10,file='empiricCDF.txt')
    do i=1,samplesSize
        write(10,*) samplesInvCDF(i), samples(i)
    end do
    close(10)

    open(10,file='theoryCDF.txt')
    do i=1,samplesSize
        write(10,*) samplesInvCDF(i), theoreticalCDF(i)
    end do
    close(10)

    open(10,file='theoreticalPDF.txt')
    do i=1,samplesSize
        write(10,*) samplesInvCDF(i), pdf(a,b,mx,sx,samplesInvCDF(i),x0,t)
    end do

    call system('gnuplot\wgnuplot plot.plt')
    call system('gnuplot\wgnuplot plot2.plt')
end subroutine
