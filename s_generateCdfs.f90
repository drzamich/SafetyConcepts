subroutine generateCdfs(samplesSize,a,b,mx,sx,q,x0,t,barSize)
    implicit none
    integer samplesSize
    real samples(samplesSize), samplesInvCDF(samplesSize)
    real theoreticalCDF(samplesSize), theoreticalPDF(samplesSize)
    integer i, j
    real u
    real x0,q
    real a,b,mx,sx, inv_cdf,cdf, pdf
    integer t
    real bottomRange, topRange, rangeWidth, mediumRange, bottomRangeTemp
    integer histogramBarsInt
    real histogramBarsReal
    real barSize
    real,allocatable:: histogramData(:,:)
    real valueCounter


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
        write(10,*) samplesInvCDF(i), real(i)/samplesSize
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
    close(10)

    !
    ! GENERATING DATA FOR HISTOGRAMS
    !
    bottomRange = samplesInvCDF(1)
    topRange = samplesInvCDF(samplesSize)
    rangeWidth = topRange-bottomRange
    histogramBarsReal = 1/barSize
    histogramBarsInt=int(histogramBarsReal)
    write(*,*) histogramBarsInt


    allocate(histogramData(histogramBarsInt,2))

    bottomRangeTemp=bottomRange
    do i=1,histogramBarsInt
        valueCounter=0
        topRange = bottomRange+barSize*rangeWidth
        mediumRange = (topRange+bottomRange)*0.5
        histogramData(i,1) = mediumRange
        write(*,*) "step",i
        write(*,*) "bottom range", bottomRange
        write(*,*) "top range", topRange
        do j=1,samplesSize
            if(j.ne.samplesSize) then
                if((samplesInvCDF(j).ge.bottomRange).and.(samplesInvCDF(j).lt.topRange)) then
                    valueCounter = valueCounter+1
                end if
            else
                if((samplesInvCDF(j).ge.bottomRange).and.(samplesInvCDF(j).le.topRange)) then
                    valueCounter = valueCounter+1
                end if
            end if
        end do
        write(*,*) "items in range:", valueCounter
        histogramData(i,2)=valueCounter/samplesSize
        bottomRangeTemp=bottomRange
        bottomRange=topRange
    end do

    open(10,file='histograms.txt')
    do i=1,histogramBarsInt
        write(10,*) histogramData(i,1), histogramData(i,2)
    end do
    close(10)
    call system('gnuplot\wgnuplot plot.plt')
    call system('gnuplot\wgnuplot plot2.plt')
end subroutine
