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
    character(len=100) distributionName

    select case(t)
        case(1)      !uniform distribution
            distributionName = "Uniform distribution"
        case(2)      !standard normal distribution
            distributionName ="Standard normal distribution"
        case(3)      !normal distribution
            distributionName ="Normal distribution"
        case(4)      !log-normal distribution
            distributionName ="Log-Normal distribution"
        case(5)       !exponential distribution
            distributionName ="Exponential distribution"
    end select


    do i=1,samplesSize
        call random_number(u)
        samples(i) = u
    end do

    call sortSamples(samples,samplesSize)

    do i=1,samplesSize
        samplesInvCDF(i) = inv_cdf(a,b,mx,sx,samples(i),x0,t)
    end do

    do i=1,samplesSize
        theoreticalCDF(i) = cdf(a,b,mx,sx,samplesInvCDF(i),x0,t)
    end do

    open(10,file='data/empiricCDF.txt')
    do i=1,samplesSize
        write(10,*) samplesInvCDF(i), real(i)/samplesSize
    end do
    close(10)

    open(10,file='data/theoryCDF.txt')
    do i=1,samplesSize
        write(10,*) samplesInvCDF(i), theoreticalCDF(i)
    end do
    close(10)

    open(10,file='data/theoreticalPDF.txt')
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
    !write(*,*) histogramBarsInt


    allocate(histogramData(histogramBarsInt,2))

    open(10,file='data/histograms.txt')
    bottomRangeTemp=bottomRange
    do i=1,histogramBarsInt
        valueCounter=0
        topRange = bottomRange+barSize*rangeWidth
        mediumRange = (topRange+bottomRange)*0.5
        histogramData(i,1) = mediumRange
        !write(*,*) "step",i
        !write(*,*) "bottom range", bottomRange
        !write(*,*) "top range", topRange
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
        !write(*,*) "items in range:", valueCounter
        write(10,*) bottomRange, valueCounter/samplesSize*10
        write(10,*) topRange, valueCounter/samplesSize*10
        bottomRangeTemp=bottomRange
        bottomRange=topRange
    end do
    close(10)

    !open(10,file='data/histograms.txt')
    !do i=1,histogramBarsInt
     !   write(10,*) histogramData(i,1), histogramData(i,2)
    !end do
    !close(10)

    open(10,file='data/plot_cdf.plt')
    write(10,*) 'set title "Empiric vs theoretical CDF - ',trim(distributionName),'"'
    write(10,*) 'set xlabel "x"'
    write(10,*) 'set ylabel "F(x)"'
    write(10,*) 'set terminal png'
    write(10,*) 'set style data lines'
    write(10,'(A,I2.2,A)') 'set output "data/cdfs_',t,'.png"'
    write(10,*) 'plot "data/empiricCDF.txt" title "Empiric CDF", "data/theoryCDF.txt" title "Theoretical CDF"'
    close(10)

    open(10,file='data/plot_pdf.plt')
    write(10,*) 'set title "Empiric vs theoretical PDF - ',trim(distributionName),'"'
    write(10,*) 'set xlabel "x"'
    write(10,*) 'set ylabel "f(x)"'
    write(10,*) 'set terminal png'
    write(10,*) 'set style data lines'
    write(10,'(A,I2.2,A)') 'set output "data/pdfs_',t,'.png"'
    write(10,*) 'plot "data/histograms.txt" title "Histograms", "data/theoreticalPDF.txt" title "Theoretical PDF"'
    !write(10,*) 'plot "data/theoreticalPDF.txt" title "Theoretical PDF", "data/histograms.txt" title "Histograms" with histograms'
    close(10)

    call system('gnuplot\wgnuplot data\plot_cdf.plt')
    call system('gnuplot\wgnuplot data\plot_pdf.plt')
end subroutine
