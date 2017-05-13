subroutine sortSamples(samples,samplesSize)
    implicit none
    integer samplesSize
    real samples(samplesSize), samplesTemp(samplesSize)
    integer i,j
    integer samplesOrder(samplesSize)

    samplesTemp=samples
    samplesOrder=1

    do i=1,samplesSize
        do j=1,samplesSize
            if(j.ne.i) then
                if(samples(i).gt.samples(j)) then
                    samplesOrder(i) = samplesOrder(i) + 1
                end if
            end if
        end do
    end do


    do i=1,samplesSize
        samplesTemp(samplesOrder(i)) = samples(i)
    end do

    samples = samplesTemp


end subroutine
