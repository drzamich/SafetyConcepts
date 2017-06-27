subroutine failureMonteCarlo

use variables
implicit none

real x1,x2
integer size1,i,count
real u
real inv_cdf, pdf
real failureFunction
real charts
integer pointsPerChart
real, allocatable:: allPoints(:,:)
real x1max, x2max, x1min,x2min
real probability, maxProbability
integer timeStart, timeEnd, k

write(*,*) "Estiamting the failure probability using the Monte Carlo approach."
write(*,*) "Sampled number of points: ", samplingMonteCarlo
call system_clock(count=timeStart)
count=0 !number of points in the failure region
maxProbability=0.0 !value of max. probability for the convergence chart
allocate(allPoints(samplingMonteCarlo,2))

pointsPerChart=samplingMonteCarlo/chartsMonteCarlo

!clearing the data files

        open(20,file='data/pointsSurvival.txt',access='sequential')
        write(20,*) "x1", "x2"
        close(20)

        open(21,file='data/pointsFailure.txt',access='sequential')
        write(21,*) "x1", "x2"
        close(21)

        open(21,file='data/exVal.txt',access='sequential')
        write(21,*) mx1, mx2
        close(21)

        open(30,file='data/convergence.txt',access='sequential')
        write(30,*) "points", "failure probability"


        open(20,file='data/pointsSurvival.txt',access='append')
		open(21,file='data/pointsFailure.txt',access='append')


do i=1,samplingMonteCarlo
	call random_number(u)
	x1= inv_cdf(a1,b1,mx1,sx1,u,x01,t1)

	call random_number(u)
	x2= inv_cdf(a2,b2,mx2,sx2,u,x02,t2)

	if(i==1) then
    x1max=x1
    x2max=x2
    x1min=x1
    x2min=x2
    endif

    allPoints(i,1) = x1
    allPoints(i,2) = x2

    if(x1.gt.x1max) x1max=x1
    if(x2.gt.x2max) x2max=x2
    if(x1.lt.x1min) x1min=x1
    if(x2.lt.x2min) x2min=x2
enddo

open(10,file='data/chartRange.txt')
write(10,*) x1max
write(10,*) x1min
write(10,*) x2max
write(10,*) x2min
write(10,*) c1_r
write(10,*) c0
close(10)


k=1
do i=1,samplingMonteCarlo
    x1=allPoints(i,1)
    x2=allPoints(i,2)
	if(failureFunction(x1,x2).lt.0.0) then
		count=count+1
		write(21,*) x1,x2  !writing data to failure file
    else
        if((int(k/montecarloPointsRatio)).eq.i) then
            write(20,*) x1,x2   !writing data to survival file
            k=k+1
		end if
	end if

	if(mod(i,pointsPerChart)==0.0) then
        call plotMonteCarlo(i,x1max,x2max,x1min,x2min)
        write(*,*) i
    end if

    probability = real(count)/real(i)
    write(30,*) i, probability   !writing data to convergence file
    if(probability.gt.maxProbability) maxProbability=probability
end do

		close(20)
		close(21)
        close(30)
        close(31)

call plotConvergence(maxProbability)

write(*,*) "Failure probability:", real(count)/real(samplingMonteCarlo)

call system_clock(count=timeEnd)

write(*,*) "Time needed [s]: ", (timeEnd-timeStart)/1000
end subroutine
