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

write(*,*) "Estiamting the failure probability using the Monte Carlo approach."
write(*,*) "Sampled number of points: ", samplingMonteCarlo

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
close(10)



do i=1,samplingMonteCarlo
    x1=allPoints(i,1)
    x2=allPoints(i,2)
	if(failureFunction(x1,x2).lt.0.0) then
		count=count+1
		write(21,*) x1,x2  !writing data to failure file
    else
		write(20,*) x1,x2   !writing data to survival file
	end if

	if(mod(i,pointsPerChart)==0.0) then
        call plotFailure(i,x1max,x2max,x1min,x2min)
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

end subroutine


real function failureFunction(x1,x2)
    use variables
	implicit none
	real x1,x2

failureFunction = c1*x1+c2*x2

end function

subroutine plotFailure(i)
    implicit none
    integer i
    real x1max,x1min,x2max,x2min,c1_r

open(10,file='data/chartRange.txt')
read(10,*) x1max
read(10,*) x1min
read(10,*) x2max
read(10,*) x2min
read(10,*) c1_r
close(10)

    open(10,file='data/plotPoints.plt')
    write(10,*) 'set terminal pngcairo enhanced font "Verdana,10"'
    write(10,*) 'set style line 1 lt 1 lc rgb "#0000FF" lw 1'
    write(10,'(A,I7.7,A)') 'set output "data/failure',i,'.png"'
    write(10,*) 'g(x) = ',c1_r,'*x'
    write(10,*) 'set xlabel "X1"'
    write(10,*) 'set ylabel "X2"'
    write(10,*) 'set xrange [',x1min-0.2*abs(x1min),':',x1max+0.2*abs(x1max),']'
    write(10,*) 'set yrange [',x2min-0.2*abs(x2min),':',x2max+0.2*abs(x2max),']'
    write(10,*) 'plot g(x) title "g(x1,x2)=0" ls 1,\'
    write(10,*) '"data/pointsSurvival.txt" u 1:2 title "Survival",\'
    write(10,*) '"data/pointsFailure.txt" u 1:2 title "Failure"'
    close(10)

    call system('gnuplot data/plotPoints.plt')
    write(*,*) "chart plotted for i=",i
end subroutine

subroutine plotConvergence(maxProbability)
    implicit none
    real maxProbability

    open(40,file='data/plotConvergence.plt')
    write(40,*) 'set terminal pngcairo enhanced font "Verdana,10"'
    write(40,*) 'set output "data/convergence.png"'
    write(40,*) 'set xlabel "Sampled points"'
    write(40,*) 'set ylabel "Pf"'
    write(40,*) 'set yrange [0:',maxProbability*1.1,']'
    write(40,*) 'plot "data/convergence.txt" u 1:2 title "Failure Probability" with lines'
    close(40)

    call system('gnuplot data/plotConvergence.plt')
    !write(*,*) "Convergence plotted"
end subroutine
