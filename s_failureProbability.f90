subroutine failureProbability
implicit none

real a1,a2,b1,b2,mx1,mx2,sx1,sx2,x1,x2,x01,x02,q1,q2
integer size1,i,count, t1,t2
real u
!real x1,x2
real inv_cdf, pdf
real failureFunction
real charts

count=0
size1=100000  !total number of sampled points
charts=5


! parameters of x1
!default parameters
t1=5 !exponential
a1= 1.0
b1= 10.0
mx1= 0.4
sx1= 0.2
x1= 0.1
x01= 0.3
q1= 0.3


!actual parameters
t1=6 !ex-max type1
mx1=50.0
sx1=5.0


!parameters of x2

!default parameters
t2=1
a2= 2.0
b2= 2.5
mx2= 1.0
sx2= 0.4
x2= 2.3
x02= 2.3
q2= 9.3


!actual parameters
t2=4
mx2=28.8*10**4
sx2=2.64*10**4
x02=19.9*10**4

charts = size1/charts
!clearing the data files

        open(20,file='data/pointsSurvival.txt',access='sequential')
        write(20,*) "x1", "x2"
        close(20)

        open(21,file='data/pointsFailure.txt',access='sequential')
        write(21,*) "x1", "x2"
        close(21)

        open(30,file='data/convergence.txt',access='sequential')
        write(30,*) "points", "failure probability"

        open(31,file='data/pdf.txt',access='sequential')
        write(30,*) "x1", "x2", "f(x1,x2)"

        open(20,file='data/pointsSurvival.txt',access='append')
		open(21,file='data/pointsFailure.txt',access='append')

do i=1,size1
	call random_number(u)
	x1= inv_cdf(a1,b1,mx1,sx1,u,x01,t1)

	call random_number(u)
	x2= inv_cdf(a2,b2,mx2,sx2,u,x02,t2)

	if(failureFunction(x1,x2).lt.0.0) then
		count=count+1
		write(21,*) x1,x2  !failure
		write(*,*) "failure prob2.: TEST pdf1 of failure:", pdf(a1,b1,mx1,sx1,72.138,x01,t1)
    else
		write(20,*) x1,x2   !survival
	end if

	if(mod(real(i),charts)==0.0) then
        call plotFailure(i)
    end if

    write(30,*) i, real(count)/real(i)
    write(31,*) x1, x2, pdf(a1,b1,mx1,sx1,x1,x01,t1)*pdf(a2,b2,mx2,sx2,x2,x02,t2)
end do

		close(20)
		close(21)
        close(30)
        close(31)

call plotConvergence

write(*,*) "Failure probability:", real(count)/real(size1)

end subroutine


real function failureFunction(x1,x2)
	implicit none
	real x1,x2
	real l1,l2,wpl,cpl

	l1=3.0
	l2=2.0
	wpl=3.66*10**(-4)

	!cpl=(-(l1*l2)/(l1+l2))*(1/wpl)
    cpl=-3278.688525

failureFunction = x2+cpl*x1

end function

subroutine plotFailure(i)
    implicit none
    integer i

    open(10,file='data/plotPoints.plt')
    write(10,*) 'set terminal pngcairo enhanced font "Verdana,10"'
    write(10,*) 'set style line 1 lt 1 lc rgb "#0000FF" lw 1'
    write(10,'(A,I6.6,A)') 'set output "data/failure',i,'.png"'
    write(10,*) 'g(x) = 3278.66*x'
    write(10,*) 'set xlabel "X1"'
    write(10,*) 'set ylabel "X2"'
    !write(10,*) 'set xrange [0:1.5]'
    !write(10,*) 'set yrange [0:420000]'
    write(10,*) 'plot g(x) title "g(x1,x2)=0" ls 1,\'
    write(10,*) '"data/pointsSurvival.txt" u 1:2 title "Survival",\'
    write(10,*) '"data/pointsFailure.txt" u 1:2 title "Failure"'
    close(10)

    call system('gnuplot data/plotPoints.plt')
    write(*,*) "chart plotted for i=",i
end subroutine

subroutine plotConvergence
    open(40,file='data/plotConvergence.plt')
    write(40,*) 'set terminal pngcairo enhanced font "Verdana,10"'
    write(40,*) 'set output "data/convergence.png"'
    write(40,*) 'set xlabel "Sampled points"'
    write(40,*) 'set ylabel "Pf"'
    write(40,*) 'set yrange [-0.001:0.001]'
    write(40,*) 'plot "data/convergence.txt" u 1:2 title "Failure Probability" with lines'
    close(40)

    call system('gnuplot data/plotConvergence.plt')
    write(*,*) "Convergence plotted"
end subroutine

subroutine plotPdf
    open(40,file='data/plotPdf.plt')
    write(40,*) 'set terminal pngcairo enhanced font "Verdana,10"'
    write(40,*) 'set output "data/pdf.png"'
    write(40,*) 'set xlabel "X1"'
    write(40,*) 'set ylabel "X2"'
    write(40,*) 'set pm3d map impl'
    write(40,*) 'set contour'
    write(40,*) 'set style increment user'
    write(40,*) 'do for [i=1:18] { set style line i lc rgb "black"}'
    write(40,*) 'set cntrparam levels incr -0.3,0.1,0.5'
    write(40,*) 'set autoscale fix'
    write(40,*) 'splot "test.txt" w pm3d notitle'

    write(40,*) 'plot "data/convergence.txt" u 1:2 title "Failure Probability" with lines'
    close(40)

    call system('gnuplot data/plotConvergence.plt')
    write(*,*) "Convergence plotted"
end subroutine
