subroutine importanceSampling

implicit none

real a1,a2,b1,b2,mx1,mx2,sx1,sx2,x1,x2,x01,x02,q1,q2
integer size1,i,count, t1,t2
real u
!real x1,x2
real inv_cdf, pdf
real failureFunction
real charts
real init1x1, init1x2, init2x1, init2x2, temp1,temp2
real omega1, omega2
integer samplingPoints, samplingPointsx1, samplingPointsx2
real inv_cdfn
integer j,k,l
integer maxIter, iter
real x1x1, x1x2, x2x1, x2x2, x1x1_temp,x1x2_temp,x2x1_temp,x2x2_temp
real pdf_cur, pdf_max
real curx1, curx2
real u1

real, allocatable:: sampledPoints(:,:)

samplingPoints=10
maxIter=1

allocate(sampledPoints(samplingPoints,2))

!default parameters of x1
t1=5 !exponential
a1= 1.0
b1= 10.0
mx1= 0.4
sx1= 0.2
x1= 0.1
x01= 0.3
q1= 0.3

!default parameters of x2
t2=1
a2= 2.0
b2= 2.5
mx2= 1.0
sx2= 0.4
x2= 2.3
x02= 2.3
q2= 9.3


!actual parameters
t1=6 !ex-max type1
mx1=50.0
sx1=5.0

t2=4
mx2=28.8*10**4
sx2=2.64*10**4
x02=19.9*10**4


!initial points
init1x1=70
init1x2=3278*init1x1
init2x1=init1x1+0.5*sx1
init2x2=3278*init2x1

x1x1=init1x1  !1st point horizontal
x1x2=init1x2  !1st point vertical
x2x1=init2x1  !2nd point horizontal
x2x2=init2x2  !2nd point vertical

!clearing the data files
open(10,file='data/sampledPoints.txt',access='sequential')
write(10,*)
close(10)

open(10,file='data/rootPoints.txt',access='sequential')
write(10,*)
close(10)


do iter=1,maxIter

open(10,file='data/rootPoints.txt',access='append')
write(10,*) x1x1, x1x2
write(10,*) x2x1, x2x2
close(10)

if((pdf(a1,b1,mx1,sx1,init1x1,x01,t1)*pdf(a2,b2,mx2,sx2,init1x2,x02,t2)).lt.&
   (pdf(a1,b1,mx1,sx1,init2x1,x01,t1)*pdf(a2,b2,mx2,sx2,init2x2,x02,t2))) then
   write(*,*) "Changing"
   temp1=init1x1
   temp2=init1x2
end if

write(*,*) failureFunction(init1x1,init1x2), failureFunction(init2x1,init2x2)

omega1=(pdf(a1,b1,mx1,sx1,x1x1,x01,t1)*pdf(a2,b2,mx2,sx2,x1x2,x02,t2))/&
        (pdf(a1,b1,mx1,sx1,x1x1,x01,t1)*pdf(a2,b2,mx2,sx2,x1x2,x02,t2)+&
        pdf(a1,b1,mx1,sx1,x2x1,x01,t1)*pdf(a2,b2,mx2,sx2,x2x2,x02,t2))

omega2=(pdf(a1,b1,mx1,sx1,x2x1,x01,t1)*pdf(a2,b2,mx2,sx2,x2x2,x02,t2))/&
        (pdf(a1,b1,mx1,sx1,x1x1,x01,t1)*pdf(a2,b2,mx2,sx2,x1x2,x02,t2)+&
        pdf(a1,b1,mx1,sx1,x2x1,x01,t1)*pdf(a2,b2,mx2,sx2,x2x2,x02,t2))

write(*,*) "omega1:", omega1
write(*,*) "omega2:", omega2

!problem here: if the values of pdfs are too small, the omega isn't calculated
!therefore the choice of the starting point is important

write(*,*) "omega1+2", omega1+omega2
samplingPointsx1=nint(samplingPoints*omega1)  !rounds to the nearest whole number
samplingPointsx2=nint(samplingPoints*omega2)

write(*,*) "Sampling points x1: ", samplingPointsx1
write(*,*) "Sampling points x2: ", samplingPointsx2


do i=1,samplingPointsx1
    call random_number(u)
    sampledPoints(i,1) = inv_cdfn(x1x1,sx1,u)

    call random_number(u)
    sampledPoints(i,2) = inv_cdfn(x1x2,sx2,u)
end do

do i=1,samplingPointsx2
    call random_number(u)
    sampledPoints(i+samplingPointsx1,1) = inv_cdfn(x2x1,sx1,u)

    call random_number(u)
    sampledPoints(i+samplingPointsx1,2) = inv_cdfn(x2x2,sx2,u)
end do


!writing just sampled points in the file
open(11,file='data/sampledPoints.txt',access='append')

do j=1,samplingPoints
    write(11,*) sampledPoints(j,1), sampledPoints(j,2)
end do




!looking in the sampled points one, with maximal pdf value, which will be the initial point for the next iteration step
pdf_max=0
x1x1_temp=0.0
x1x2_temp=0.0
do k=1,samplingPoints
    curx1=sampledPoints(k,1)
    curx2=sampledPoints(k,2)
    if(failureFunction(curx1,curx2).le.0.0) then
        pdf_cur=pdf(a1,b1,mx1,sx1,curx1,x01,t1)*pdf(a2,b2,mx2,sx2,curx2,x02,t2)
        if(pdf_cur.gt.pdf_max) then
            pdf_max=pdf_cur
            x1x1_temp=curx1
            x1x2_temp=curx2
        end if
    end if
end do

if((x1x1_temp==0.0).and.(x1x2_temp==0.0)) then !if no suitable point was found
    call random_number(u1)
    do k=1,10000
        if(u1.lt.omega1) then !generation of new point accoring to point x1

            call random_number(u)
            curx1 = inv_cdfn(x1x1,sx1,u)

            call random_number(u)
            curx2 = inv_cdfn(x1x2,sx2,u)

        else !generation of new point accoring to point x2
            call random_number(u)
            curx1 = inv_cdfn(x2x1,sx1,u)

            call random_number(u)
            curx2 = inv_cdfn(x2x2,sx2,u)
        end if

        write(11,*) curx1,curx2  !writing the sampled points to the file
        !new point generated. checking if its okay with our prerequisites
        if(failureFunction(curx1,curx2).le.0.0) then
            x1x1_temp=curx1
            x1x2_temp=curx2
            exit
        end if
    end do
end if

x1x1=x1x1_temp
x1x2=x1x2_temp

!looking in the sampled points a point for x2
pdf_max=0
do k=1,samplingPoints
    curx1=sampledPoints(k,1)
    curx2=sampledPoints(k,2)
    if((curx1.ne.x1x1).and.(curx2.ne.x1x2)) then  !excluding the already chosen point
        if(failureFunction(curx1,curx2).le.0.0) then  !searching only in the failure domain
            if((abs(curx1-x1x1).gt.0.9*sx1).and.(abs(curx2-x1x2).gt.0.9*sx2)) then  !searching only for points outside the cuboid
                pdf_cur=pdf(a1,b1,mx1,sx1,curx1,x01,t1)*pdf(a2,b2,mx2,sx2,curx2,x02,t2)
                if(pdf_cur.gt.pdf_max) then
                    pdf_max=pdf_cur
                    x2x1=curx1
                    x2x2=curx2
                end if
            end if
        end if
    end if
end do

close(11) !samplingPoints
call plotSampling(iter)
end do
end subroutine

subroutine plotSampling(iter)
    implicit none
    integer iter

open(40,file='data/plotSampling.plt')
    write(40,*) 'set terminal pngcairo enhanced font "Verdana,10"'
    write(40,'(A,I6.6,A)') 'set output "data/sampling',iter,'.png"'
    write(40,*) 'set style line 1 lt 1 lc rgb "#0000FF" lw 1'
    write(40,*) 'set xrange [35:90]'
    write(40,*) 'set yrange [100000:500000]'
    write(40,*) 'g(x) = 3278.66*x'
    write(40,*) 'set xlabel "X1"'
    write(40,*) 'set ylabel "X2"'
    write(40,*) 'plot g(x) title "g(x1,x2)=0" ls 1,\'
    write(40,*) '"data/rootPoints.txt" u 1:2 title "Root points",\'
    write(40,*) '"data/sampledPoints.txt" u 1:2 title "Sampled points"'
close(40)

    call system('gnuplot data/plotSampling.plt')

end subroutine

real function h(omega1,omega2,x1x1,x1x2,x2x1,x2x2,mx1_x1,mx2_x1,mx1_x2,mx2_x2,sx1,sx2)
    real omega1,omega2,x1x1,x1x2,x2x1,x2x2,mx1_x1,mx2_x1,mx1_x2,mx2_x2,sx1,sx2

    h=omega1*pdfn(mx1_x1,sx1,x1x1)*pdfn(mx2_x1,sx2,x1x2)+&
      omega2*pdfn(mx1_x2,sx1,x2x1)*pdfn(mx2_x2,sx2,x2x2)
end function
