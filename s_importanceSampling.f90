subroutine importanceSampling

implicit none

real a1,a2,b1,b2,mx1,mx2,sx1,sx2,x1,x2,x01,x02,q1,q2
integer size1,i, t1,t2
real u
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
real part1, part2, part3
real h
real indicator, indicator_val, pdf_xi1, pdf_xi2
real haha

real, allocatable:: iterationSamples(:)  !array which keeps number of sampled points in each iteration
real, allocatable:: iterationProbability(:) !array which keeps the value of probability in each iteration
real, allocatable:: allSampledPoints(:,:,:)
real maxPoints !maximal number of points to be generated in each iteration
real, allocatable:: h_params(:,:)

real, allocatable:: sampledPoints(:,:)
integer totalSampled
integer sampledThisIter

samplingPoints=10  !points sampled in one iteration
maxIter=200  !maximal number of iterations
maxPoints=1000

allocate(iterationSamples(maxIter))
allocate(iterationProbability(maxIter))
allocate(allSampledPoints(maxIter,maxPoints,2))
allocate(h_params(maxIter,6))


iterationSamples=0.0
iterationProbability=0.0
allSampledPoints=0.0
totalSampled=0


!actual parameters
t1=6 !ex-max type1
mx1=50.0
sx1=5.0


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

t2=4 !log-normal distr
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
open(10,file='data/sampledPoints.txt',access='sequential')  !file for keeping points sampled in current iteration
write(10,*)
close(10)

open(10,file='data/rootPointsX1.txt',access='sequential')
write(10,*)
close(10)

open(10,file='data/rootPointsX2.txt',access='sequential')
write(10,*)
close(10)

!writing the initial sampled points in the file
open(10,file='data/rootPointsX1.txt',access='append')
write(10,*) x1x1, x1x2
close(10)

open(10,file='data/rootPointsX2.txt',access='append')
write(10,*) x2x1, x2x2
close(10)


call plotSampling(0)
!allocate(sampledPoints(samplingPoints,2))
!
!
!  START OF THE ITERATION
!
!

do iter=1,maxIter
    write(*,*) iter
sampledThisIter=0

if((pdf(a1,b1,mx1,sx1,init1x1,x01,t1)*pdf(a2,b2,mx2,sx2,init1x2,x02,t2)).lt.&
   (pdf(a1,b1,mx1,sx1,init2x1,x01,t1)*pdf(a2,b2,mx2,sx2,init2x2,x02,t2))) then
   write(*,*) "Changing"
   temp1=init1x1
   temp2=init1x2
end if

omega1=(pdf(a1,b1,mx1,sx1,x1x1,x01,t1)*pdf(a2,b2,mx2,sx2,x1x2,x02,t2))/&
        (pdf(a1,b1,mx1,sx1,x1x1,x01,t1)*pdf(a2,b2,mx2,sx2,x1x2,x02,t2)+&
        pdf(a1,b1,mx1,sx1,x2x1,x01,t1)*pdf(a2,b2,mx2,sx2,x2x2,x02,t2))

omega2=(pdf(a1,b1,mx1,sx1,x2x1,x01,t1)*pdf(a2,b2,mx2,sx2,x2x2,x02,t2))/&
        (pdf(a1,b1,mx1,sx1,x1x1,x01,t1)*pdf(a2,b2,mx2,sx2,x1x2,x02,t2)+&
        pdf(a1,b1,mx1,sx1,x2x1,x01,t1)*pdf(a2,b2,mx2,sx2,x2x2,x02,t2))

h_params(iter,1)=omega1
h_params(iter,2)=omega2
h_params(iter,3) = x1x1
h_params(iter,4) = x1x2
h_params(iter,5) = x2x1
h_params(iter,6) = x2x2



write(*,*) "omega1:", omega1
write(*,*) "omega2:", omega2

!problem here: if the values of pdfs are too small, the omega isn't calculated
!therefore the choice of the starting point is important

write(*,*) "omega1+2", omega1+omega2
samplingPointsx1=nint(samplingPoints*omega1)  !rounds to the nearest whole number
samplingPointsx2=nint(samplingPoints*omega2)

write(*,*) "Sampling points x1: ", samplingPointsx1
write(*,*) "Sampling points x2: ", samplingPointsx2


!generating sampling points

do i=1,samplingPointsx1
    call random_number(u)
    allSampledPoints(iter,i,1) = inv_cdfn(x1x1,sx1,u)

    call random_number(u)
    allSampledPoints(iter,i,2) = inv_cdfn(x1x2,sx2,u)
end do

do i=1,samplingPointsx2
    call random_number(u)
    allSampledPoints(iter,i+samplingPointsx1,1) = inv_cdfn(x2x1,sx1,u)

    call random_number(u)
    allSampledPoints(iter,i+samplingPointsx1,2) = inv_cdfn(x2x2,sx2,u)
end do

totalSampled = totalSampled+samplingPointsx1+samplingPointsx2
sampledThisIter=samplingPointsx1+samplingPointsx2
iterationSamples(iter) = sampledThisIter

if(iter.ne.1) then

part3=0.0

do u=1,iter
        part2=0
        do i=1,iterationSamples(u)
        part1=0.0
        do j=1,iter
                        call h_test(h_params(j,1),h_params(j,2),h_params(j,3),h_params(j,4),h_params(j,5),h_params(j,6),&
                        sx1,sx2,allSampledPoints(u,i,1),allSampledPoints(u,i,2),haha)
                        part1=part1+(iterationSamples(j)/totalSampled)*haha
                        !h(h_params(j,1),h_params(j,2),h_params(j,3),h_params(j,4),h_params(j,5),h_params(j,6),&
                        !sx1,sx2,allSampledPoints(u,i,1),allSampledPoints(u,i,2))



                        !write(*,*) "parameters of h function"
                        !write(*,*) "omega1", h_params(j,1)
                        !write(*,*) "omega2", h_params(j,2)
                        !write(*,*) "x1x1", h_params(j,3)
                        !write(*,*) "x1x2", h_params(j,4)
                        !write(*,*) "x2x1", h_params(j,5)
                        !write(*,*) "x2x2", h_params(j,6)
                        !write(*,*) "Checked point:"
                        !write(*,*) "x1:", allSampledPoints(u,i,1)
                        !write(*,*) "x2:", allSampledPoints(u,i,2)
                        !write(*,*) "h:", h(h_params(j,1),h_params(j,2),h_params(j,3),h_params(j,4),h_params(j,5),h_params(j,6),&
                        !sx1,sx2,allSampledPoints(u,i,1),allSampledPoints(u,i,2))

                        !write(*,*) "h:", haha

        end do
        !write(*,*) "part1: ", part1
            !write(*,*) "x1:", allSampledPoints(u,i,1)
            !write(*,*) "x2:", allSampledPoints(u,i,2)
            indicator_val=indicator(allSampledPoints(u,i,1),allSampledPoints(u,i,2))
            pdf_xi1=pdf(a1,b1,mx1,sx1,allSampledPoints(u,i,1),x01,t1)
            pdf_xi2=pdf(a2,b2,mx2,sx2,allSampledPoints(u,i,2),x02,t2)
            if(isnan(pdf_xi1)) then
                pdf_xi1=0.0
            end if
            if(isnan(pdf_xi2)) then
                pdf_xi2=0.0
            end if
            !write(*,*) "x1: ", allSampledPoints(u,i,1)
            !write(*,*) "x2: ", allSampledPoints(u,i,2)
            !write(*,*) "indicator:", indicator_val
            !write(*,*) "pdf1:", pdf_xi1
            !write(*,*) "pdf2:", pdf_xi2
        part2=part2+((indicator_val*pdf_xi1*pdf_xi2)/part1)
        !write(*,*) "part2: ", part2

    end do
    part3=part3+part2
end do

iterationProbability(iter) = part3/totalSampled

write(*,*) "probability: ", iterationProbability(iter)
end if


!looking in the sampled points one, with maximal pdf value, which will be the initial point for the next iteration step
pdf_max=0.0
curx1=0.0
curx2=0.0
do k=1,sampledThisIter
    curx1=allSampledPoints(iter,k,1)
    curx2=allSampledPoints(iter,k,2)
    if(failureFunction(curx1,curx2).le.0.0) then
        pdf_cur=pdf(a1,b1,mx1,sx1,curx1,x01,t1)*pdf(a2,b2,mx2,sx2,curx2,x02,t2)
        if(pdf_cur.gt.pdf_max) then
            pdf_max=pdf_cur
            x1x1_temp=curx1
            x1x2_temp=curx2
        end if
    end if
end do


!if no suitable point was found
!if((x1x1_temp==0.0).and.(x1x2_temp==0.0)) then
!    write(*,*) 'No suitable point was found. Increase the number of sampled points in the configuration file and try again.'
!    return !exiting the subroutine
!end if

!if no suitable point was found
if((x1x1_temp==0.0).and.(x1x2_temp==0.0)) then
    call random_number(u1)
    do k=1,10000
        sampledThisIter=sampledThisIter+1
        totalSampled=totalSampled+1
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

        allSampledPoints(iter,sampledThisIter,1)=curx1
        allSampledPoints(iter,sampledThisIter,2)=curx2  !writing the sampled points to the file sampledPoints.txt
        !new point generated. checking if its okay with our prerequisites
        if(failureFunction(curx1,curx2).le.0.0) then
            x1x1_temp=curx1
            x1x2_temp=curx2
            exit
        end if
    end do
end if

!assigning the new initial points to its variable
x1x1=x1x1_temp
x1x2=x1x2_temp

!looking in the sampled points a point for x2
pdf_max=0
x2x1_temp=0.0
x2x2_temp=0.0
do k=1,sampledThisIter
    curx1=allSampledPoints(iter,k,1)
    curx2=allSampledPoints(iter,k,2)
    if((curx1.ne.x1x1).and.(curx2.ne.x1x2)) then  !excluding the already chosen point
        if(failureFunction(curx1,curx2).le.0.0) then  !searching only in the failure domain
            if((abs(curx1-x1x1).gt.0.3*sx1).and.(abs(curx2-x1x2).gt.0.3*sx2)) then  !searching only for points outside the cuboid
                pdf_cur=pdf(a1,b1,mx1,sx1,curx1,x01,t1)*pdf(a2,b2,mx2,sx2,curx2,x02,t2)
                if(pdf_cur.gt.pdf_max) then
                    pdf_max=pdf_cur
                    x2x1_temp=curx1
                    x2x2_temp=curx2
                end if
            end if
        end if
    end if
end do

!if no suitable point was found
if((x2x1_temp==0.0).and.(x2x2_temp==0.0)) then
    call random_number(u1)
    do k=1,10000
        sampledThisIter=sampledThisIter+1
        totalSampled=totalSampled+1
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

        allSampledPoints(iter,sampledThisIter,1)=curx1
        allSampledPoints(iter,sampledThisIter,2)=curx2  !writing the sampled points to the file sampledPoints.txt
        !new point generated. checking if its okay with our prerequisites
        if(failureFunction(curx1,curx2).le.0.0) then
            if((abs(curx1-x1x1).gt.0.9*sx1).and.(abs(curx2-x1x2).gt.0.9*sx2)) then  !searching only for points outside the cuboid of the 1st point sampled before
                x2x1_temp=curx1
                x2x2_temp=curx2
                exit
            end if
        end if
    end do
end if

x2x1 = x2x1_temp
x2x2 = x2x2_temp

open(10,file='data/rootPointsX1.txt',access='append')
write(10,*) x1x1, x1x2
close(10)
open(10,file='data/rootPointsX2.txt',access='append')
write(10,*) x2x1, x2x2
close(10)

open(11,file='data/sampledPoints.txt',access='append')
do j=1,samplingPoints
    write(11,*)  allSampledPoints(iter,j,1), allSampledPoints(iter,j,2)
end do
close(11) !sampledPoints.txt

iterationSamples(iter) = sampledThisIter


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
    write(40,*) '"data/sampledPoints.txt" u 1:2 title "Sampled points",\'
    write(40,*) '"data/rootPointsX1.txt" u 1:2 title "X1 root points",\'
    write(40,*) '"data/rootPointsX2.txt" u 1:2 title "X2 root points"'
close(40)

    call system('gnuplot data/plotSampling.plt')

end subroutine

real function h(omega1,omega2,x1x1,x1x2,x2x1,x2x2,sx1,sx2,x1,x2)
    real omega1,omega2,x1x1,x1x2,x2x1,x2x2,sx1,sx2,x1,x2
    real pdfn1, pdfn2, pdfn3, pdfn4
    real h_temp

    pdfn1=pdfn(x1x1,sx1,x1)
    pdfn2=pdfn(x1x2,sx2,x2)
    pdfn3=pdfn(x2x1,sx1,x1)
    pdfn4=pdfn(x2x2,sx2,x2)
    if(isnan(pdfn1)) then
        pdfn1=0.0
    end if
    if(isnan(pdfn2)) then
        pdfn2=0.0
    end if
    if(isnan(pdfn3)) then
        pdfn2=0.0
    end if
    if(isnan(pdfn4)) then
        pdfn2=0.0
    end if

    h_temp=omega1*pdfn1*pdfn2+omega2*pdfn3*pdfn4
    if(isnan(h_temp)) then
        h=0.0
    else
        h=t_temp
    end if

end function

subroutine h_test(omega1,omega2,x1x1,x1x2,x2x1,x2x2,sx1,sx2,x1,x2,test_val)
    real omega1,omega2,x1x1,x1x2,x2x1,x2x2,sx1,sx2,x1,x2
    real pdfn1, pdfn2, pdfn3, pdfn4
    real h_temp
    real test_val
    real h

    pdfn1=pdfn(x1x1,sx1,x1)
    pdfn2=pdfn(x1x2,sx2,x2)
    pdfn3=pdfn(x2x1,sx1,x1)
    pdfn4=pdfn(x2x2,sx2,x2)

    !write(*,*) "pdfn1:", pdfn1
    !write(*,*) "pdfn2:", pdfn2
    !write(*,*) "pdfn3:", pdfn3
    !write(*,*) "pdfn4:", pdfn4
    !write(*,*) "omega1:", omega1
    !write(*,*) "omega2:", omega2
    if(isnan(pdfn1)) then
        pdfn1=0.0
    end if
    if(isnan(pdfn2)) then
        pdfn2=0.0
    end if
    if(isnan(pdfn3)) then
        pdfn2=0.0
    end if
    if(isnan(pdfn4)) then
        pdfn2=0.0
    end if

    h_temp=omega1*pdfn1*pdfn2+omega2*pdfn3*pdfn4
    !write(*,*) "h_temp from test subroutine:", h_temp
    if(isnan(h_temp)) then
        h=0.0
    else
        h=t_temp
    end if

    test_val=h_temp

end subroutine


real function indicator(x1,x2)
    real x1, x2, failureFunction
    if(failureFunction(x1,x2).le.0.0) then
        indicator=1.0
    else
        indicator=0.0
    end if
end function
