subroutine failureImportanceSampling
use variables
implicit none

real u
real inv_cdf, pdf
real failureFunction
real init1x1, init1x2, init2x1, init2x2, temp1,temp2
real omega1, omega2
integer samplingPointsx1, samplingPointsx2
real inv_cdfn
integer i,j,k,l
integer iter
real x1x1, x1x2, x2x1, x2x2, x1x1_temp,x1x2_temp,x2x1_temp,x2x2_temp
real pdf_cur, pdf_max
real curx1, curx2
real u1
real part1, part2, part3
real indicator, indicator_val, pdf_xi1, pdf_xi2
real h_val
real x1max,x1min, x2max,x2min

real, allocatable:: iterationSamples(:)  !array which keeps number of sampled points in each iteration
real, allocatable:: iterationProbability(:) !array which keeps the value of probability in each iteration
real, allocatable:: allSampledPoints(:,:,:)
real, allocatable:: h_params(:,:)

real, allocatable:: sampledPoints(:,:)
integer totalSampled
integer sampledThisIter
real previousProbability

integer passInARow

allocate(iterationSamples(maxIter))
allocate(iterationProbability(maxIter))
allocate(allSampledPoints(maxIter,maxPoints,2))
allocate(h_params(maxIter,6))

write(*,*) "Calculating the failure according to Importance Sampling Apporach"
iterationSamples=0.0
iterationProbability=0.0
allSampledPoints=0.0
totalSampled=0
passInARow=0

!initial points
x1x1=x1_init
x1x2=c1_r*0.99*x1x1
x2x1=x1x1+0.9*sx1
x2x2=c1_r*0.99*x2x1

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



call plotSampling(0,c1_r)

!
!
!  START OF THE ITERATION
!
!

do iter=1,maxIter
sampledThisIter=0

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

!therefore the choice of the starting point is important

samplingPointsx1=nint(samplingPoints*omega1)  !rounds to the nearest whole number
samplingPointsx2=nint(samplingPoints*omega2)


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


if(iter==1) then
    part1=0.0
do i=1,sampledThisIter
    j=1
        call h(h_params(j,1),h_params(j,2),h_params(j,3),h_params(j,4),h_params(j,5),h_params(j,6),&
                        sx1,sx2,allSampledPoints(u,i,1),allSampledPoints(u,i,2),h_val)
        pdf_xi1=pdf(a1,b1,mx1,sx1,allSampledPoints(1,i,1),x01,t1)
        pdf_xi2=pdf(a2,b2,mx2,sx2,allSampledPoints(1,i,2),x02,t2)
        indicator_val=indicator(allSampledPoints(1,i,1),allSampledPoints(1,i,2))
        part1=part1+(indicator_val*pdf_xi1*pdf_xi2)/h_val
end do

iterationProbability(iter) = part1/sampledThisIter

end if


if(iter.ne.1) then

part3=0.0

do u=1,iter
        part2=0
        do i=1,iterationSamples(u)
        part1=0.0
        do j=1,iter
            call h(h_params(j,1),h_params(j,2),h_params(j,3),h_params(j,4),h_params(j,5),h_params(j,6),&
            sx1,sx2,allSampledPoints(u,i,1),allSampledPoints(u,i,2),h_val)
            part1=part1+(iterationSamples(j)/totalSampled)*h_val
        end do
        indicator_val=indicator(allSampledPoints(u,i,1),allSampledPoints(u,i,2))
        pdf_xi1=pdf(a1,b1,mx1,sx1,allSampledPoints(u,i,1),x01,t1)
        pdf_xi2=pdf(a2,b2,mx2,sx2,allSampledPoints(u,i,2),x02,t2)

        if(isnan(pdf_xi1)) pdf_xi1=0.0
        if(isnan(pdf_xi2)) pdf_xi2=0.0

        part2=part2+((indicator_val*pdf_xi1*pdf_xi2)/part1)
    end do
    part3=part3+part2
end do

iterationProbability(iter) = part3/totalSampled

end if

!write(*,*) "probability: ", iterationProbability(iter)

previousProbability=iterationProbability(iter-1)
if((previousProbability-eps*previousProbability).lt.iterationProbability(iter)&
    .and.(previousProbability+eps*previousProbability).gt.iterationProbability(iter)) then
    passInARow=passInARow+1
    if(passInARow==5) then
        write(*,*) "Failure probability:", iterationProbability(iter)
        write(*,*) "Sampled points: ", totalSampled
        exit
    end if
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
            if((abs(curx1-x1x1).gt.0.9*sx1).and.(abs(curx2-x1x2).gt.0.9*sx2)) then  !searching only for points outside the cuboid
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


call plotSampling(iter,c1_r)
end do
end subroutine

subroutine plotSampling(iter)
    implicit none
    integer iter
    real x1max,x1min,x2max,x2min
    real c1_r

open(10,file='data/chartRange.txt')
read(10,*) x1max
read(10,*) x1min
read(10,*) x2max
read(10,*) x2min
read(10,*) c1_r
close(10)

open(40,file='data/plotSampling.plt')
    write(40,*) 'set terminal pngcairo enhanced font "Verdana,10"'
    write(40,'(A,I6.6,A)') 'set output "data/sampling',iter,'.png"'
    write(40,*) 'set style line 1 lt 1 lc rgb "#0000FF" lw 1'
    write(40,*) 'set xrange [',x1min-0.1*abs(x1min),':',x1max+0.1*abs(x1max),']'
    write(40,*) 'set yrange [',x2min-0.1*abs(x2min),':',x2max+0.1*abs(x2max),']'
    write(40,*) 'g(x) = ',c1_r,'*x'
    write(40,*) 'set xlabel "X1"'
    write(40,*) 'set ylabel "X2"'
    write(40,*) 'plot g(x) title "g(x1,x2)=0" ls 1,\'
    write(40,*) '"data/sampledPoints.txt" u 1:2 title "Sampled points",\'
    write(40,*) '"data/rootPointsX1.txt" u 1:2 title "X1 root points",\'
    write(40,*) '"data/rootPointsX2.txt" u 1:2 title "X2 root points"'
close(40)

    call system('gnuplot data/plotSampling.plt')

end subroutine


subroutine h(omega1,omega2,x1x1,x1x2,x2x1,x2x2,sx1,sx2,x1,x2,res)
    real omega1,omega2,x1x1,x1x2,x2x1,x2x2,sx1,sx2,x1,x2
    real pdfn1, pdfn2, pdfn3, pdfn4
    real res, res_temp

    pdfn1=pdfn(x1x1,sx1,x1)
    pdfn2=pdfn(x1x2,sx2,x2)
    pdfn3=pdfn(x2x1,sx1,x1)
    pdfn4=pdfn(x2x2,sx2,x2)

    if(isnan(pdfn1)) pdfn1=0.0
    if(isnan(pdfn2)) pdfn2=0.0
    if(isnan(pdfn3)) pdfn2=0.0
    if(isnan(pdfn4)) pdfn2=0.0

    res_temp=omega1*pdfn1*pdfn2+omega2*pdfn3*pdfn4

    if(isnan(h_temp)) then
        res_temp=0.0
    end if

    res=res_temp

end subroutine


real function indicator(x1,x2)
    real x1, x2, failureFunction
    if(failureFunction(x1,x2).le.0.0) then
        indicator=1.0
    else
        indicator=0.0
    end if
end function
