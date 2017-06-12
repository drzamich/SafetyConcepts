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

real x1x1, x1x2, x2x1, x2x2

samplingPoints=20

count=0
size1=1000
charts=5


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

x1x1=init1x1
x1x2=init1x2
x2x1=init2x1
x2x2=init2x2

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
end subroutine

real function h(omega1,omega2,x1x1,x1x2,x2x1,x2x2,mx1_x1,mx2_x1,mx1_x2,mx2_x2,sx1,sx2)
    real omega1,omega2,x1x1,x1x2,x2x1,x2x2,mx1_x1,mx2_x1,mx1_x2,mx2_x2,sx1,sx2

    h=omega1*pdfn(mx1_x1,sx1,x1x1)*pdfn(mx2_x1,sx2,x1x2)+&
      omega2*pdfn(mx1_x2,sx1,x2x1)*pdfn(mx2_x2,sx2,x2x2)
end function
