recursive real function phi(x) result(phi_res)
    implicit none

    real x
    real p,b1,b2,b3,b4,b5
    real t
    real z
    real pi

    p=0.2316419
    b1=0.319381530
    b2=-0.356563782
    b3=1.781477937
    b4=-1.821255978
    b5=1.330274429
    pi=3.1415926535

    t=1.0/(1+p*x)


    if (x.ge.0.0) then
        z=((1.0)/(sqrt(2*pi)))*exp((x*x*(-1))/2)
        phi_res=1-z*(b1*t+b2*t**2+b3*t**3+b4*t**4+b5*t**5)
    else if (x.lt.0.0) then
        phi_res=1-phi(x*(-1))
    end if
end function

recursive real function phi_inv(q) result(phi_inv_res)
    implicit none

    real q
    real u
    real c0,c1,c2
    real d1,d2,d3

    c0=2.515517
    c1=0.802853
    c2=0.010328

    d1=1.432788
    d2=0.189269
    d3=0.001308

    if((q.gt.0).and.(q.le.0.5)) then
        u=sqrt(log(q**(-2)))
        phi_inv_res=(-1)*u+(c0+c1*u+c2*u**2)/(1+d1*u+d2*u**2+d3*u**3)
    else if((q.gt.0.5).and.(q.le.1.0)) then
        phi_inv_res=(-1)*phi_inv(1-q)
    end if
end function

