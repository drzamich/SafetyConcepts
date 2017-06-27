real function cdf(a,b,mx,sx,x,x0,t)
    implicit none

    real x,x0
    real a,b,mx,sx
    real mu,su
    integer t
    real phi


    select case(t)
        case(1)      !uniform distribution
            cdf = (x-a)/(b-a)
        case(2)      !standard normal distribution
            cdf=phi(x)
        case(3)      !normal distribution
            cdf=phi((x-mx)/sx)
        case(4)      !log-normal distribution
            su=sqrt(log((sx/(mx-x0))**2+1))
            mu=log(mx-x0)-(su**2)*0.5
            cdf=phi((log(x-x0)-mu)/su)
        case(5)       !exponential distribution
            cdf=1-exp(-1*b*(x-x0))
		case(6) 		!ex-max type 1
			cdf=exp(-1*exp(-1*a*(x-b)))
    end select


end function


real function inv_cdf(a,b,mx,sx,q,x0,t)
    implicit none

    real x,q
    real phi_inv
    real a,b,x0
    real mx,sx,mu,su
    integer t

    select case(t)
        case(1)      !uniform distribution
            inv_cdf=a+q*(b-a)
        case(2)      !standard normal distribution
            inv_cdf=phi_inv(q)
        case(3)      !normal distribution
            inv_cdf=mx+sx*phi_inv(q)
        case(4)      !log-normal distribution
            su=sqrt(log((sx/(mx-x0))**2+1))
            mu=log(mx-x0)-(su**2)*0.5
            inv_cdf=x0+exp(mu+su*phi_inv(q))
        case(5)       !exponential distribution
			inv_cdf=x0-(1/b)*log(1-q)
		case(6) !ex-max type 1 distribution
			a=(3.141/sqrt(6.0))*(1/sx)
			b=mx-(0.577216/a)
			inv_cdf=b-(1/a)*log(-1*log(q))
    end select

end function


real function pdf(a,b,mx,sx,x,x0,t)
    implicit none

    real x,q
    real phi_inv, phi_small
    real a,b,x0
    real mx,sx,mu,su
    integer t
    real pi

    pi=3.1415926535

    select case(t)
        case(1)      !uniform distribution
            pdf=1/(b-a)
        case(2)      !standard normal distribution
            pdf=phi_small(x)
        case(3)      !normal distribution
            pdf=(1/(sx*sqrt(2*pi)))*exp(-0.5*((x-mx)/sx)**2)
        case(4)      !log-normal distribution
            su=sqrt(log((sx/(mx-x0))**2+1))
            mu=log(mx-x0)-(su**2)*0.5
            pdf=(1/(su*(x-x0)))*phi_small((log(x-x0)-mu)/(su))
        case(5)       !exponential distribution
            pdf=b*exp(-b*(x-x0))
        case(6) !ex-max type I
            a=(3.141/sqrt(6.0))*(1/sx)
			b=mx-(0.577216/a)
            pdf=a*exp(-1*a*(x-b)-exp(-1*a*(x-b)))
    end select

end function

real function pdfn(mx,sx,x)
    implicit none
    real mx,sx,x,pi
    pi=3.1415926535
    pdfn=(1/(sx*sqrt(2*pi)))*exp(-0.5*((x-mx)/sx)**2)
end function


real function inv_cdfn(mx,sx,q)
    real mx,sx,q
    inv_cdfn=mx+sx*phi_inv(q)
end function


real function failureFunction(x1,x2)
    use variables
	implicit none
	real x1,x2

failureFunction = c0+c1*x1+c2*x2

end function


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

    !some changes
    !some more changes
    !even more chang
    !more more more
end function

real function phi_small(x)
    implicit none
    real x
    real pi
    pi=3.1415926535

    phi_small=(1/sqrt(2*pi))*exp((-1)*(x**2)*0.5)
end function
