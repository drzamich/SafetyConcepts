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
            pdf=a*exp(-1*a*(x-b)-exp(-1*a*(x-b)))
    end select

end function

real function pdfn(mx,sx,x)
    implicit none
    real mx,sx,x,pi
    pi=3.1415926535
    pdfn=(1/(sx*sqrt(2*pi)))*exp(-0.5*((x-mx)/sx)**2)
end function
