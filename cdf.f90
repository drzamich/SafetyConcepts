real function cdf(a,b,mx,sx,t,x,x0)
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
    end select


end function


real function inv_cdf(a,b,mx,sx,t,q,x0)
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
    end select

end function

