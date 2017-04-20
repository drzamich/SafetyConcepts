subroutine testFunctions
    implicit none
    real a,b,mx,sx,q,x0,x
    real cdf,inv_cdf
    integer t
    real test,mu,su,phi_inv,phi

    a=-500
    b=5
    mx=0.6
    sx=0.3
    t=1
    x=0.6
    x0=0.5
    q=0.4


    t=1
    write(*,*) "Uniform distribution"
    write(*,*) "Value of CDF:", cdf(a,b,mx,sx,t,x,x0)
    write(*,*) "Value of CDF-1:", inv_cdf(a,b,mx,sx,t,q,x0)
    write(*,*)

    t=2
    write(*,*) "Standard normal distribution"
    write(*,*) "Value of CDF:", cdf(a,b,mx,sx,t,x,x0)
    write(*,*) "Value of CDF-1:", inv_cdf(a,b,mx,sx,t,q,x0)
    write(*,*)

    t=3
    write(*,*) "Normal distribution"
    write(*,*) "Value of CDF:", cdf(a,b,mx,sx,t,x,x0)
    write(*,*) "Value of CDF-1:", inv_cdf(a,b,mx,sx,t,q,x0)
    write(*,*)

    t=4
    write(*,*) "Log-normal distribution"
    write(*,*) "Value of CDF:", cdf(a,b,mx,sx,t,x,x0)
    write(*,*) "Value of CDF-1:", inv_cdf(a,b,mx,sx,t,q,x0)
    write(*,*)

    t=5
    write(*,*) "Exponential distribution"
    write(*,*) "Value of CDF:", cdf(a,b,mx,sx,t,x,x0)
    write(*,*) "Value of CDF-1:", inv_cdf(a,b,mx,sx,t,q,x0)
    write(*,*)

end subroutine
