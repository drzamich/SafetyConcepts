subroutine testFunctions(a,b,mx,sx,q,x0,x)
    implicit none

    real a,b,mx,sx,q,x0,x
    real cdf,inv_cdf
    integer t


    t=1
    write(*,*) "Uniform distribution"
    write(*,*) "Value of CDF:", cdf(a,b,mx,sx,x,x0,t)
    write(*,*) "Value of CDF-1:", inv_cdf(a,b,mx,sx,q,x0,t)
    write(*,*)

    t=2
    write(*,*) "Standard normal distribution"
    write(*,*) "Value of CDF:", cdf(a,b,mx,sx,x,x0,t)
    write(*,*) "Value of CDF-1:", inv_cdf(a,b,mx,sx,q,x0,t)
    write(*,*)

    t=3
    write(*,*) "Normal distribution"
    write(*,*) "Value of CDF:", cdf(a,b,mx,sx,x,x0,t)
    write(*,*) "Value of CDF-1:", inv_cdf(a,b,mx,sx,q,x0,t)
    write(*,*)

    t=4
    write(*,*) "Log-normal distribution"
    write(*,*) "Value of CDF:", cdf(a,b,mx,sx,x,x0,t)
    write(*,*) "Value of CDF-1:", inv_cdf(a,b,mx,sx,q,x0,t)
    write(*,*)

    t=5
    write(*,*) "Exponential distribution"
    write(*,*) "Value of CDF:", cdf(a,b,mx,sx,x,x0,t)
    write(*,*) "Value of CDF-1:", inv_cdf(a,b,mx,sx,q,x0,t)
    write(*,*)


end subroutine

subroutine testFunction(a,b,mx,sx,q,x0,x,t)
    implicit none

    real a,b,mx,sx,q,x0,x
    real cdf,inv_cdf
    integer t

    write(*,*) "Value of CDF:", cdf(a,b,mx,sx,x,x0,t)
    write(*,*) "Value of CDF-1:", inv_cdf(a,b,mx,sx,q,x0,t)

end subroutine

