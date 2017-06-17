module variables2
    real a,b,mx,sx,q,x0,x
    real cdf,inv_cdf
    integer t
    real test,mu,su,phi_inv,phi
    real pi

    integer samplesType, samplesSize
    real barSize
end module


module variables
    real a1,b1,mx1,sx1,x01
    real a2,b2,mx2,sx2,x02
    real eps
    integer maxIter, maxPoints
    real c1,c2 !parameters of the failure function
    real x1_init
    real c1_r !paramater next to x1 in a rearanged equation c1*x1+c2*x2=0
    integer samplingPoints

    integer t1,t2

end module
