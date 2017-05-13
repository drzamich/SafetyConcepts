program main
    use variables
    implicit none

    call readFile

    !call testFunctions(a,b,mx,sx,q,x0,x)
    !call testFunction(a,b,mx,sx,q,x0,x,t)
    call generateCdfs(samplesSize,a,b,mx,sx,q,x0,t)
end
