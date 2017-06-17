program main
    use variables
    implicit none
    integer i
    call readFile

    !call testFunctions(a,b,mx,sx,q,x0,x)
    !call testFunction(a,b,mx,sx,q,x0,-0.546529412,t)

    !do i=1,5
     !  call generateCdfs(samplesSize,a,b,mx,sx,q,x0,i,barSize)
    !enddo

    !call generateCdfs(samplesSize,a,b,mx,sx,q,x0,t,barSize)


    call readFileNew

    call failureMonteCarlo

    call failureImportanceSampling
end
