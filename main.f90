program main
    use variables
    implicit none
    real phi

    write(*,*) phi(-3.898015)
    call readFile

    !call failureMonteCarlo

    call failureImportanceSampling

end
