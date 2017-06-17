subroutine readFile
    use variables2
    implicit none

!Reading values from input file
open (10, file='input.txt')

read(10,*)
read(10,*) a
read(10,*) b
read(10,*) mx
read(10,*) sx
read(10,*) x
read(10,*) x0
read(10,*) q
read(10,*) t

read(10,*)
read(10,*)
read(10,*) samplesSize
read(10,*) barSize
close(10)


end subroutine

subroutine readFileNew
    use variables
    implicit none

!Reading values from input file
open (10, file='inputFailure.txt')

read(10,*)
read(10,*) t1
read(10,*) a1
read(10,*) b1
read(10,*) mx1
read(10,*) sx1
read(10,*) x01


read(10,*)
read(10,*)
read(10,*) t2
read(10,*) a2
read(10,*) b2
read(10,*) mx2
read(10,*) sx2
read(10,*) x02

read(10,*)
read(10,*) !properties of monte carlo apporach
read(10,*) samplingMonteCarlo
read(10,*) chartsMonteCarlo

read(10,*)
read(10,*) !properties of adaptive importance apporach
read(10,*) samplingPoints
read(10,*) x1_init
read(10,*) eps
read(10,*) maxIter
read(10,*) maxPoints

read(10,*)
read(10,*) !properties of failure function
read(10,*) c1
read(10,*) c2

close(10)


c1_r = -(c1/c2)

end subroutine
