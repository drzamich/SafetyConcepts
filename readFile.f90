subroutine readFile
    use variables
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
