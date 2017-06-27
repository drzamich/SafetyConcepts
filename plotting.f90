subroutine plotImportanceSampling(iter)
    implicit none
    integer iter
    real x1max,x1min,x2max,x2min
    real c1_r,c0

open(10,file='data/chartRange.txt')
read(10,*) x1max
read(10,*) x1min
read(10,*) x2max
read(10,*) x2min
read(10,*) c1_r
read(10,*) c0
close(10)

open(40,file='data/plotImportanceSampling.plt')
    write(40,*) 'set terminal pngcairo enhanced font "Verdana,10"'
    write(40,'(A,I6.6,A)') 'set output "chartsSampling/sampling',iter,'.png"'
    write(40,*) 'set style line 1 lt 1 lc rgb "#0000FF" lw 1'
    write(40,*) 'set xrange [',x1min-1.1*abs(x1min),':',x1max+1.1*abs(x1max),']'
    write(40,*) 'set yrange [',x2min-1.1*abs(x2min),':',x2max+1.1*abs(x2max),']'
    write(40,*) 'g(x) = ',c1_r,'*x-',c0
    write(40,*) 'set xlabel "X1"'
    write(40,*) 'set ylabel "X2"'
    write(40,*) 'plot g(x) title "g(x1,x2)=0" ls 1,\'
    write(40,*) '"data/sampledPoints.txt" u 1:2 title "Sampled points",\'
    write(40,*) '"data/rootPointsX1.txt" u 1:2 title "X1 root points",\'
    write(40,*) '"data/rootPointsX2.txt" u 1:2 title "X2 root points",\'
    write(40,*) '"data/exVal.txt" u 1:2 title "Expected value"'
close(40)

    call system('gnuplot data/plotImportanceSampling.plt')

end subroutine

subroutine plotMonteCarlo(i)
    implicit none
    integer i
    real x1max,x1min,x2max,x2min,c1_r,c0

open(10,file='data/chartRange.txt')
read(10,*) x1max
read(10,*) x1min
read(10,*) x2max
read(10,*) x2min
read(10,*) c1_r
read(10,*) c0
close(10)

open(10,file='data/plotMonteCarlo.plt')
    write(10,*) 'set terminal pngcairo enhanced font "Verdana,10"'
    write(10,*) 'set style line 1 lt 1 lc rgb "#0000FF" lw 1'
    write(10,'(A,I7.7,A)') 'set output "chartsMonteCarlo/failure',i,'.png"'
    write(10,*) 'g(x) = ',c1_r,'*x-',c0
    write(10,*) 'set xlabel "X1"'
    write(10,*) 'set ylabel "X2"'
    write(10,*) 'set xrange [',x1min-0.2*abs(x1min),':',x1max+0.2*abs(x1max),']'
    write(10,*) 'set yrange [',x2min-0.2*abs(x2min),':',x2max+0.2*abs(x2max),']'
    write(10,*) 'plot g(x) title "g(x1,x2)=0" ls 1,\'
    write(10,*) '"data/pointsSurvival.txt" u 1:2 title "Survival",\'
    write(10,*) '"data/pointsFailure.txt" u 1:2 title "Failure",\'
    write(10,*) '"data/exVal.txt" u 1:2 title "Expected value"'
    close(10)

    call system('gnuplot data/plotMonteCarlo.plt')
    write(*,*) "chart plotted for i=",i
end subroutine

subroutine plotConvergence(maxProbability)
    implicit none
    real maxProbability

    open(40,file='data/plotConvergence.plt')
    write(40,*) 'set terminal pngcairo enhanced font "Verdana,10"'
    write(40,*) 'set output "chartsMonteCarlo/convergence.png"'
    write(40,*) 'set xlabel "Sampled points"'
    write(40,*) 'set ylabel "Pf"'
    write(40,*) 'set yrange [0:',maxProbability*1.1,']'
    write(40,*) 'plot "data/convergence.txt" u 1:2 title "Failure Probability" with lines'
    close(40)

    call system('gnuplot data/plotConvergence.plt')
    !write(*,*) "Convergence plotted"
end subroutine

