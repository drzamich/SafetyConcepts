set title "Empiric vs theoretical CDF"
set xlabel "x"
set ylabel "F(x)"
set terminal png
set style data lines
set output "cdfs.png"
plot "empiricCDF.txt" title "Empiric CDF", "theoryCDF.txt" title "Theoretical CDF"
