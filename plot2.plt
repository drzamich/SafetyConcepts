set title "Theoretical PDF vs historgram"
set xlabel "x"
set ylabel "f(x)"
set terminal png
set style data lines
set output "pdfs".png"
plot "theoreticalPDF.txt" title "Theoretical PDF"

