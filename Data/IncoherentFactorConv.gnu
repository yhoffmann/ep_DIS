reset session

set terminal wxt size 600,600

set termoption font ',8'

#set view 64,159

set size square 0.98,0.98
set key bottom left box

set grid

set title "Cross section first order analytical divided by first order numerical result"

set xlabel "|{/Symbol D}|"
set ylabel "d{/Symbol s}/dt-longitudinal [nb GeV^{-2}]"

#set xrange [1.8:2.5]
#set yrange [1e-5:100]


#set palette defined (0 "white", 0.5 "#ffd000" , 1.6 "red", 2.7 "blue", 3.8 "black")
#set palette rgb 21,22,23

#set logscale y 10
#set pm3d
#splot "FirstOrderResultsAnalyticalVsNumerical.txt" using 1:2:($4/$6) notitle# pt 2 ps 2 lw 2 title "a"

#plot for [col=0:8:1] "DeltaAndQ_Coherent_Test.txt" using 2:(int($0)%col==1?$4:1/0) title "strcol($1)" with lines

#plot [0:1.4] "FirstOrderResultsAnalyticalVsNumerical.txt" using 2:($4/$6) with lines
plot 'FirstTermDividedBySecondTermRange10.txt' u ($2**2.0):($3) lw 2 w lines title "sqrt(10)"
replot 'FirstTermDividedBySecondTermRange20.txt' u ($2**2.0):($3) lw 2 w lines title "sqrt(20)"
replot 'FirstTermDividedBySecondTermRange40.txt' u ($2**2.0):($3) lw 2 w lines title "sqrt(40)"
replot 'FirstTermDividedBySecondTermRange80.txt' u ($2**2.0):($3) lw 2 w lines title "sqrt(80)"
