reset session

set terminal wxt size 600,600

set termoption font ',8'

#set view 64,159

set size square 0.98,0.98
set key top right box

set grid

set title "Longitudinal cross section"

set xlabel "|t| (={/Symbol D}²)"
set ylabel "d{/Symbol s}/dt-longitudinal [nb GeV^{-2}]"

#set xrange [0.0:2.0]


#set palette defined (0 "white", 0.5 "#ffd000" , 1.6 "red", 2.7 "blue", 3.8 "black")
#set palette rgb 21,22,23

set logscale y 10
#set pm3d
#splot "FirstOrderResultsAnalyticalVsNumerical.txt" using 1:2:($4/$6) notitle# pt 2 ps 2 lw 2 title "a"

#plot for [col=0:8:1] "DeltaAndQ_Coherent_Test.txt" using 2:(int($0)%col==1?$4:1/0) title "strcol($1)" with lines

plot for [IDX=0:4] 'FirstOrderResultsAnalyticalVsNumerical.txt' i IDX u ($2**2):($6) lw 2  title columnheader(1)
#plot [1.8:2.5] for [IDX=0:1] 'FirstOrderResultsAnalyticalVsNumerical.txt' i IDX u 2:($4/$6) lw 2  title columnheader(1)
