reset session

set terminal wxt size 600,600

set termoption font ',8'

set lmargin 2
set rmargin 2
set tmargin 2
set bmargin 2

set view 64,159

set size square 0.98,0.98
set key top right box

set grid

set xlabel "Q"
set ylabel "|{/Symbol D}|"
set zlabel "d{/Symbol s}/dt-longitudinal [nb]" rotate parallel offset -3.5

set title "Q and Delta Tests"

set pm3d #map interpolate 20,20

#set palette defined (0 "white", 0.5 "#ffd000" , 1.6 "red", 2.7 "blue", 3.8 "black")
#set palette rgb 21,22,23

#set logscale z 10

#splot "DeltaAndQ_Coherent_Test.txt" using 1:2:3 notitle # transverse
splot "DeltaAndQ_Coherent_Test.txt" using 1:2:($4/$6) notitle # longitudinal