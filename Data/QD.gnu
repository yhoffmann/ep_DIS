reset session

set terminal wxt size 500,500

set termoption font ',7'

set lmargin 5
set rmargin 5
set tmargin 5
set bmargin 5

set view 64,159

#set size square 1.0,1.0
set key top right box

set grid

set xlabel "Q"
set ylabel "|{/Symbol D}|"
set zlabel "d{/Symbol s}/dt-longitudinal [nb]" rotate parallel offset -3.5

set title "Q and Delta Tests"

set pm3d

#set palette defined (0 "white", 0.5 "#ffd000" , 1.6 "red", 2.7 "blue", 3.8 "black")

set logscale z 10

#splot "DeltaAndQ_Coherent_Test.txt" using 1:2:3 notitle # transverse
splot "DeltaAndQ_Coherent_Test.txt" using 1:2:($4*1e+7) notitle # longitudinal
