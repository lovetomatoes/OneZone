rese
se key autotitle columnhead
se log
set key right top font ',25,Bold'
#se ytics 0.5
#se mytics 2
set tics font ", 30,Bold"
set xtics offset 0,-1.5
# 坐标刻度格式
set format y "10^{%L}"
se format x "10^{%L}"
set xlabel "eV" font ",35,Bold"
#set ylabel "c" font ",35,Bold"

se xlabel offset 0,-2
se ylabel offset -10,0
#Y軸の余白
set lmargin 18
se rmargin 4
#X軸の余白
set bmargin 6
set tmargin 3

#se xrange[:13.6]
c = 3.e8; hp = 6.63e-27; eV = 1.6e-12; A = 1.e-10
ltoeV = c*hp/A/eV;
f = 'spec_send/rspec_pop3_cont_fe00.txt'
p f u (ltoeV/$1):2 w l lw 2
linef = 'spec_send/rspec_pop3_line_fe00.txt'
rep linef u (ltoeV/$1):2 w l lw 2
