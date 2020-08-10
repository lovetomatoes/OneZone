main: main.cpp evol.o class_gas.o LE_iso.o LE_adb.o read_aTree.o class_halo.o dyn.o thermo.o kpd.o reaction.o Newton.o my_linalg.o gsl_inverse.o PARA.o RK4.o
	g++ main.cpp -L/usr/local/lib evol.o class_gas.o LE_iso.o LE_adb.o read_aTree.o class_halo.o dyn.o thermo.o kpd.o reaction.o Newton.o my_linalg.o gsl_inverse.o RK4.o -lgsl -lgslcblas -lm -o main

# 用 make evol 生成cc.so以供evol.py 
evol.o: evol.cpp class_gas.o LE_iso.o LE_adb.o read_aTree.o class_halo.o dyn.o thermo.o kpd.o reaction.o Newton.o my_linalg.o gsl_inverse.o PARA.o RK4.o
	g++ evol.cpp -L/usr/local/lib class_gas.o LE_iso.o LE_adb.o read_aTree.o class_halo.o dyn.o thermo.o kpd.o reaction.o Newton.o my_linalg.o gsl_inverse.o RK4.o -lgsl -lgslcblas -lm -o cc.so -shared -fPIC
	g++ evol.cpp -L/usr/local/lib -Wall -I/usr/local/include -c class_gas.o LE_iso.o read_aTree.o class_halo.o dyn.o thermo.o kpd.o reaction.o Newton.o my_linalg.o gsl_inverse.o RK4.o -shared -fPIC

LE_iso.o: LE_iso.cpp RK4.o class_halo.o dyn.o PARA.o
	g++ -c LE_iso.cpp -o LE_iso.o -shared -fPIC 
	#g++ -c LE_iso.cpp RK4.o class_halo.o -o LE_iso.o -shared -fPIC


class_gas.o: class_gas.cpp class_halo.o kpd.o reaction.o Newton.o my_linalg.o gsl_inverse.o thermo.o dyn.o PARA.o read_aTree.o RK4.o
	g++ -c class_gas.cpp read_aTree.o -shared -fPIC

read_aTree.o: read_aTree.cpp PARA.o dyn.o LE_adb.o
	g++ -c read_aTree.cpp dyn.o LE_adb.o -shared -fPIC

LE_adb.o: LE_adb.cpp RK4.o class_halo.o dyn.o PARA.o
	g++ -c LE_adb.cpp -o LE_adb.o -shared -fPIC 

RK4.o: RK4.cpp my_linalg.o
	g++ -c RK4.cpp my_linalg.o -shared -fPIC


class_halo.o: class_halo.cpp PARA.o
	g++ -c class_halo.cpp -shared -fPIC

reaction.o: reaction.cpp PARA.o
	g++ -c reaction.cpp -shared -fPIC

kpd.o: kpd.cpp PARA.o
	g++ -c kpd.cpp -shared -fPIC

dyn.o: dyn.cpp PARA.o
	g++ -c dyn.cpp -shared -fPIC

thermo.o: thermo.cpp PARA.o
	g++ -c thermo.cpp -shared -fPIC

PARA.o: PARA.cpp
	g++ -c PARA.cpp -shared -fPIC

Newton.o: Newton.cpp gsl_inverse.o
	g++ -c  Newton.cpp gsl_inverse.o -shared -fPIC

my_linalg.o: my_linalg.cpp
	g++ -c my_linalg.cpp -shared -fPIC

gsl_inverse.o: gsl_inverse.cpp
	g++ -Wall -I/usr/local/include -c gsl_inverse.cpp -shared -fPIC
#	g++ -L/usr/local/lib gsl_inverse.o -lgsl -lgslcblas -lm

clean: # ./以免文件名前缀奇怪 e.g. "-" 不识别; -f force 忽略不存在的文件无警告
	rm -f ./*.o *.txt *.10 *.so *.out *.plt main