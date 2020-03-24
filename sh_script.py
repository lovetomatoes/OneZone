import os

""" os.system("g++ -Wall -I/usr/local/include -c gsl_inverse.cpp")
os.system("g++ -L/usr/local/lib gsl_inverse.o -lgsl -lgslcblas -lm")
os.system("./a.out")
os.system("python gsl_inverse.py") """

""" os.system("g++ -c pointer.cpp")
os.system("g++ pointer.o")
os.system("./a.out") """

os.system("g++ -c reaction.cpp")
os.system("g++ class_gas.cpp reaction.o -o aaa")
os.system("./aaa")