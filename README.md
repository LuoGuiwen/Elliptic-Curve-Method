# Elliptic-Curve-Method
Elliptic Curve Method for integer factorization.

compile:
g++ -g -O2 ECM.cpp -o ECM2.0 -lgmpxx -lgmp

ATTENTION THAT the programme is just suitable for finding
factor whoes length is no more than 50-digit.
smoothness parameter B1 must be no more than 43000000,
B2 must be no more than 2.4e11.

RECOMMANDED PARAMETERS

digits  B1  		B2     			        d
15		2000		  1.5e5(147396)		    2*3*5*7*11;
20    11000    	1.9e6(1873422) 		  2*3*5*7*11;
25		50000		  1.3e7(12746592)   	2*3*5*7*11*13;
30		250000		1.3e8(128992510)  	2*3*5*7*11*13;
35  	1000000		1.0e9(1045563762)	  2*3*5*7*11*13*17;
40		3000000		1.0e9(1045563762)	  2*3*5*7*11*13*17;
45		3000000		5.3e9(5706890290) 	2*3*5*7*11*13*17;
45		1100000		3.5e10(35133391030)	2*3*5*7*11*13*17;
50 		43000000	2.4e11(240491351116)2*3*5*7*11*13*17*19;

--edit by GUIWEN LUO--
    --20180806--
