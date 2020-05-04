This Matlab code was developed for the paper:

Peter C., Salina Borello E., Dalman R., Karamitopoulos P., Busschers F., Sacchi Q., Verga F. (2919) 
"Improved lithology prediction in channelized reservoirs by integrating stratigraphic forward modelling: 
towards improved model calibration in a case study of the Holocene Rhine-Meuse fluvio-deltaic system"
Computers and Geosciences. Elsevier (in press)

The code implement the Neighborhood algorithm to be applied to Geological basin inversion, based on:
- Sambridge, M. (1999). Geophysical inversion with a neighborhood algorithm–Searching a parameters space: Geophysical Journal International, 138, 479–494, doi:10.1046/j.1365-246X.1999.00876.x.
- Sambridge, M. (2001). Finding acceptable models in nonlinear inverse problems using a neighborhood algorithm: Inverse Problems, 17, 387–403, doi:10.1088/0266-5611/17/3/302.
- Sambridge, M., and K. Mosegaard (2002). Monte Carlo methods in geophysical inverse problems: Reviews of Geophysics, 40, 1–29. doi: 10.1029/2000RG000089.

The code is given with a set of standard test functions('Multi','Holder', 'Goldstein-Price','Booth')
but can be applied to custom functions provided that the user defines:
 - ra : range of each variable, ex. ra{1}=[min_x1, max_x1];
 - fo : fitness function, written as anonymous function
 - tol_cal: the fo value under which the corresponding variable set is considered as calibrated
 - header : the header of the output file 'population.txt'
 - formato_out : the printing format of the output file 'population.txt'

To launch the code with one of the proposed test function, just uncomment the corresponding line in the code. 
Ex. 

[header, formato_out,ra,fo,tol_cal]= input_test('Multi');
%[header, formato_out,ra,fo,tol_cal]= input_test('Holder');
%[header, formato_out,ra,fo,tol_cal]= input_test('Goldstein-Price');
%[header, formato_out,ra,fo,tol_cal]= input_test('Booth');
