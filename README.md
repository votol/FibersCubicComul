This project is a calculation programm to model quantum properties of light propagation
along an array of single mode waveguides with cubic nonliniarity. Calculation is performed 
based on comulant approach. Project contains analytic computing part which allows to write a
closed system of ordinary differential equations based on Heizenberg approach.
The numeric solution of the given system of equetions is performed here by using 4th order
Runge–Kutta method.

To build code use 
    
    make
    
To run the result
    
    make run_test
    
or better use 
    
    make release
    
after which in workDir will appear a directory with name revN, where N is some number form 1.
In this directory will appear binary called CalcProj, it requaires .yaml file as it only input.
See "test" folder for example, also please see description.yaml file for description of what should be placed 
in "parameters" and "properties" section of this input .yaml file.

If you have any questions about usage of the code, please contact martvo@ipfran.ru
