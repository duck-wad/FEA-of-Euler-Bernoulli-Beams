# FEA of Euler-Bernoulli Beams

- This repository contains simple code to create a beam with boundary conditions and loads, and analyze the displacements and reactions using the Finite Element Method for Euler-Bernoulli beams.
- This was made as a simple exercise to improve my coding skills and learn about the FEA process!
- Python scripts are used to generate a text input file which is read by a C++ program to perform the FEA. Python code then reads the C++ output from a CSV and performs post-processing exercises like plotting deformations and making BMD/SFD.

## Generating the Input
1. In the run.py script, create a beam with a specified length and material parameters.
2. Add boundary conditions using the "add_BC" class method, and specifying a type (ex. clamp, pin) and a location on the beam.
3. Apply loads using the "add_load" class method, specifying a type (ex. distributed, moment, force), location, and magnitude. For distributed loads, start and end location and magnitude is required.
4. Discretize the beam into elements using the "discretize" method. Input an approximate element length, and the program will discretize the beam into approximately those sizes depending on the location of the nodes defined in step 2 and 3.
5. Run the "create_infile" method to generate the input file.
6. Run the "run" method to call the C++ executable.

## Performing FEA
- The C++ program reads the "INPUT.txt" file generated from the Python script. 
- The mesh is discretized according to the node and element input from the infile. Elemental stiffness matrices and force vectors are developed and assembled into the global stiffness matrix and global force vector.
- The system of equations is solved using Gaussian Elimination.
- Currently nodal displacements and reactions are output from the program. In the future, moment and shear will also be output.

## Displaying the Output
- The run.py script calls "interpreter.py" which runs the post-processing activities.
- The deformed shape of the beam is plotted along with the calculated reactions using matplotlib.
- BMD and SFD are constructed as derivatives of a CubicHermiteSpline of the nodal displacement. 
- Results are exported to an output Excel file for further analysis if needed.

### Example Beam Deformed Shape
![image](https://github.com/user-attachments/assets/b77f1da0-0b13-4e36-a0ae-d9ec9a6efc89)

### Example BMD/SFD
![image](https://github.com/user-attachments/assets/6068127e-f01d-412e-b497-154ba712e519)
