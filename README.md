Low Temperature Air Plasma README

    Prerequisites:
        NOTE - The low temperature air plasma model was developed in Eilmer Version 4.1.0, if you wish to convert the files to Eilmer Version 5.0.0 (lmr), you will need to make the conversions yourself (sorry)
        It is presumed you have Eilmer installed on your device and all of its prerequisites, the below prerequisites are additional to the base Eilmer model:
            - Gas-dynamic library
                $ cd gdtk/src/gas 
                $ make install

    To run the unit tests:
        Copy the kinetics file into gdtk/src/kinetics
        Copy the gas file into gdtk/src/gas
        Copy the .lua files in this repository (in sample-data and sample-input) into their respective gas and kinetics sample-data and sample-input directories, the unit tests will refer to them. Ensure you choose your kinetics ODE method, only Forward Euler and Runge Kutta 4 are implemented, a Backward Euler was not implemented purely since time was of the essence in this project and we have a state vector with 12 variables
        Choose your case (relaxation or energy investigation)
        In your linux environment, move to the kinetics directory ($ cd gdtk/src/kinetics) and run all kinetics tests ($ make test) which will run through all of the kinetics file tests

    To use the plotting scripts:
        If you used the relaxation test case, run the low_temp_air_plasma_relaxation_plotting.py script
        If you used the energy investigation test case, run the low_temp_air_plasma_energy_plotting.py script
    
    Version updates:
        Version 1.0.0
            - Very initial kinetics and gas module with the very first few unit test results
        Version 1.1.0
            - The slightly cleaned up and improved kinetics file which fixed some typos and human error when it came to the hard coding of the reaction rates
        Version 1.2.0
            - Implemented a better management of the two energy modes in the kinetics file
        Version 1.2.1
            - Implemented the correct reference vibrational energy in the gas model and kinetics file, did not attach the corresponding figures as they were unpleasant
        Version 2.0.0
            - Optional total energy which is influenced from the electron beam power Qb was added to better inform the simple unit test cases I was doing to present reliable results
              - This needs to be handled with care, the PFE jobs and source term files account for this additional energy from the electron beam power themselves, so any points where Qb is added to the total energy should have the Qb term removed. There should be enough comments in my kinetics file to direct users to do this
            - Better care has been taken when calculating the number density of NO+, ensuring we don't get any negative values
            - The GasState is now updated after each step in the Forward Euler and Runge Kutta 4 methods to improve stability
            - Changed the elastic energy exchange method to better reflect the Appleton-Bray expression as defined by Nick Gibbons in two_temperature_air_kinetics.d
            - The elastic energy exchange was kept outside of the state vector derivatives, i.e. an external source, since the electron temperature is maintained by the reduced electric field strength rather than being maintained by the vibroelectronic energy. This energy exchange is implemented into the ODE loops. The electron velocity for the elastic energy exchange was also changed from the bulk electron velocity to the average thermal velocity from the Boltzmann distribution in statistical mechanics
            - Removed the reliance on the bulk electron velocity in the update_reaction_rates function due to the update above
            - Fixed some silly unit errors
            - Included a git diff file for Eilmer developers to better implement my model into the gdtk source code
            - Included a git log file to better inform which version of Eilmer my model was made in
            - Cleaned up unnecessary code to save on some computational cost