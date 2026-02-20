title: "Low Temperature Air Plasma README"

    Prerequisites:
        NOTE - The low temperature air plasma model was developed in Eilmer Version 4.1.0, if you wish to convert the files to Eilmer Version 5.0.0 (lmr),
            you will need to make the conversions yourself (sorry)
        It is presumed you have Eilmer installed on your device and all of its prerequisites, the below prerequisites are additional to the base Eilmer model:
            - Gas-dynamic library
                $ cd gdtk/src/gas 
                $ make install

    To run the unit tests:
        Copy the kinetics file into gdtk/src/kinetics
        Copy the gas file into gdtk/src/gas
        Copy the .lua files in this repository (in sample-data and sample-input) into their respective gas and kinetics sample-data and sample-input
            directories, the unit tests will refer to them
        Choose your case (relaxation or energy investigation)
        In your linux environment, move to the kinetics directory ($ cd gdtk/src/kinetics) and run all kinetics tests ($ make test)
            which will run through all of the kinetics file tests

    To use the plotting scripts:
        If you used the relaxation test case, run the low_temp_air_plasma_relaxation_plotting.py script
        If you used the energy investigation test case, run the low_temp_air_plasma_energy_plotting.py script
    
