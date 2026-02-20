title: "Low Temperature Air Plasma README"

    Prerequisites:
        NOTE - The low temperature air plasma model was developed in Eilmer Version 4.1.0, if you wish to convert the files to Eilmer Version 5.0.0 (lmr),
            you will need to make the conversions yourself (sorry)
        It is presumed you have Eilmer installed on your device and all of its prerequisites, the below prerequisites are additional to the base Eilmer model:
            - Gas-dynamic library
                ```
                $ cd gdtk/src/gas 
                $ make install
                ```
    
    Additions to the source code:
        In gdtk/src/eilmer/globalconfig.d 
            line 1778:
                import gas.low_two_temperature_reacting_air;
            line 1826:
                if (cast(TwoTemperatureReactingAir)gm) { multiT = true; }
        
        In gdtk/src/gas/gas-package-test.tcl 
            line 130-132:
                test two-temperature-reacting-air-test {Testing Tom Condon's two-T low temperature reacting air gas model} -body {
                    exec ./two_temperature_reacting_air_test
                } -result {} -returnCodes {0}
        
        In gdtk/src/gas/gas_files.mk 
            line 34:
                $(GAS_DIR)/low_two_temperature_reacting_air.d
        
        In gdtk/src/gas/init_gas_model.d 
            line 29:
                import gas.low_two_temperature_reacting_air: TwoTemperatureReactingAir;
            line 127-129:
                case "TwoTemperatureReactingAir":
                    gm = new TwoTemperatureReactingAir(L);
                    break;
        
        In gdtk/src/gas/makefile
            line 39:
                two_temperature_reacting_air_test
            line 205:
                - rm -f ./two-temperature-reacting-air-model.lua
            line 597-603:
                two_temperature_reacting_air_test: $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) \
                $(KINETICS_FILES) $(KINETICS_LUA_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
                $(LIBEQC) $(LIBLUA)
                $(DMD) $(OF)$@ $(DFLAGS) $(DVERSION)$@ \
                    $(GAS_FILES) $(GAS_LUA_FILES) $(EQC_SRC_FILES) $(KINETICS_FILES) $(KINETICS_LUA_FILES) \
                    $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
                    $(LIBEQC) $(LIBLUA) $(DLINKFLAGS)
        
        In gdtk/src/kinetics/init_thermochemical_reacctor.d 
            line 33:
                import gas.low_two_temperature_reacting_air;
            line 55:
                import kinetics.low_two_temperature_reacting_air_kinetics;
            line 105-107:
                if ((cast(TwoTemperatureReactingAir) gmodel) !is null) {
                    reactor = new LowTwoTemperatureAirKinetics(fileName1, gmodel);
                }
            
        In gdtk/src/kinetics/kinetics-package-test.tcl 
            line 35-37:
                test two-temperature-reacting-air-kinetics-test {Testing Tom Cond's two-temperature reacting air reaction mechanism.} -body {
                    exec ./two_temperature_reacting_air_kinetics_test
                } -result {} -returnCodes {0}
        
        In gdtk/src/kinetics/kinetics_files.mk 
            line 28:
                $(KINETICS_DIR)/low_two_temperature_reacting_air_kinetics.d
        
        In gdtk/src/kinetics/makefile
            line 27:
                two_temperature_reacting_air_kinetics_test
            line 82-86:
                - rm -f two_temperature_reacting_air_kinetics_test_results.data
                - rm -f two_temperature_reacting_air_kinetics_test_results_v0.data
                - rm -f two_temperature_reacting_air_kinetics_test_results_v1.data
                - rm -f two_temperature_reacting_air_kinetics_test_results_v2.data
                - rm -f two_temperature_reacting_air_kinetics_test_results_v3.data
            line 145-151:
                two_temperature_reacting_air_kinetics_test:  $(KINETICS_FILES) $(GAS_FILES) $(EQC_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
                $(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
                $(LIBEQC) $(LIBLUA)
                ldmd2 -of$@ -debug -g -dip1008 -version=$@ \
                    $(KINETICS_FILES) $(GAS_FILES) $(EQC_SRC_FILES) $(UTIL_FILES) $(NM_FILES) $(NTYPES_FILES) \
                    $(KINETICS_LUA_FILES) $(GAS_LUA_FILES) \
                    $(LIBEQC) $(LIBLUA) $(DLINKFLAGS)

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
    