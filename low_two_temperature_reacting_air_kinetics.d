/**
 * low_two_temperature_reacting_air_kinetics.d
 *
 * Low two-temperature reacting air based on the model of:
 * Bernard Parent and Mikhail N. Shneider and Sergey O. Macheret
 * "Detailed Modeling of Plasmas for Computational Aerodynamics"
 * AIAA Journal, Vol. 54, No. 3, 898 (2016); doi 10.2514/1.J054624
 *
 * Authors: Thomas Condon
 * Inspired By: Daniel Smith, Rory Kelly and Peter J.
 * Version: 14-January-2026: initial development
 */

module kinetics.low_two_temperature_reacting_air_kinetics;

import std.stdio : writeln, writefln;
import std.format;
import std.math;
import std.conv;
import ntypes.complex;
import nm.number;
import nm.bbla;

import gas;
import gas.low_two_temperature_reacting_air;
import gas.physical_constants;
import util.lua;
import util.lua_service;
import kinetics.thermochemical_reactor;


final class LowTwoTemperatureAirKinetics : ThermochemicalReactor {
    // Add public member variables which are relevant for the reaction rates, for users to populate
    public double E_mag_p = 0.0; // Electric Field (V/m)
    public double v_elec_p = 0.0; // Electron Velocity (m/s)
    public double B_mag_p = 0.0; // Magnetic Field (T)
    public double Qb_p = 0.0; // (W/m^3)

    this(string fname, GasModel gmodel)
    {
        super(gmodel); // Hang on to a reference to the gas model
        _reactingairModel = cast(TwoTemperatureReactingAir) gmodel;
        _mu.length = _gmodel.n_species;
        foreach (isp; 0 .. _gmodel.n_species) _mu[isp].length = _gmodel.n_species;
        foreach (isp; 0 .. _gmodel.n_species) {
            foreach (jsp; 0 .. _gmodel.n_species) {
                double M_isp = _gmodel.mol_masses[isp];
                double M_jsp = _gmodel.mol_masses[jsp];
                _mu[isp][jsp] = (M_isp*M_jsp)/(M_isp + M_jsp);
                _mu[isp][jsp] *= 1000.0; // convert kg/mole to g/mole
            }
        }
        // Map species to their respective index in Eilmer
        i_e   = gmodel.species_index("e-");
        i_O   = gmodel.species_index("O");
        i_Op  = gmodel.species_index("O+");
        i_O2  = gmodel.species_index("O2");
        i_O2p = gmodel.species_index("O2+");
        i_O2m = gmodel.species_index("O2-");
        i_N   = gmodel.species_index("N");
        i_Np  = gmodel.species_index("N+");
        i_N2  = gmodel.species_index("N2");
        i_N2p = gmodel.species_index("N2+");
        i_NO  = gmodel.species_index("NO");
        i_NOp = gmodel.species_index("NO+");

        // Following 2T Argon kinetics, We need to pick a number of pieces out of the gas-model file.
        // Although some already exist in the GasModel object, they are private. Also need to load the Lua config
        auto L = init_lua_State();
        doLuaFile(L, fname);
        lua_getglobal(L, "TwoTemperatureReactingAir");
        if (!lua_isnil(L, -1)) 
        {
            _ion_tol = getDoubleWithDefault(L, -1, "ion_tol", 0.0);
            _integration_method = getStringWithDefault(L, -1, "integration_method", "Forward_Euler");
            _T_min_for_reaction = getDoubleWithDefault(L, -1, "T_min_for_reaction", 200.0); // set to 200.0K from 3000.0K
            _Te_default = getDoubleWithDefault(L, -1, "Te_default", 10000.0); // User can set to 200.0K for low temperature applications
            _n_step_suggest = getIntWithDefault(L, -1, "n_step_suggest", 10);
            _max_iter_newton = getIntWithDefault(L, -1, "max_iter_newton", 30);
        }
        lua_close(L);
    }

    // It is expected that params[0] is the electric field magnitude E, params[1] is the electron velocity in the electric field v_elec
    // params[2] is the magnetic field magnitude B, and params[3] is the electron beam power per unit volume Qb
    // These positions are an interface requirement for this 'self-contained' hard-coded model for explicit uses
    // Further changes will be made to generalise this model such that Carolyn J. can use it easily (hopefully)
    // For generalisation, simply just set Qb = 0.0 in params and everything should still run smoothly ?

    // The above has been changed after running into issues with having to pass in the params values constantly into different functions

    @nogc
    override void opCall(ref GasState Q, double tInterval, ref double dtSuggest, ref number[maxParams] params) 
    {
        if (dtSuggest < 0.0) {
            // A negative value indicated that the flow solver did not have
            // a suggested time-step size
            dtSuggest = tInterval/_n_step_suggest;
        }
        // There are changes only if the gas is hot enough
        if (Q.T > _T_min_for_reaction) {

            // Pull in the values from the public 'this' function thing
            number E_mag = to!number(this.E_mag_p); // Electric Field (V/m)
            number v_elec = to!number(this.v_elec_p); // Electron Velocity (m/s)
            number B_mag = to!number(this.B_mag_p); // Magnetic Field (T)
            number Qb = to!number(this.Qb_p); // (W/m^3)

            // Determine the current number densities for each species (/m^3)
            number[12] numden;
            _gmodel.massf2numden(Q, numden);
            number n_e   = numden[i_e]; // number density of e-
            number n_O   = numden[i_O]; // number density of O
            number n_Op  = numden[i_Op]; // number density of O+
            number n_O2  = numden[i_O2]; // number density of O2
            number n_O2p = numden[i_O2p]; // number density of O2+
            number n_O2m = numden[i_O2m]; // number density of O2-
            number n_N   = numden[i_N]; // number density of N
            number n_Np  = numden[i_Np]; // number density of N+
            number n_N2  = numden[i_N2]; // number density of N2
            number n_N2p = numden[i_N2p]; // number density of N2+
            number n_NO  = numden[i_NO]; // number density of NO
            number n_NOp = numden[i_NOp]; // number density of NO+

            // Determine the total number density of the plasma (/m^3), which is the sum of the number density of each species
            // since this is a model of an isolated reactor
            number N_plasma = 0.0;
            foreach (n; 0 .. _gmodel.n_species) {
                N_plasma += numden[n];
            }

            u_total = Q.u + Q.u_modes[0];

            // Give the amount of energy, we have a limit to the number of ions that can be formed before driving the
            // energy of the heavy particles too low, as per the constructor we will put a temperature limit of the heavy particles
            _u_min_heavy_particles = 3.0 / 2.0 * R * _T_min_for_reaction;
            number u_available = u_total - _u_min_heavy_particles;
            // number alpha_max = fmin(1.0, u_available / (R * (Q.T_modes[0] + _theta_ion))); -> We have no _theta_ion here ?
            // n_e_max = alpha_max * N_plasma;

            // Retain a copy of the initial state, in case the integration is restarted with smaller time steps
            number initial_n_e = n_e;
            number initial_n_O = n_O;
            number initial_n_Op = n_Op;
            number initial_n_O2 = n_O2;
            number initial_n_O2p = n_O2p;
            number initial_n_O2m = n_O2m;
            number initial_n_N = n_N;
            number initial_n_Np = n_Np;
            number initial_n_N2 = n_N2;
            number initial_n_N2p = n_N2p;
            number initial_n_NO = n_NO;
            number initial_n_NOp = n_NOp;
            number initial_translational_energy = Q.u;
            number initial_vibrational_energy = Q.u_modes[0];
            number initial_translational_energy_T = Q.T;

            // Calculate necessary characteristics of the plasma for the rate coefficients
            // E* = |E + Ve x B|/N
            number E_star = (E_mag + (v_elec * B_mag)) / N_plasma; // Reduced effective electric field in electron reference frame (Vm^2)
            // Qb* = Qb/N
            number Qb_star = Qb / N_plasma; // Ratio between electron beam power per unit volume and total number density of plasma (W)

            // Start with the suggested time step size
            int NumberSteps = cast(int) fmax(floor(tInterval/dtSuggest), 1.0);
            number[13] S;
            int integration_attempt = 0;
            bool finished_integration = false;

            while (!finished_integration) {
                ++integration_attempt;
                try {
                    // Integrate across the interval with fixed time steps
                    _chem_dt = tInterval/NumberSteps;

                    // Pack the state vector, ready for the integrator
                    S[0] = initial_n_e;
                    S[1] = initial_n_O;
                    S[2] = initial_n_Op;
                    S[3] = initial_n_O2;
                    S[4] = initial_n_O2p;
                    S[5] = initial_n_O2m;
                    S[6] = initial_n_N;
                    S[7] = initial_n_Np;
                    S[8] = initial_n_N2;
                    S[9] = initial_n_N2p;
                    S[10] = initial_n_NO;
                    S[11] = initial_translational_energy;
                    S[12] = initial_vibrational_energy;

                    switch (_integration_method) {
                    case "Forward_Euler":
                        foreach (n; 0 .. NumberSteps) {
                            number[13] myF = F(S, Q, E_star, Qb_star, v_elec, N_plasma);
                            foreach (i; 0 .. 13) {
                                S[i] += _chem_dt * myF[i];
                            }
                            if (!state_vector_is_within_limits(S)) {
                                string msg = "State vector not within limits";
                                debug { msg ~= format("\n   y=[%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g]", S[0], S[1], S[2], S[3], S[4], S[5], S[6], S[7], S[8], S[9], S[10], S[11], S[12]); }
                                throw new ThermochemicalReactorUpdateException(msg);
                            }
                        } // end foreach n
                        break;
                    case "RK4": // Same as argon RK4 but with more care around the energy exchange
                        number [13] k1, k2_in, k2, k3_in, k3, k4_in, k4;
                        foreach (n; 0 .. NumberSteps) {
                            number[13] myF = F(S, Q, E_star, Qb_star, v_elec, N_plasma);
                            foreach (i; 0 .. 13) { k1[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 11) { k2_in[i] = fmax(S[i] + k1[i]/2.0, 0.0); }
                            k2_in[11] = S[11] + k1[11]/2.0;
                            k2_in[12] = S[12] + k1[12]/2.0;
                            //foreach (i; 0 .. 13) { k2_in[i] = S[i] + k1[i]/2.0; }
                            myF = F(k2_in, Q, E_star, Qb_star, v_elec, N_plasma);
                            foreach (i; 0 .. 13) { k2[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 11) { k3_in[i] = fmax(S[i] + k2[i]/2.0, 0.0); }
                            k3_in[11] = S[11] + k2[11]/2.0;
                            k3_in[12] = S[12] + k2[12]/2.0;
                            //foreach (i; 0 .. 13) { k3_in[i] = S[i] + k2[i]/2.0; }
                            myF = F(k3_in, Q, E_star, Qb_star, v_elec, N_plasma);
                            foreach (i; 0 .. 13) { k3[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 11) { k4_in[i] = fmax(S[i] + k3[i], 0.0); }
                            k4_in[11] = S[11] + k3[11];
                            k4_in[12] = S[12] + k3[12];
                            //foreach (i; 0 .. 13) { k4_in[i] = S[i] + k3[i]; }
                            myF = F(k4_in, Q, E_star, Qb_star, v_elec, N_plasma);
                            foreach (i; 0 .. 13) { k4[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 11) { S[i] += fmax(1.0/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]), 0.0); }
                            S[11] += 1.0/6.0*(k1[11]+2.0*k2[11]+2.0*k3[11]+k4[11]);
                            S[12] += 1.0/6.0*(k1[12]+2.0*k2[12]+2.0*k3[12]+k4[12]);
                            //foreach (i; 0 .. 13) { S[i] += 1.0/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]); }
                            if (!state_vector_is_within_limits(S)) {
                                string msg = "State vector not within limits";
                                debug { msg ~= format("\n    y=[%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g]", S[0], S[1], S[2], S[3], S[4], S[5], S[6], S[7], S[8], S[9], S[10], S[11], S[12]); }
                                throw new ThermochemicalReactorUpdateException(msg);
                            }
                        } // end foreach n
                        break;
                    default:
                        throw new Exception("Invalid ODE update selection.");
                    } // end switch _integration_method
                    // Unpack the results, clipping round-off errors
                    n_e = fmax(S[0], 0.0);
                    n_O   = fmax(S[1], 0.0);
                    n_Op  = fmax(S[2], 0.0);
                    n_O2  = fmax(S[3], 0.0);
                    n_O2p = fmax(S[4], 0.0);
                    n_O2m = fmax(S[5], 0.0);
                    n_N   = fmax(S[6], 0.0);
                    n_Np  = fmax(S[7], 0.0);
                    n_N2  = fmax(S[8], 0.0);
                    n_N2p = fmax(S[9], 0.0);
                    n_NO  = fmax(S[10], 0.0);
                    Q.u_modes[0] = fmin(S[12], u_total);
                    Q.u = fmax(u_total - Q.u_modes[0], 0.0); // Energy conservation

                    // Reconstruct the other parts of the flow state
                    // Utilise charge neutrality
                    //n_NOp = (-1.0*n_e) + n_Op + n_O2p - n_O2m + n_Np + n_N2p;
                    n_NOp = n_e + n_O2m - n_Op - n_O2p - n_Np - n_N2p;
                
                    // Now to update the mass fractions, the electron temperature, thermodyanmic behaviour of the gas, and the sound speed !
                    numden[i_e] = n_e; numden[i_O] = n_O; numden[i_Op] = n_Op; numden[i_O2] = n_O2;
                    numden[i_O2p] = n_O2p; numden[i_O2m] = n_O2m; numden[i_N] = n_N; numden[i_Np] = n_Np;
                    numden[i_N2] = n_N2; numden[i_N2p] = n_N2p; numden[i_NO] = n_NO; numden[i_NOp] = n_NOp;
                    // Let's just make sure the density is going to be correct, having some issues after going through unit tests
                    number rho_calc = 0.0;
                    foreach (i; 0 .. _gmodel.n_species) {
                        rho_calc += numden[i] * _gmodel.mol_masses[i] / AvN;
                    }
                    Q.rho = rho_calc;
                    _gmodel.numden2massf(numden, Q);
                    // Ignore the below comment, I was oblivious to the fact that in a 2T model, the vibroelectronic temperature is not just the electron temperature
                    //Q.T_modes[0] = calc_Te_HOP(Q, E_star); // We actually calculate electron temperature here in the kinetics file, as opposed to in the gas file like usual
                    _gmodel.update_thermo_from_rhou(Q);
                    _gmodel.update_sound_speed(Q);

                    // Let's ensure we can capture some of the transient effects that may occur during massive/violent reactions
                    //number frac_change = abs(Q.T - initial_translational_energy_T) / initial_translational_energy_T; // We don't want non-sensical instant temp differences
                    //if (frac_change > 0.1) {
                    //    // Scale the change back a bit !
                    //    number scale = 0.1 / frac_change;
                    //    Q.u = initial_translational_energy + (Q.u - initial_translational_energy) * scale;
                    //    Q.u_modes[0] = fmax(u_total - Q.u, 0.0);
                    //    // Run the updates again
                    //    _gmodel.update_thermo_from_rhou(Q);
                    //    _gmodel.update_sound_speed(Q);
                    //}

                    finished_integration = true;
                    //if (integration_attempt == 1) {
                    //    // Success on the first attempt so increase the time step a little for the next call
                    //    _chem_dt *= 1.01;
                    //}
                } catch(Exception e) {
                    if (integration_attempt < 10) {
                        // We will allow a retry with a smaller time step
                        NumberSteps *= 3;
                    } else {
                        // We have taken too many attempts and have continued to fail
                        string msg = "Thermochemical update fails after several attempts.";
                        debug { msg ~= text("\n Previous exception message: ", e.msg); }
                        throw new ThermochemicalReactorUpdateException(msg);
                    }
                }
            } // End while !finished_integration

            // Remember the size of the successful time step for the next call
            dtSuggest = _chem_dt;
        } // end if the neutrals temperature is less than the temperature required for reactions
    } // end opCall()

    @nogc
    override void eval_source_terms(GasModel gmodel, ref GasState Q, ref number[] source)
    {
        // Pull in the values from the public 'this' function thing
        number E_mag = to!number(this.E_mag_p); // Electric Field (V/m)
        number v_elec = to!number(this.v_elec_p); // Electron Velocity (m/s)
        number B_mag = to!number(this.B_mag_p); // Magnetic Field (T)
        number Qb = to!number(this.Qb_p); // (W/m^3)

        // Determine the current number densities for each species (/m^3)
        number[12] numden;
        _gmodel.massf2numden(Q, numden); // I feel like I am calling this excessively, hopefully it is not too computationally expensive :/
        number n_e   = numden[i_e]; // number density of e-
        number n_O   = numden[i_O]; // number density of O
        number n_Op  = numden[i_Op]; // number density of O+
        number n_O2  = numden[i_O2]; // numder density of O2
        number n_O2p = numden[i_O2p]; // number density of O2+
        number n_O2m = numden[i_O2m]; // number density of O2-
        number n_N   = numden[i_N]; // number density of N
        number n_Np  = numden[i_Np]; // number density of N+
        number n_N2  = numden[i_N2]; // number density of N2
        number n_N2p = numden[i_N2p]; // number density of N2+
        number n_NO  = numden[i_NO]; // number density of NO
        number n_NOp = numden[i_NOp]; // number density of NO+

        // Determine the total number density of the plasma (/m^3), which is the sum of the number density of each species
        // since this is a model of an isolated reactor
        number N_plasma = 0.0;
        foreach (n; 0 .. gmodel.n_species) {
            N_plasma += numden[n];
        }

        // Determine the reduced versions like before
        number E_star = (E_mag + (v_elec * B_mag)) / N_plasma;
        number Qb_star = Qb / N_plasma;

        // Determine the necessary rates to update the source terms
        number[13] myF = update_reaction_rates(Q, E_star, Qb_star, v_elec);
        number energy_leaving_trans = myF[11];
        number energy_entering_vib = myF[12];

        
        // Map the local enum indicies to Eilmer's dynamic indices and convert dn/dt
        // to the species production rates (kg/(s*m^3)) if the temperature is hot enough
        if (Q.T > _T_min_for_reaction) {
            number production_rate_sum = 0.0;
            source[i_e]   = myF[SRT_e]   * _gmodel.mol_masses[i_e]   / AvN;
            production_rate_sum += source[i_e];
            source[i_O]   = myF[SRT_O]   * _gmodel.mol_masses[i_O]   / AvN;
            production_rate_sum += source[i_O];
            source[i_Op]  = myF[SRT_Op]  * _gmodel.mol_masses[i_Op]  / AvN;
            production_rate_sum += source[i_Op];
            source[i_O2]  = myF[SRT_O2]  * _gmodel.mol_masses[i_O2]  / AvN;
            production_rate_sum += source[i_O2];
            source[i_O2p] = myF[SRT_O2p] * _gmodel.mol_masses[i_O2p] / AvN;
            production_rate_sum += source[i_O2p];
            source[i_O2m] = myF[SRT_O2m] * _gmodel.mol_masses[i_O2m] / AvN;
            production_rate_sum += source[i_O2m];
            source[i_N]   = myF[SRT_N]   * _gmodel.mol_masses[i_N]   / AvN;
            production_rate_sum += source[i_N];
            source[i_Np]  = myF[SRT_Np]  * _gmodel.mol_masses[i_Np]  / AvN;
            production_rate_sum += source[i_Np];
            source[i_N2]  = myF[SRT_N2]  * _gmodel.mol_masses[i_N2]  / AvN;
            production_rate_sum += source[i_N2];
            source[i_N2p] = myF[SRT_N2p] * _gmodel.mol_masses[i_N2p] / AvN;
            production_rate_sum += source[i_N2p];
            source[i_NO]  = myF[SRT_NO]  * _gmodel.mol_masses[i_NO]  / AvN;
            production_rate_sum += source[i_NO];
            source[i_NOp] = -1.0 * production_rate_sum; // The sum of all species production rates must be 0 due to mass conservation.
            source[12] = energy_leaving_trans;
            source[13] = energy_entering_vib;
        } else {
            source[i_e]   = 0.0;
            source[i_O]   = 0.0;
            source[i_Op]  = 0.0;
            source[i_O2]  = 0.0;
            source[i_O2p] = 0.0;
            source[i_O2m] = 0.0;
            source[i_N]   = 0.0;
            source[i_Np]  = 0.0;
            source[i_N2]  = 0.0;
            source[i_N2p] = 0.0;
            source[i_NO]  = 0.0;
            source[i_NOp] = 0.0;
            source[12] = 0.0;
            source[13] = 0.0;
        }
    }


    private:
    @nogc
    number[13] F(ref const(number[13]) S, ref GasState Q, number E_star, number Qb_star, number v_elec, number N_plasma) {
        update_Q_from_state_vector(S, Q, E_star, N_plasma);
        return update_reaction_rates(Q, E_star, Qb_star, v_elec);
    }

    @nogc
    void update_Q_from_state_vector(ref const(number[13]) S, ref GasState Q, number E_star, number N_plasma) {
        // Definition of the state vector is:
        // [0] -> [10] number density of e-, O, O+, ... , NO+
        // [11] Translational energy of heavy particles
        number[12] S_dash; S_dash[0] = 0.0; S_dash[1] = 0.0; S_dash[2] = 0.0;
        S_dash[3] = 0.0; S_dash[4] = 0.0; S_dash[5] = 0.0; S_dash[6] = 0.0;
        S_dash[7] = 0.0; S_dash[8] = 0.0; S_dash[9] = 0.0; S_dash[10] = 0.0;
        S_dash[11] = 0.0;

        // Update the GasState from the state vector, that way we always work the rate calculations
        // from a physically realisable state
        number n_e = S[0]; number n_O = S[1]; number n_Op = S[2]; number n_O2 = S[3];
        number n_O2p = S[4]; number n_O2m = S[5]; number n_N = S[6]; number n_Np = S[7];
        number n_N2 = S[8]; number n_N2p = S[9]; number n_NO = S[10];
        // Utilise charge neutrality
        number n_NOp = n_e + n_O2m - n_Op - n_O2p - n_Np - n_N2p;
        if (n_e < 0.0) {
            // Do not let the number of electrons go negative
            string msg = "Electron number density tried to go negative.";
            throw new Exception(msg);
        }
        // Clip off rounding errors
        //Q.u = fmin(S[11], u_total);
        //Q.u_modes[0] = fmax(u_total - Q.u, 0.0);
        Q.u_modes[0] = fmin(S[12], u_total);
        Q.u = fmax(u_total - Q.u_modes[0], 0.0);
        
        number Te = calc_Te_HOP(Q, E_star);
        number T = Q.T;

        // Now update the mass fractions, then update the GasState from rho and u
        number[12] nden;
        foreach (ref n; nden) n = to!number(0.0);
        nden[i_e] = n_e; nden[i_O] = n_O; nden[i_Op] = n_Op; nden[i_O2] = n_O2;
        nden[i_O2p] = n_O2p; nden[i_O2m] = n_O2m; nden[i_N] = n_N; nden[i_Np] = n_Np;
        nden[i_N2] = n_N2; nden[i_N2p] = n_N2p; nden[i_NO] = n_NO; nden[i_NOp] = n_NOp;
        // Let's just make sure the density is going to be correct, having some issues after going through unit tests
        number rho_calc = 0.0;
        foreach (i; 0 .. _gmodel.n_species) {
            rho_calc += nden[i] * _gmodel.mol_masses[i] / AvN;
        }
        Q.rho = rho_calc;
        _gmodel.numden2massf(nden, Q);
        _gmodel.update_thermo_from_rhou(Q);
        _gmodel.update_sound_speed(Q);
    }

    @nogc
    number calc_Te_HOP(in GasState Q, number E_star) {
        //====================================================
        // Electron temperature calculation taken from table 3
        // E* = |E| / N (V*m^2)
        // Te = max{T, exp[sum(k_n * (lnE*)^n)]}
        // Avoid log domain errors by fitting Te = T when E_star is too small
        if (E_star <= to!number(1.0e-25)) return Q.T;

        number lnEs = log(E_star);

        number poly = to!number(0.0);

        // Loop through table 3 coefficients with Horner's method (recommended for speed supposedly)
        foreach_reverse(coef; k_Te_coef) {
            poly = poly * lnEs +  to!number(coef);
        }

        number Te_calc = exp(poly);
        if (Te_calc < Q.T) return Q.T;
        return Te_calc;
    }

    @nogc
    number calculate_tau_vt(in GasState Q, int isp) {
        // Millikan-White formula (equation taken from Gnoffo et al. 1989)
        number sum = 0.0;
        number totalND = 0.0;
        // Gnoffo et al. Table 1
        double A_MW; // A parameter in Millikan-White expression
        if (isp == i_O2) {
            A_MW = 129.0; // O2
        } else if (isp == i_N2) {
            A_MW = 220.0; // N2
        } else if (isp == i_NO) {
            A_MW = 168.0; // NO
        } else {
            A_MW = 0.0;  // Should never happen, but safe fallback
        }

        foreach (csp; 0 .. _gmodel.n_species) {
            if (csp == i_e) continue; // Electrons don't contribute

            number nd = Q.massf[csp] * Q.rho / _gmodel.mol_masses[csp] * AvN;
            sum += nd * exp(A_MW * (pow(Q.T, -1.0/3.0) - 0.015*pow(_mu[isp][csp], 0.25)) - 18.42);
            totalND += nd;
        }

        number tau_MW = sum / totalND;  // p*tau in (atm*s)
        tau_MW *= 101325.0 / Q.p;       // Convert to seconds

        //  Park high-temp correction
        number ndP = Q.massf[isp] * Q.rho / _gmodel.mol_masses[isp] * AvN;
        number cBar = sqrt(8.0 * k_b * Q.T / (to!double(PI) * _gmodel.mol_masses[isp] / AvN));
        number sigma = 1.0e-20; // m^2
        number tau_P = 1.0 / (sigma * cBar * ndP);

        return tau_MW + tau_P;
    }

    @nogc
    number[13] update_reaction_rates(in GasState Q, number E_star, number Qb_star, number v_elec) {
        // Compute the rate of change of the state vector for the reactions
        // Definition of this new rate of change state vector is:
        // [0] -> [10] dn/dt of e-, O, O+, ... , NO (excluding NO+ since it is a dependent variable)
        // [11] energy leaving translational mode
        // [12] energy entering vibrational mode
        number[13] S_dash; S_dash[0] = 0.0; S_dash[1] = 0.0; S_dash[2] = 0.0;
        S_dash[3] = 0.0; S_dash[4] = 0.0; S_dash[5] = 0.0; S_dash[6] = 0.0;
        S_dash[7] = 0.0; S_dash[8] = 0.0; S_dash[9] = 0.0; S_dash[10] = 0.0;
        S_dash[11] = 0.0; S_dash[12] = 0.0;

        // Determine the current number densities for each species (/m^3)
        number[12] numden;
        _gmodel.massf2numden(Q, numden);
        number n_e   = numden[i_e]; // number density of e-
        number n_O   = numden[i_O]; // number density of O
        number n_Op  = numden[i_Op]; // number density of O+
        number n_O2  = numden[i_O2]; // numder density of O2
        number n_O2p = numden[i_O2p]; // number density of O2+
        number n_O2m = numden[i_O2m]; // number density of O2-
        number n_N   = numden[i_N]; // number density of N
        number n_Np  = numden[i_Np]; // number density of N+
        number n_N2  = numden[i_N2]; // number density of N2
        number n_N2p = numden[i_N2p]; // number density of N2+
        number n_NO  = numden[i_NO]; // number density of NO
        number n_NOp = numden[i_NOp]; // number density of NO+
        number n_ions = n_Op + n_O2p - n_O2m + n_Np + n_N2p + n_NOp; // HOW DOES THE NEGATIVE ION COME INTO PLAY FOR ENERGY EXCHANGE?? - Just gets substituted?

        // Determine the total number density of the plasma (/m^3), which is the sum of the number density of each species
        // since this is a model of an isolated reactor
        number N_plasma = 0.0;
        foreach (n; 0 .. _gmodel.n_species) {
            N_plasma += numden[n];
        }
        number Te = calc_Te_HOP(Q, E_star);
        number T = Q.T;

        double lnEs = (E_star.re > 1e-30) ? log(E_star.re): -69.0;

        //=========================================================================================
        // Rate coefficients for determining the mass production per unit volume of the kth species
        // Table 1 from Parent, 2016 - converted to SI units for m^3/s

        // ========== TOWNSEND IONIZATION ==========

        // 1a: e- + N2 -> N2+ + e- + e-
        number k1a = exp(-0.0105809 * lnEs*lnEs - 2.40411e-75 * pow(lnEs, 46)); // cm3/s
        k1a *= 1.0e-6;

        // 1b: e- + O2 -> O2+ + e- + e-
        number k1b = exp(-0.0102785 *lnEs*lnEs - 2.42260e-75 * pow(lnEs, 46)); // cm3/s
        k1b *= 1.0e-6;

        // 1c: e- + NO -> NO+ + e- + e-
        // Parent, 2021 provided this additional ionization reaction for NO
        number k1c = exp(-5.9890e-6 * pow(lnEs, 4) + 2.5988e-84 * pow(lnEs, 51)); // cm3/s
        k1c *= 1.0e-6;

        // 1d: e- + O -> O+ + e- + e-
        // Parent, 2021 provided this additional ionization reaction for O
        number k1d = (3.6e31/AvN) * pow(Te, -2.91) * exp(-158000.0/Te); // cm3/s
        k1d *= 1.0e-6;

        //1e: e- + N -> N+ + e- + e-
        // Parent, 2021 provided this additional ionization reaction for N
        number k1e = (1.1e32/AvN) * pow(Te, -3.14) * exp(-169000.0/Te); // cm3/s
        k1e *= 1.0e-6;

        {
        // Need to look for the below rate coefficients ? :
        // e- + N2 -> N + N + e-
        // No this is the same as regular N2 dissociation anyway
        }

        // ========== ELECTRON-ION RECOMBINATION ==========

        // 2a: e- + O2+ -> O + O
        // Parent, 2024 updated to the below equation from 2.0e-7*pow(300.0/Te, 0.7)
        number k2a = 2.4e-7 * pow(300.0/Te, 0.7); // cm3/s
        k2a *= 1.0e-6;

        // 2b: e- + N2+ -> N + N
        number k2b = 2.8e-7 * sqrt(300.0/Te); // cm3/s
        k2b*= 1.0e-6;

        // 2c: e- + NO+ -> N + O
        // Parent, 2024 provided this additional reaction
        number k2c = 3.0e-7 * pow(300.0/Te, 0.56); // cm3/s
        k2c *= 1.0e-6;

        // 2d: e- + e- + N+ -> e- + N (3-Body)
        // Parent, 2024 provided this additional reaction
        number k2d = 2.2e40 * pow(AvN, -2) * pow(Te, -4.5); // cm6/s
        k2d *= 1.0e-12;

        // 2e: e- + e- + O+ -> e- + O (3-Body)
        // Parent, 2024 provided this additional reaction
        number k2e = 2.2e40 * pow(AvN, -2) * pow(Te, -4.5); // cm6/s
        k2e *= 1.0e-12;

        // ========== ION-ION RECOMBINATION ==========

        // 3a: O2- + N2+ -> O2 + N2
        number k3a = 2.0e-7 * sqrt(300.0/T); // cm3/s
        k3a *= 1.0e-6;

        // 3b: O2- + O2+ -> O2 + O2
        number k3b = 2.0e-7 * sqrt(300.0/T); // cm3/s
        k3b *= 1.0e-6;

        // 4a: O2- + N2+ + N2 -> O2 + N2 + N2 (3-Body)
        number k4a = 2.0e-25 * pow(300.0/T, 2.5); // cm6/s
        k4a *= 1.0e-12;

        // 4b: O2- + O2+ + N2 -> O2 + O2 + N2 (3-Body)
        number k4b = 2.0e-25 * pow(300.00/T, 2.5); // cm6/s
        k4b *= 1.0e-12;

        // 4c: O2- + N2+ + O2 -> O2 + N2 + O2 (3-Body)
        number k4c = 2.0e-25 *  pow(300.0/T, 2.5); //cm6/s
        k4c *= 1.0e-12;

        // 4d: O2- + O2+ + O2 -> O2 + O2 + O2 (3-Body)
        number k4d = 2.0e-25 * pow(300.0/T, 2.5); // cm6/s
        k4d *= 1.0e-12;

        // ========== ELECTRON ATTACHMENT ==========

        // 5a: e- + O2 + O2 -> O2- + O2
        number k5a = 1.4e-29 * (300.0/Te) * exp(-600.0/T) * exp(700.0 * (Te - T)/(Te*T)); // cm6/s
        k5a *= 1.0e-12;

        // 5b: e- + O2 + N2 -> O2- + N2
        number k5b = 1.07e-31 * pow(300.0/Te,2) * exp(-70.0/T) * exp(1500.0 * (Te - T)/(Te*T)); // cm6/s
        k5b *= 1.0e-12;

        // 6: O2- + O2 -> e- + O2 + O2
        number k6 = 8.6e-10 * exp(-6030.0/T) * (1 - exp(-1570.0/T)); // cm6/s
        k6 *= 1.0e-12;

        // ========== ELECTRON BEAM IONIZATION ==========

        // 7a: O2 -> e- + O2+
        number k7a = 2.0e17 * Qb_star; // /s
                
        // 7b: N2 -> e- + N2+
        number k7b = 1.8e17 * Qb_star; // /s
        
        // The following reaction possibilities 8-10 will not provide much contribution
        // to the concentration of species in low temperature plasma due to their reliance
        // on the gas temperature T. Many of these reactions were exempt from Parent, 2016 
        // because of this. The ones which weren't included are labeled as being taken from
        // a Parent, 2024 paper which expected hotter gas temperatures than this model

        // ========== OXYGEN DISSOCIATION ==========

        // 8a: O2 + O2 -> 2O + O2
        number k8a = 3.7e-8 * exp(-59380.0/T) * (1 - exp(-2240.0/T)); // cm3/s
        k8a *= 1.0e-6;

        // 8b: O2 + N2 -> 2O + N2
        number k8b = 9.3e-9 * exp(-59380.0/T) * (1 - exp(-2240.0/T)); // cm3/s
        k8b *= 1.0e-6;

        // 8c: O2 + O -> 3O
        number k8c = 1.3e-7 * exp(-59380.0/T) * (1 - exp(-2240.0/T)); // cm3/s
        k8c *= 1.0e-6;

        // 8d: O2 + N → O + O + N
        // Parent, 2024 provided this additional reaction
        number k8d = 3.0e18 / (AvN*T) * exp(-59400.0/T); // cm3/s
        k8d *= 1.0e-6;

        // 8e: O2 + NO → O + O + NO
        // Parent, 2024 provided this additional reaction
        number k8e = 3.0e18 / (AvN*T) * exp(-59400.0/T); // cm3/s
        k8e *= 1.0e-6;

        // 8f: O2 + NO+ → O + O + NO+
        // Parent, 2024 provided this additional reaction
        number k8f = 3.0e18 / (AvN*T) * exp(-59400.0/T); // cm3/s
        k8f *= 1.0e-6;

        // ========== NITROGEN DISSOCIATION ==========

        // 8g: N2 + O2 -> 2N + O2
        number k8g = 5.0e-8 * exp(-113200.0/T) * (1 - exp(-3354.0/T)); // cm3/s
        k8g *= 1.0e-6;

        // 8h: N2 + N2 -> 2N + N2
        number k8h = 5.0e-8 * exp(-113200.0/T) * (1 - exp(-3354.0/T)); // cm3/s
        k8h *= 1.0e-6;

        // 8i: N2 + O -> 2N + O
        number k8i = 1.1e-7 * exp(-113200.0/T) * (1 - exp(-3354/T)); // cm3/s
        k8i *= 1.0e-6;

        // 8j: N2 + N -> 3N
        // Parent, 2024 provided this additional reaction
        number k8j = 1.3e20 / (AvN*T) * exp(-113200.0/T); // cm3/s
        k8j *= 1.0e-6;

        // 8k: N2 + NO → N + N + NO
        // Parent, 2024 provided this additional reaction
        number k8k = 1.9e19 / (AvN*T) * exp(-113200.0/T); // cm3/s
        k8k *= 1.0e-6;

        // 8l: N2 + NO+ → N + N + NO+
        // Parent, 2024 provided this additional reaction
        number k8l = 1.9e19 / (AvN*T) * exp(-113200.0/T); // cm3/s
        k8l *= 1.0e-6;

        // ========== NITRIC OXIDE DISSOCIATION ==========

        // 8m: NO + O → N + O + O
        // Parent, 2024 provided this additional reaction
        number k8m = 2.4e17 / AvN / sqrt(T) * exp(-75500.0/T); // cm3/s
        k8m *= 1.0e-6;

        // 8n: NO + O2 → N + O + O2
        // Parent, 2024 provided this additional reaction
        number k8n = 2.4e17 / AvN / sqrt(T) * exp(-75500.0/T); // cm3/s
        k8n *= 1.0e-6;

        // 8o: NO + N → N + O + N
        // Parent, 2024 provided this additional reaction
        number k8o = 2.4e17 / AvN / sqrt(T) * exp(-75500.0/T); // cm3/s
        k8o *= 1.0e-6;

        // 8p: NO + N2 → N + O + N2
        // Parent, 2024 provice this additional reaction
        number k8p = 2.4e17 / AvN / sqrt(T) * exp(-75500.0/T); // cm3/s
        k8p *= 1.0e-6;

        // 8q: NO + NO+ → N + O + NO+
        // Parent, 2024 provided this additional reaction
        number k8q = 2.4e17 / AvN / sqrt(T) * exp(-75500.0/T); // cm3/s
        k8q *= 1.0e-6;

        // ========== NITROGEN AND OXYGEN EXCHANGE ==========

        // 8r: O + N2 -> N + NO
        // Parent, 2024 provided this additional reaction
        number k8r = 6.8e13 / AvN * exp(-37750.0/T); // cm3/s
        k8r *= 1.0e-6;

        // 8s: O + NO -> N + O2
        // Parent, 2024 provided this additional reaction
        number k8s = 4.3e7 / AvN * pow(T, 1.5) * exp(-19100.0/T); //cm3/s
        k8s *= 1.0e-6;

        // 8t: N + NO -> O + N2
        // Parent, 2024 provided this additional reaction
        number k8t = 1.5e13 / AvN; // cm3/s
        k8t *= 1.0e-6;

        // 8u: N + O2 -> O + NO
        //Parent, 2024 provided this additional reaction
        number k8u = 1.8e8 / AvN * pow(T, 1.5) * exp(-3020.0/T); // cm3/s
        k8u *= 1.0e-6;

        // ========== BACKWARD OXYGEN DISSOCIATION ==========

        // 9a: O + O + O2 -> 2O2
        number k9a = 2.45e-31 * pow(T, -0.63); // cm6/s
        k9a *= 1.0e-12;

        // 9b: O + O + N2 -> O2 + N2
        number k9b = 2.76e-34 * exp(720.0/T); // cm6/s
        k9b *= 1.0e-12;

        // 9c: O + O + O -> O2 + O
        number k9c = 8.8e-31 * pow(T, -0.63); // cm6/s
        k9c *= 1.0e-12;

        // ========== BACKWARD NITROGEN DISSOCIATION ==========

        // 9d: N + N + O2 -> N2 + O2
        number k9d = 8.27e-34 * exp(500.0/T); // cm6/s
        k9d *= 1.0e-12;

        // 9e: N + N + N2 -> 2N2
        number k9e = 8.27e-34 * exp(500.0/T); // cm6/s
        k9e *= 1.0e-12;

        // 9f: N + N + O -> N2 + O
        number k9f = 8.27e-34 * exp(500.0/T); // cm6/s
        k9f *= 1.0e-12;

        // 9g: N + N + N -> N2 + N
        number k9g = 8.27e-34 * exp(500.0/T); // cm6/s
        k9g *= 1.0e-12;

        // ========== ASSOCIATIVE IONIZATION ==========

        // 10: N + O -> NO+ + e-
        // Parent, 2024 provided this additional reaction
        number k10 = 1.3e8 / AvN * T * exp(-31900.0/T); // cm3/s
        k10 *= 1.0e-6;

        //============================================================================================
        // Determining the rates of change of the reactants based on their rate coefficients over time
        // Requires the stoichiometric ratios (SR) from the SRT
        // dn_s/dt = sum_r(SR_s * k_r * product_j(n_j ^alpha_j,r))
        // First, determine the rates (R1a = k1a * n_e * n_N2)

        number[Num_Reactions] reaction_rates;

        reaction_rates[0] = k1a * n_e * n_N2;

        reaction_rates[1] = k1b * n_e * n_O2;

        reaction_rates[2] = k1c * n_e * n_NO;

        reaction_rates[3] = k1d * n_e * n_O;

        reaction_rates[4] = k1e * n_e * n_N;

        reaction_rates[5] = k2a * n_e * n_O2p;

        reaction_rates[6] = k2b * n_e * n_N2p;

        reaction_rates[7] = k2c * n_e * n_NOp;

        reaction_rates[8] = k2d * n_e * n_e * n_Np;

        reaction_rates[9] = k2e * n_e * n_e * n_Op;

        reaction_rates[10] = k3a * n_O2m * n_N2p;

        reaction_rates[11] = k3b * n_O2m * n_O2p;

        reaction_rates[12] = k4a * n_O2m * n_N2p * n_N2;
        
        reaction_rates[13] = k4b * n_O2m * n_O2p * n_N2;

        reaction_rates[14] = k4c * n_O2m * n_N2p * n_O2;

        reaction_rates[15] = k4d * n_O2m * n_O2p * n_O2;

        reaction_rates[16] = k5a * n_e * n_O2 * n_O2;

        reaction_rates[17] = k5b * n_e * n_O2 * n_N2;

        reaction_rates[18] = k6 * n_O2m * n_O2;

        reaction_rates[19] = k7a * n_O2;

        reaction_rates[20] = k7b * n_N2;

        reaction_rates[21] = k8a * n_O2 * n_O2;

        reaction_rates[22] = k8b * n_O2 * n_N2;

        reaction_rates[23] = k8c * n_O2 * n_O;

        reaction_rates[24] = k8d * n_O2 * n_N;

        reaction_rates[25] = k8e * n_O2 * n_NO;

        reaction_rates[26] = k8f * n_O2 * n_NOp;

        reaction_rates[27] = k8g * n_N2 * n_O2;

        reaction_rates[28] = k8h * n_N2 * n_N2;

        reaction_rates[29] = k8i * n_N2 * n_O;

        reaction_rates[30] = k8j * n_N2 * n_N;

        reaction_rates[31] = k8k * n_N2 * n_NO;

        reaction_rates[32] = k8l * n_N2 * n_NOp;

        reaction_rates[33] = k8m * n_NO * n_O;

        reaction_rates[34] = k8n * n_NO * n_O2;

        reaction_rates[35] = k8o * n_NO * n_N;

        reaction_rates[36] = k8p * n_NO * n_N2;

        reaction_rates[37] = k8q * n_NO * n_NOp;

        reaction_rates[38] = k8r * n_O * n_N2;

        reaction_rates[39] = k8s * n_O * n_NO;

        reaction_rates[40] = k8t * n_N * n_NO;

        reaction_rates[41] = k8u * n_N * n_O2;

        reaction_rates[42] = k9a * n_O * n_O * n_O2;

        reaction_rates[43] = k9b * n_O * n_O * n_N2;

        reaction_rates[44] = k9c * n_O * n_O * n_O;

        reaction_rates[45] = k9d * n_N * n_N * n_O2;

        reaction_rates[46] = k9e * n_N * n_N * n_N2;

        reaction_rates[47] = k9f * n_N * n_N * n_O;

        reaction_rates[48] = k9g * n_N * n_N * n_N;

        reaction_rates[49] = k10 * n_N * n_O;

        // Loop over each reaction
        foreach (r; 0 .. Num_Reactions) {
            number rate = reaction_rates[r];

            // Loop over each species
            foreach (s; 0 .. Num_Species - 1) {
                // Don't multiply if the SR is 0
                number coeff = SRT[r][s];
                if (coeff != 0.0) {
                    // Store the time evolution in number density of the plasma species in S_dash (as defined before)
                    S_dash[s] += (rate * coeff);
                }
            }
        }
        // Energy moving from '-1' to the electronic modes due to reactions and collisions
        // '-1' Translational/Rotational
        // '0' Vibrational (electronic)
        // '1' Free Electron (electronic)

        number Q_e_O, Q_e_N, Q_e_O2, Q_e_N2, Q_e_NO;

        // Calculate collision cross sectional areas for each of the neutrals
        // Based on Gnoffo(1989) curve fits of electron-neutral energy exchange
        Q_e_O = 1.2e-20 + 1.7e-24*Te - 2.0e-28*Te^^2; // Higher order constant defined as -2.0e-29 in Gnoffo, but Roshan changed to -2.0e-28 ?
        Q_e_O2 = 2.0e-20 + 6.0e-24*Te;
        Q_e_N = 5.0e-20;
        Q_e_N2 = 7.5e-20 + 5.5e-24*Te - 1.0e-28*Te^^2;
        Q_e_NO = 1.0e-19;

        // Calculate the collisions of each neutral with electrons, then sum them up for total electron-neutral collision, based on Imamura (2018)
        number nu_e_O = 4.0/3.0 * Q_e_O * n_O * v_elec;
        number nu_e_O2 = 4.0/3.0 * Q_e_O2 * n_O2 * v_elec;
        number nu_e_N = 4.0/3.0 * Q_e_N * n_N * v_elec;
        number nu_e_N2 = 4.0/3.0 * Q_e_N2 * n_N2 * v_elec;
        number nu_e_NO = 4.0/3.0 * Q_e_NO * n_NO * v_elec;
        number nu_e_n = nu_e_O + nu_e_O2 + nu_e_N + nu_e_N2 + nu_e_NO;

        // Calculate the collision of each ion with electrons, based on Imamura (2018)
        number  nu_e_i = (6.0 * to!number(PI)) * ((q_e^^2) / (12.0 * to!number(PI) * eps_zero * k_b * Te))^^2 * log(12.0 * to!number(PI) * pow(eps_zero * k_b / (q_e^^2), 3.0/2.0) * sqrt(pow(Te, 3.0) / n_e)) * n_ions * v_elec;
        
        // Find total collision frequency
        number nu_tot = nu_e_n + nu_e_i; // Not useful anymore :(

        // Accumulate the elastic energy exchange contributions from neutrals and ions
        number elastic_energy_constants = 3.0 * n_e * _m_e * k_b * (Te - T); // Terms not related to the species
        // The neutral species have corresponding collisions, so we must manually tally up the elastic energy for the neutrals
        number een_neutrals = (nu_e_O / (_gmodel.mol_masses[i_O]/AvN)) + (nu_e_O2 / (_gmodel.mol_masses[i_O2]/AvN)) + (nu_e_N / (_gmodel.mol_masses[i_N]/AvN)) + (nu_e_N2 / (_gmodel.mol_masses[i_N2]/AvN)) + (nu_e_NO / (_gmodel.mol_masses[i_NO]/AvN));
        // The ionic species have one overarching collision, so we can just divide by the average mass of the ions instead of tallying
        number ions_avg_mass = ((_gmodel.mol_masses[i_Op] + _gmodel.mol_masses[i_O2p] + _gmodel.mol_masses[i_O2m] + _gmodel.mol_masses[i_Np] + _gmodel.mol_masses[i_N2p] + _gmodel.mol_masses[i_NOp]) / AvN) / 6.0;
        number een_ions = nu_e_i / ions_avg_mass;
        number elastic_energy = elastic_energy_constants * (een_neutrals + een_ions); // W/m3

        // Accumulate the inelastic energy exchange contributions from the electrons
        // Q_inel = sumof_l(n_e * k_l * N_l * del_eps_l), k_l is the reaction rate coefficient for the lth electron impact process, N_l
        // is the density of the gas phase collision partner, del_eps_l is the corresponding change in the electron energy,
        // positive for gaining energy and negative for energy losses. sumof_l(n_e*k_l*del_eps_l) is equivalent to what we did before. Hence,
        //number n_neutrals = n_O + n_O2 + n_N + n_N2 + n_NO;
        //number inelastic_energy = S_dash[SRT_e] * n_neutrals; REMOVE FOR NOW, HAVING ISSUES

        // Calculate and store this elastic energy exchange rate in S_dash. We add inelastic energy rather than minus, as the change in electron
        // energy will take care of the signs for us :)
        number energy_exchange = elastic_energy;// + inelastic_energy;

        // Now to determine V-T relaxation
        number VT_relax_rate = 0.0;
        int[3] vib_species = [i_O2, i_N2, i_NO]; 
        number[3] theta_v_array = [theta_v_O2, theta_v_N2, theta_v_NO];

        foreach (idx; 0 .. 3) { // Only want O2, N2, and NO
            int isp = vib_species[idx];
            number theta_v = theta_v_array[idx];

            // Calculate relaxation time based on Millikan-White and Park high-T correction
            number tau_vt = calculate_tau_vt(Q, isp);

            // Equilibrium vib energy from transrotational temperature
            number e_v_star = (R_universal / _gmodel.mol_masses[isp]) * theta_v / (exp(theta_v / Q.T) - 1.0);

            // Current vib energy from vibroelectronic temperature
            number e_v = (R_universal / _gmodel.mol_masses[isp]) * theta_v / (exp(theta_v / Q.T_modes[0]) - 1.0);

            // Landau-Teller relaxation rate
            VT_relax_rate += Q.massf[isp] * (e_v_star - e_v) / tau_vt;
        }

        // Apply this energy exchange to the S_dash to feed into the source terms
        //S_dash[11] = energy_exchange - VT_relax_rate; // Energy leaving translation
        //S_dash[12]  = VT_relax_rate; // Energy entering vibration
        // Trying new attempt at thingsssssssss
        S_dash[11] = (-1.0 * Q.rho * VT_relax_rate) + elastic_energy;
        S_dash[12] = (Q.rho * VT_relax_rate);
        return S_dash;
    } // end rates

    @nogc
    bool state_vector_is_within_limits(ref number[13] S)
    {
        bool result = true;
        if (S[0] < 0.0) { result = false; }
        if (S[0] > _n_e_max) { result = false; }
        if (S[11] < _u_min_heavy_particles) { result = false; }
        // Allow for finite-difference perturbations when evaluating Jacobian.
        if (S[11] > u_total*1.01) { result = false; }
        return result;
    }
private:
    // Table 3: Polynomial coefficients needed to determine the electron temperature
    immutable double[] k_Te_coef = [
        -3.69167532692495882511e+08, // k0
        -6.26956713747712671757e+07, // k1
        -4.65528490607805550098e+06, // k2
        -1.97394448288739687996e+05, // k3
        -5.22784662897089219769e+03, // k4
        -8.85545617874565635930e+01, // k5
        -9.36914737923363882821e-01, // k6
        -5.66073394421067171284e-03, // k7
        -1.49535882691330832494e-05  // k8
    ];

    // Stoichiometric ratio table (SRT)
    // Reactions are represented as the rows, whilst the species are the columns
    enum {SRT_e, SRT_O, SRT_Op, SRT_O2, SRT_O2p, SRT_O2m, SRT_N, SRT_Np, SRT_N2, SRT_N2p, SRT_NO, SRT_NOp}
    enum Num_Species = 12;
    enum Num_Reactions = 50;

    // Matrix [Reaction][Species]
    static immutable double[Num_Species][Num_Reactions] SRT = [
        //          e-    O     O+    O2    O2+    O2-    N     N+    N2    N2+    NO    NO+
        /* 1a */ [ 1.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  -1.0,  1.0,  0.0,  0.0],
        /* 1b */ [ 1.0,  0.0,  0.0, -1.0,   1.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 1c */ [ 1.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0, -1.0,  1.0],
        /* 1d */ [ 1.0, -1.0,  1.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 1e */ [ 1.0,  0.0,  0.0,  0.0,   0.0,   0.0, -1.0,   1.0,   0.0,  0.0,  0.0,  0.0],
        /* 2a */ [-1.0,  2.0,  0.0,  0.0,  -1.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 2b */ [-1.0,  0.0,  0.0,  0.0,   0.0,   0.0,  2.0,   0.0,   0.0, -1.0,  0.0,  0.0],
        /* 2c */ [-1.0,  1.0,  0.0,  0.0,   0.0,   0.0,  1.0,   0.0,   0.0,  0.0,  0.0, -1.0],
        /* 2d */ [-1.0,  0.0,  0.0,  0.0,   0.0,   0.0,  1.0,  -1.0,   0.0,  0.0,  0.0,  0.0],
        /* 2e */ [-1.0,  1.0, -1.0,  0.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 3a */ [ 0.0,  0.0,  0.0,  1.0,   0.0,  -1.0,  0.0,   0.0,   1.0, -1.0,  0.0,  0.0],
        /* 3b */ [ 0.0,  0.0,  0.0,  2.0,  -1.0,  -1.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 4a */ [ 0.0,  0.0,  0.0,  1.0,   0.0,  -1.0,  0.0,   0.0,   1.0, -1.0,  0.0,  0.0],
        /* 4b */ [ 0.0,  0.0,  0.0,  2.0,  -1.0,  -1.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 4c */ [ 0.0,  0.0,  0.0,  1.0,   0.0,  -1.0,  0.0,   0.0,   1.0, -1.0,  0.0,  0.0],
        /* 4d */ [ 0.0,  0.0,  0.0,  2.0,  -1.0,  -1.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 5a */ [-1.0,  0.0,  0.0, -1.0,   0.0,   1.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 5b */ [-1.0,  0.0,  0.0, -1.0,   0.0,   1.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 6  */ [ 1.0,  0.0,  0.0,  1.0,   0.0,  -1.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 7a */ [ 1.0,  0.0,  0.0, -1.0,   1.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],         
        /* 7b */ [ 1.0,  0.0,  0.0,  0.0,   0.0,   0.0,  0.0,   0.0,  -1.0,  1.0,  0.0,  0.0],
        /* 8a */ [ 0.0,  2.0,  0.0, -1.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 8b */ [ 0.0,  2.0,  0.0, -1.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 8c */ [ 0.0,  2.0,  0.0, -1.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 8d */ [ 0.0,  2.0,  0.0, -1.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 8e */ [ 0.0,  2.0,  0.0, -1.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 8f */ [ 0.0,  2.0,  0.0, -1.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 8g */ [ 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  2.0,   0.0,  -1.0,  0.0,  0.0,  0.0],
        /* 8h */ [ 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  2.0,   0.0,  -1.0,  0.0,  0.0,  0.0],
        /* 8i */ [ 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  2.0,   0.0,  -1.0,  0.0,  0.0,  0.0],
        /* 8j */ [ 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  2.0,   0.0,  -1.0,  0.0,  0.0,  0.0],
        /* 8k */ [ 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  2.0,   0.0,  -1.0,  0.0,  0.0,  0.0],
        /* 8l */ [ 0.0,  0.0,  0.0,  0.0,   0.0,   0.0,  2.0,   0.0,  -1.0,  0.0,  0.0,  0.0],
        /* 8m */ [ 0.0,  1.0,  0.0,  0.0,   0.0,   0.0,  1.0,   0.0,   0.0,  0.0, -1.0,  0.0],
        /* 8n */ [ 0.0,  1.0,  0.0,  0.0,   0.0,   0.0,  1.0,   0.0,   0.0,  0.0, -1.0,  0.0],
        /* 8o */ [ 0.0,  1.0,  0.0,  0.0,   0.0,   0.0,  1.0,   0.0,   0.0,  0.0, -1.0,  0.0],
        /* 8p */ [ 0.0,  1.0,  0.0,  0.0,   0.0,   0.0,  1.0,   0.0,   0.0,  0.0, -1.0,  0.0],
        /* 8q */ [ 0.0,  1.0,  0.0,  0.0,   0.0,   0.0,  1.0,   0.0,   0.0,  0.0, -1.0,  0.0],
        /* 8r */ [ 0.0, -1.0,  0.0,  0.0,   0.0,   0.0,  1.0,   0.0,  -1.0,  0.0,  1.0,  0.0],
        /* 8s */ [ 0.0, -1.0,  0.0,  1.0,   0.0,   0.0,  1.0,   0.0,   0.0,  0.0, -1.0,  0.0],
        /* 8t */ [ 0.0,  1.0,  0.0,  0.0,   0.0,   0.0, -1.0,   0.0,   1.0,  0.0, -1.0,  0.0],
        /* 8u */ [ 0.0,  1.0,  0.0, -1.0,   0.0,   0.0, -1.0,   0.0,   0.0,  0.0,  1.0,  0.0],
        /* 9a */ [ 0.0, -2.0,  0.0,  1.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 9b */ [ 0.0, -2.0,  0.0,  1.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 9c */ [ 0.0, -2.0,  0.0,  1.0,   0.0,   0.0,  0.0,   0.0,   0.0,  0.0,  0.0,  0.0],
        /* 9d */ [ 0.0,  0.0,  0.0,  0.0,   0.0,   0.0, -2.0,   0.0,   1.0,  0.0,  0.0,  0.0],
        /* 9e */ [ 0.0,  0.0,  0.0,  0.0,   0.0,   0.0, -2.0,   0.0,   1.0,  0.0,  0.0,  0.0],
        /* 9f */ [ 0.0,  0.0,  0.0,  0.0,   0.0,   0.0, -2.0,   0.0,   1.0,  0.0,  0.0,  0.0],
        /* 9g */ [ 0.0,  0.0,  0.0,  0.0,   0.0,   0.0, -2.0,   0.0,   1.0,  0.0,  0.0,  0.0],
        /* 10 */ [ 1.0, -1.0,  0.0,  0.0,   0.0,   0.0, -1.0,   0.0,   0.0,  0.0,  0.0,  1.0]
    ];

    immutable double AvN = Avogadro_number; // Avogadro's number, in /mole
    immutable double R = R_universal; // Universal gas constant, in Joules/kg/mole/Kelvin
    immutable double q_e = electron_volt_energy; // Charge of an electron, in Joules
    immutable double k_b = Boltzmann_constant; // Boltzmann's constant, in Joules/Kelvin
    immutable double eps_zero = vacuum_permittivity; // Permittivity of free space
    immutable double _m_e = 9.10938e-31; // mass of electron (kg)

    // Adjustable parameters, will be set in the constructor.
    int i_e, i_O, i_Op, i_O2, i_O2p, i_O2m, i_N, i_Np, i_N2, i_N2p, i_NO, i_NOp; // Species variables
    double _relTol;
    double _absTol;
    string _integration_method;
    double _T_min_for_reaction; // degrees K
    double _Te_default;
    double _ion_tol;
    int _n_step_suggest; // number of substeps to take if dtSuggest<0.0
    int _max_iter_newton;
    // The below are needed for VT relaxation
    double theta_v_N2 = 3353; // K
    double theta_v_NO = 2719; // K
    double theta_v_O2 = 2240; // K
    double[][] _mu;
    TwoTemperatureReactingAir _reactingairModel;

    // The following items are constants for the duration
    // of the update for our abstract isolated reactor.
    number N_plasma; // number density of atoms and ions combined
    number u_total; // energy within the reactor
    // Since the energy in the reactor is limited,
    // the number of ionized particles also limited.
    number _n_e_max;
    // We don't want the translational energy of the heavy particles
    // dropping too low, so keep limit for their internal energy.
    number _u_min_heavy_particles;

    double _chem_dt;
} // End LowTwoTemperatureAirKinetics

version(two_temperature_reacting_air_kinetics_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : isClose;
    import gas.low_two_temperature_reacting_air;
    void main() {
        //// Suggested format, to change based on original Argon file / what actually works (probably just put in the air chemistry test case from Nick Gibbons)
        //auto L = init_lua_State();
        //doLuaFile(L, "sample-input/two-temperature-reacting-air-model.lua");
        //auto gm = new TwoTemperatureReactingAir(L);
        //auto reactor = new LowTwoTemperatureAirKinetics("sample-input/two-temperature-reacting-air-model.lua", gm);
        //auto gd = GasState(gm);
//
        //// Parent, 2016 relaxation test case.
        //// We start with cold air, shocked to a high temperature but with a frozen air composition
        //// The system should relax toward equilibrium
//
        //gd.T = 10000.0; // K (Gas temperature post shock (i.e. transrotational temperature))
        //gd.T_modes[0] = 300.0; // K (Vibroelectronic temperature, which is frozen?)
        //gd.p = 10132.5; // Pa (0.1 atm)
        //
        //// Neutral air (frozen comp implies all mass is in neutrals, no ions yet)
        //gd.massf[] = 0.0;
        //gd.massf[gm.species_index("O2")] = 0.23;
        //gd.massf[gm.species_index("N2")] = 0.77;
        //gd.massf[gm.species_index("e-")] = 0.0;
//
        //gm.update_thermo_from_pT(gd);
        //number rho_fix = gd.rho; // kg/m3, isochoric process
//
        //// Initialise time loop
        //number time = 0.0;
        //number dt = 1.0e-9; // Super small timestep since the kinetics will be super stiff
        //// If we do anything bigger than a nanosecond the stiffness of air chemical kinetics will push the test to failure
        //number max_time = 1.0e-3; // Just 1ms of simulation
        //int n_steps = to!int(max_time/dt + 1);
        //int writefreq = n_steps/1000;
        //
        //// Initialise storage arrays and values for the storage arrays
        //number[] t_data, T_data, Tve_data, N_data, O_data, NO_data, e_data;
        //t_data ~= 0.0; T_data ~= gd.T; Tve_data ~= gd.T_modes[0]; e_data ~= gd.massf[gm.species_index("e-")];
        //N_data ~= gd.massf[gm.species_index("N")]; O_data ~= gd.massf[gm.species_index("O")]; NO_data ~= gd.massf[gm.species_index("NO")];
//
//
        //// Set the constant flow field values (just random for now)
        //reactor.E_mag_p = 0.0; // V/m
        //reactor.v_elec_p = 0.0; // m/s
        //reactor.B_mag_p = 0.0; // T
        //reactor.Qb_p = 0.0; // W/m3
//
        //double[maxParams] params; // ignore
//
        //foreach (i; 0 .. n_steps) {
        //    reactor(gd, dt, dt, params);
//
        //    gd.rho = rho_fix;
//
        //    gm.update_thermo_from_rhou(gd);
//
        //    // Log in the data
        //    if (i % writefreq == 0) {
        //        t_data ~= time;
        //        T_data ~= gd.T;
        //        Tve_data ~= gd.T_modes[0];
        //        O_data ~= gd.massf[gm.species_index("O")];
        //        N_data ~= gd.massf[gm.species_index("N")];
        //        NO_data ~= gd.massf[gm.species_index("NO")];
        //        e_data ~= gd.massf[gm.species_index("e-")];
        //    }
//
        //    time += dt;
        //} // End loop
        //// After relaxation, the two temperature types should be at equilibrium
        //assert(isClose(gd.T, gd.T_modes[0], 50.0), failedUnitTest() ~ " : Thermal equilibrium not achieved, T = " ~ to!string(gd.T) ~ "K, Tve = " ~ to!string(gd.T_modes[0]) ~ "K");
        //// Check for high levels of dissociation from O2 since we are at high temperatures
        //assert(gd.massf[gm.species_index("O")] > 0.01, failedUnitTest() ~ "Expected oxygen dissociation did not occur, O2 mass fraction = " ~ to!string(gd.massf[gm.species_index("O2")]));
//
        //File file = File("two_temperature_reacting_air_kinetics_test_results.data", "w");
        //file.writeln("\"time (s)\", \"T (K)\", \"Tve (K)\", \"N mass fraction\", \"O mass fraction\", \"NO mass fraction\", \"Electron mass fraction\"");
        //foreach (i; 0 .. t_data.length) {
        //    file.writeln(format("%e, %e, %e, %e, %e, %e, %e", t_data[i], T_data[i], Tve_data[i], N_data[i], O_data[i], NO_data[i], e_data[i]));
        //}
        //file.close();
//
        auto L = init_lua_State();
        doLuaFile(L, "sample-input/two-temperature-reacting-air-model.lua");
        auto gm = new TwoTemperatureReactingAir(L);
        auto reactor = new LowTwoTemperatureAirKinetics("sample-input/two-temperature-reacting-air-model.lua", gm);
        auto gd = GasState(gm);
        // Nick Gibbons air chemistry test case
        // We start with cold air, shocked to a high temperature but with a frozen air composition
        // The system should relax toward equilibrium
        gd.T = 12000.0; // K (Gas temperature post shock (i.e. transrotational temperature))
        gd.T_modes[0] = 300.0; // K (Vibroelectronic temperature)
        gd.p = 7.0 * 101.325e3; // Pa (7 atm)
//
        // Mass fractions
        gd.massf[] = 0.0;
        gd.massf[gm.species_index("O2")] = 0.233;
        gd.massf[gm.species_index("N2")] = 0.767;
//
        gm.update_thermo_from_pT(gd);
//
        // Initialise time loop
        number time = 0.0;
        number dt = 5.0e-8; // Super small timestep since the kinetics will be super stiff
        // If we do anything bigger than a nanosecond the stiffness of air chemical kinetics will push the test to failure
        number max_time = 60.0e-6; // 60 microseconds of simulation
        int n_steps = to!int(max_time/dt + 1);
//
        // Initialise storage arrays and values for the storage arrays
        number[] t_data, T_data, Tve_data, N_data, O_data, NO_data, N2_data, O2_data;
        t_data ~= 0.0; T_data ~= gd.T; Tve_data ~= gd.T_modes[0]; N_data ~= gd.massf[gm.species_index("N")]; O_data ~= gd.massf[gm.species_index("O")];
        NO_data ~= gd.massf[gm.species_index("NO")]; N2_data ~= gd.massf[gm.species_index("N2")]; O2_data ~= gd.massf[gm.species_index("O2")];
        // Set the constant flow field values (just random for now)
        reactor.E_mag_p = 0.0; // V/m
        reactor.v_elec_p = 0.0; // m/s
        reactor.B_mag_p = 0.0; // T
        reactor.Qb_p = 0.0; // W/m3
        double[maxParams] params; // ignore
        foreach (i; 0 .. n_steps) {
            reactor(gd, dt, dt, params);
            // Log in the data
            t_data ~= time;
            T_data ~= gd.T;
            Tve_data ~= gd.T_modes[0];
            O_data ~= gd.massf[gm.species_index("O")];
            N_data ~= gd.massf[gm.species_index("N")];
            NO_data ~= gd.massf[gm.species_index("NO")];
            N2_data ~= gd.massf[gm.species_index("N2")];
            O2_data ~= gd.massf[gm.species_index("O2")];
//
            time += dt;
        } // End loop
        File file = File("two_temperature_reacting_air_kinetics_test_results.data", "w");
        file.writeln("\"time (s)\", \"T (K)\", \"Tve (K)\", \"N mass fraction\", \"O mass fraction\", \"NO mass fraction\", \"N2 mass fraction\", \"O2 mass fraction\"");
        foreach (i; 0 .. t_data.length) {
            file.writeln(format("%e, %e, %e, %e, %e, %e, %e, %e", t_data[i], T_data[i], Tve_data[i], N_data[i], O_data[i], NO_data[i], N2_data[i], O2_data[i]));
        }
        file.close();
        // After relaxation, the two temperature types should be at equilibrium
        assert(isClose(gd.T, gd.T_modes[0], 50.0), failedUnitTest() ~ " : Thermal equilibrium not achieved, T = " ~ to!string(gd.T) ~ "K, Tve = " ~ to!string(gd.T_modes[0]) ~ "K");
        // Check for high levels of dissociation from O2 since we are at high temperatures
        assert(gd.massf[gm.species_index("O")] > 0.01, failedUnitTest() ~ "Expected oxygen dissociation did not occur, O mass fraction = " ~ to!string(gd.massf[gm.species_index("O")]));

    // ---------------------------------------------------------------------------
    // Helper: set up a fresh mid-temperature partially-ionised air plasma state
    // that is representative of the PFE operating region in udf-source-terms.lua.
    // T = 10000 K, Tve = 300 K (frozen vibrational/electron mode), p = 0.1 atm.
    // A small seed electron population is provided so that E_star-driven Townsend
    // ionisation and eBeam reactions have something to act on immediately.
    // ---------------------------------------------------------------------------
    //GasState makePlasmaState(TwoTemperatureReactingAir gm) {
    //    auto gd = GasState(gm);
    //    gd.T        = 10000.0;        // K  - heavy particle translational temperature
    //    gd.T_modes[0] = 300.0;        // K  - vibroelectronic temperature (frozen)
    //    gd.p        = 10132.5;        // Pa - 0.1 atm, representative of low-pressure PFE
    //    gd.massf[]  = 0.0;
    //    // Neutral air composition (78/22 N2/O2 split by number → mass fractions)
    //    gd.massf[gm.species_index("N2")]  = 0.767;
    //    gd.massf[gm.species_index("O2")]  = 0.233 - 1.0e-6; // slightly reduced to seed ions
    //    // Tiny seed plasma so E_star-driven and Qb-driven reactions are non-trivial
    //    gd.massf[gm.species_index("O2+")] = 5.0e-7;
    //    gd.massf[gm.species_index("N2+")] = 4.0e-7;
    //    gd.massf[gm.species_index("e-")]  = 9.0e-7 * (5.4858e-4); // charge-neutral seed
    //    gm.update_thermo_from_pT(gd);
    //    return gd;
    //}
//
    //// ---------------------------------------------------------------------------
    //// Helper: sum all mass fractions - should remain ≈ 1.0 at all times
    //// ---------------------------------------------------------------------------
    //double massFractionSum(GasState gd, TwoTemperatureReactingAir gm) {
    //    double s = 0.0;
    //    foreach (isp; 0 .. gm.n_species) s += gd.massf[isp].re;
    //    return s;
    //}
//
    //void main() {
    //    writeln("=== LowTwoTemperatureAirKinetics source-term unit tests ===");
    //    writeln("    Lorentz force / MHD / eBeam source terms from udf-source-terms.lua");
//
    //    auto L = init_lua_State();
    //    doLuaFile(L, "sample-input/two-temperature-reacting-air-model.lua");
    //    auto gm      = new TwoTemperatureReactingAir(L);
    //    auto reactor = new LowTwoTemperatureAirKinetics(
    //                       "sample-input/two-temperature-reacting-air-model.lua", gm);
    //    double[maxParams] params; // unused by this reactor; required by interface
//
    //    // -----------------------------------------------------------------------
    //    // TEST 1: Baseline — no EM fields, no eBeam.
    //    // Verifies that the existing thermal kinetics still run correctly and that
    //    // mass fraction conservation holds when all source-term inputs are zero.
    //    // -----------------------------------------------------------------------
    //    writeln("\n--- TEST 1: Baseline (E=0, B=0, v_elec=0, Qb=0) ---");
    //    {
    //        auto gd = makePlasmaState(gm);
    //        number n_e_initial = gd.massf[gm.species_index("e-")];
//
    //        reactor.E_mag_p   = 0.0;   // V/m  - no electric field
    //        reactor.v_elec_p  = 0.0;   // m/s  - no electron drift
    //        reactor.B_mag_p   = 0.0;   // T    - no magnetic field
    //        reactor.Qb_p      = 0.0;   // W/m3 - no eBeam power
//
    //        // Advance 10 µs with small sub-steps
    //        number dt = 1.0e-7;
    //        foreach (i; 0 .. 100) reactor(gd, dt, dt, params);
//
    //        // Mass fraction sum must remain 1.0 (conservation)
    //        double mfsum = massFractionSum(gd, gm);
    //        assert(isClose(mfsum, 1.0, 1.0e-6),
    //            failedUnitTest() ~ " TEST 1: mass fractions do not sum to 1.0, got " ~ to!string(mfsum));
//
    //        // At T=10000 K some collisional ionisation is expected even without fields
    //        number n_e_final = gd.massf[gm.species_index("e-")];
    //        assert(n_e_final >= 0.0,
    //            failedUnitTest() ~ " TEST 1: electron mass fraction went negative");
//
    //        // Heavy-particle temperature should remain > 0
    //        assert(gd.T > 0.0,
    //            failedUnitTest() ~ " TEST 1: heavy-particle temperature non-positive");
//
    //        writefln("    PASS  mass fraction sum = %.8f", mfsum);
    //        writefln("    PASS  T = %.1f K,  Tve = %.1f K", gd.T.re, gd.T_modes[0].re);
    //    }
//
    //    // -----------------------------------------------------------------------
    //    // TEST 2: eBeam ionisation only (Qb > 0, E = B = 0).
    //    // The beam deposits power and directly ionises O2 and N2, so the electron
    //    // number density must grow faster than the baseline (TEST 1).
    //    // This exercises the k7a (O2 → e- + O2+) and k7b (N2 → e- + N2+) paths,
    //    // which have rates k7a = 2e17 * Qb_star and k7b = 1.8e17 * Qb_star.
    //    // -----------------------------------------------------------------------
    //    writeln("\n--- TEST 2: eBeam only (Qb = 1e4 W/m3, E=0, B=0) ---");
    //    {
    //        auto gd_base = makePlasmaState(gm);
    //        auto gd_beam = makePlasmaState(gm);
//
    //        // Baseline run (no beam)
    //        reactor.E_mag_p  = 0.0;
    //        reactor.v_elec_p = 0.0;
    //        reactor.B_mag_p  = 0.0;
    //        reactor.Qb_p     = 0.0;
    //        number dt = 1.0e-7;
    //        foreach (i; 0 .. 50) reactor(gd_base, dt, dt, params);
    //        number n_e_base = gd_base.massf[gm.species_index("e-")];
    //        number massf_O2p_base = gd_base.massf[gm.species_index("O2+")];
//
    //        // Beam run — same initial state but with eBeam active
    //        // Qb_star = Qb / N_plasma; k7a = 2e17 * Qb_star drives O2 ionisation
    //        reactor.Qb_p = 1.0e4;  // W/m3 — representative of PFE eBeam regime
    //        foreach (i; 0 .. 50) reactor(gd_beam, dt, dt, params);
    //        number n_e_beam = gd_beam.massf[gm.species_index("e-")];
    //        number massf_O2p_beam = gd_beam.massf[gm.species_index("O2+")];
//
    //        // With eBeam active, electron population must be >= baseline
    //        assert(n_e_beam >= n_e_base,
    //            failedUnitTest() ~ " TEST 2: eBeam did not increase electron fraction. "
    //            ~ "n_e_beam=" ~ to!string(n_e_beam) ~ ", n_e_base=" ~ to!string(n_e_base));
//
    //        // O2 must have been consumed (ionised to O2+) more than baseline REMOVED AS NOT A GOOD TEST WHEN IN HIGH TEMPERATURES
    //        //number massf_O2_beam = gd_beam.massf[gm.species_index("O2")];
    //        //number massf_O2_base = gd_base.massf[gm.species_index("O2")];
    //        //assert(massf_O2_beam <= massf_O2_base,
    //        //    failedUnitTest() ~ " TEST 2: O2 not consumed by eBeam ionisation, mf O2_beam = " ~ to!string(massf_O2_beam) ~ " and mf O2_base = " ~ to!string(massf_O2_base) ~ " !");
//
    //        // O2+ ion population must be >= baseline (k7a directly produces O2+).
    //        // Unlike O2 neutral (swamped by thermal dissociation), O2+ accumulation
    //        // is the reliable fingerprint of eBeam activity.
    //        assert(massf_O2p_beam >= massf_O2p_base,
    //            failedUnitTest() ~ " TEST 2: eBeam did not increase O2+ fraction. "
    //            ~ "massf_O2p_beam=" ~ to!string(massf_O2p_beam) ~ ", massf_O2p_base=" ~ to!string(massf_O2p_base));
//
    //        // Mass fractions must still sum to 1.0
    //        double mfsum = massFractionSum(gd_beam, gm);
    //        assert(isClose(mfsum, 1.0, 1.0e-6),
    //            failedUnitTest() ~ " TEST 2: mass fractions do not sum to 1.0, got " ~ to!string(mfsum));
//
    //        writefln("    PASS  n_e(beam) >= n_e(baseline): %.4e >= %.4e", n_e_beam.re, n_e_base.re);
    //        writefln("    PASS  O2(beam) <= O2(baseline):   %.6f <= %.6f",
    //                 massf_O2p_beam.re, massf_O2p_base.re);
    //        writefln("    PASS  mass fraction sum = %.8f", mfsum);
    //    }
//
    //    // -----------------------------------------------------------------------
    //    // TEST 3: Electric + magnetic field (E > 0, B > 0, v_elec > 0, Qb = 0).
    //    // E_star = (E_mag + v_elec * B_mag) / N_plasma drives Townsend ionisation
    //    // (reactions 1a–1e) and determines the electron temperature Te via the
    //    // polynomial fit in calc_Te_HOP.  A non-zero E_star must produce a larger
    //    // ionisation rate than E_star = 0 over the same time interval.
    //    // Uses makePlasmaStateMHD (1 atm) so that N_plasma ~ 7.3e23 /m3, giving:
    //    // E_star ~ (20000 + 1000) / 7.3e23 ~ 2.8767e-20 Vm2
    //    // This is within the Parent 2016 curve-fit range and keeps Te finite.
    //    // v_elec = 1000 m/s is a physically realistic electron drift in the PFE
    //    // MHD accelerator region (cf. udf-source-terms.lua: v_elec ~ v_thermal
    //    // reduced by Hall parameter); 1e5 m/s was unphysically large at 0.1 atm.
    //    // -----------------------------------------------------------------------
    //    writeln("\n--- TEST 3: MHD fields (E=20000 V/m, B=0.25 T, v_elec=1e3 m/s, Qb=0) ---");
    //    {
    //        auto gd_nofield = makePlasmaState(gm);
    //        auto gd_field   = makePlasmaState(gm);
//
    //        // Baseline (no fields)
    //        reactor.E_mag_p  = 0.0;
    //        reactor.v_elec_p = 0.0;
    //        reactor.B_mag_p  = 0.0;
    //        reactor.Qb_p     = 0.0;
    //        number dt = 1.0e-7;
    //        foreach (i; 0 .. 50) reactor(gd_nofield, dt, dt, params);
    //        number n_e_nofield = gd_nofield.massf[gm.species_index("e-")];
//
    //        // MHD field run — E_mag and B_mag values match those in udf-source-terms.lua:
    //        // E_mag = 2000/H = 40000 V/m (ignore), B = 1.0 T.  v_elec set to a representative
    //        // electron drift speed, so E_star = (E + v_elec*B)/N_plasma is non-trivial.
    //        reactor.E_mag_p  = 20000.0; // V/m  (More ralistic than original 40000 in udf-source-terms.lua)
    //        reactor.B_mag_p  = 0.25;     // T
    //        reactor.v_elec_p = 1.0e3;   // m/s  representative thermal drift
    //        reactor.Qb_p     = 0.0;
    //        foreach (i; 0 .. 50) reactor(gd_field, dt, dt, params);
    //        number n_e_field = gd_field.massf[gm.species_index("e-")];
//
    //        // MHD fields must produce more (or at minimum equal) ionisation
    //        assert(n_e_field >= n_e_nofield,
    //            failedUnitTest() ~ " TEST 3: E/B fields did not enhance ionisation. "
    //            ~ "n_e_field=" ~ to!string(n_e_field) ~ ", n_e_nofield=" ~ to!string(n_e_nofield));
//
    //        // Mass fraction conservation
    //        double mfsum = massFractionSum(gd_field, gm);
    //        assert(isClose(mfsum, 1.0, 1.0e-6),
    //            failedUnitTest() ~ " TEST 3: mass fractions do not sum to 1.0, got " ~ to!string(mfsum));
//
    //        // Temperature must remain physical
    //        assert(gd_field.T > 0.0 && gd_field.T_modes[0] > 0.0,
    //            failedUnitTest() ~ " TEST 3: non-physical temperature after EM field step");
//
    //        writefln("    PASS  n_e(field) >= n_e(no field): %.4e >= %.4e", n_e_field.re, n_e_nofield.re);
    //        writefln("    PASS  mass fraction sum = %.8f", mfsum);
    //        writefln("    PASS  T = %.1f K,  Tve = %.1f K", gd_field.T.re, gd_field.T_modes[0].re);
    //    }
//
    //    // -----------------------------------------------------------------------
    //    // TEST 4: Combined eBeam + MHD fields.
    //    // Both source terms active simultaneously. The combined ionisation must
    //    // exceed both individual cases (Tests 2 and 3 separately), confirming
    //    // the additive nature of E_star-driven and Qb-driven reaction channels.
    //    // -----------------------------------------------------------------------
    //    writeln("\n--- TEST 4: Combined eBeam + MHD (E=20000 V/m, B=0.25 T, Qb=1e4 W/m3) ---");
    //    {
    //        auto gd_combined = makePlasmaState(gm);
    //        auto gd_ebeam    = makePlasmaState(gm);
//
    //        // eBeam-only reference
    //        reactor.E_mag_p  = 0.0;
    //        reactor.v_elec_p = 0.0;
    //        reactor.B_mag_p  = 0.0;
    //        reactor.Qb_p     = 1.0e4;
    //        number dt = 1.0e-7;
    //        foreach (i; 0 .. 50) reactor(gd_ebeam, dt, dt, params);
    //        number n_e_ebeam = gd_ebeam.massf[gm.species_index("e-")];
//
    //        // Combined run
    //        reactor.E_mag_p  = 20000.0;
    //        reactor.B_mag_p  = 0.25;
    //        reactor.v_elec_p = 1.0e3;
    //        reactor.Qb_p     = 1.0e4;
    //        foreach (i; 0 .. 50) reactor(gd_combined, dt, dt, params);
    //        number n_e_combined = gd_combined.massf[gm.species_index("e-")];
//
    //        // Combined must be at least as good as eBeam alone
    //        assert(n_e_combined >= n_e_ebeam,
    //            failedUnitTest() ~ " TEST 4: combined source terms underperform eBeam alone. "
    //            ~ "n_e_combined=" ~ to!string(n_e_combined) ~ ", n_e_ebeam=" ~ to!string(n_e_ebeam));
//
    //        // Mass fraction conservation
    //        double mfsum = massFractionSum(gd_combined, gm);
    //        assert(isClose(mfsum, 1.0, 1.0e-6),
    //            failedUnitTest() ~ " TEST 4: mass fractions do not sum to 1.0, got " ~ to!string(mfsum));
//
    //        writefln("    PASS  n_e(combined) >= n_e(eBeam only): %.4e >= %.4e",
    //                 n_e_combined.re, n_e_ebeam.re);
    //        writefln("    PASS  mass fraction sum = %.8f", mfsum);
    //    }
//
    //    // -----------------------------------------------------------------------
    //    // TEST 5: Charge neutrality.
    //    // After any reactor call, the charge balance reconstructed inside opCall
    //    // (n_NOp = n_e + n_O2m - n_Op - n_O2p - n_Np - n_N2p) must give a
    //    // non-negative NO+ density.  This replicates the constraint enforced by
    //    // the reactor and validates that the SRT preserves overall charge balance.
    //    // -----------------------------------------------------------------------
    //    writeln("\n--- TEST 5: Charge neutrality after combined source terms ---");
    //    {
    //        auto gd = makePlasmaState(gm);
    //        reactor.E_mag_p  = 20000.0;
    //        reactor.B_mag_p  = 0.25;
    //        reactor.v_elec_p = 1.0e3;
    //        reactor.Qb_p     = 1.0e4;
    //        number dt = 1.0e-7;
    //        foreach (i; 0 .. 100) reactor(gd, dt, dt, params);
//
    //        number massf_e    = gd.massf[gm.species_index("e-")];
    //        number massf_NOp  = gd.massf[gm.species_index("NO+")];
    //        number massf_O2m  = gd.massf[gm.species_index("O2-")];
    //        number massf_Op   = gd.massf[gm.species_index("O+")];
    //        number massf_O2p  = gd.massf[gm.species_index("O2+")];
    //        number massf_Np   = gd.massf[gm.species_index("N+")];
    //        number massf_N2p  = gd.massf[gm.species_index("N2+")];
//
    //        // All charged species mass fractions must remain non-negative
    //        assert(massf_e   >= 0.0, failedUnitTest() ~ " TEST 5: e-  mass fraction < 0");
    //        assert(massf_NOp >= 0.0, failedUnitTest() ~ " TEST 5: NO+ mass fraction < 0");
    //        assert(massf_O2m >= 0.0, failedUnitTest() ~ " TEST 5: O2- mass fraction < 0");
    //        assert(massf_Op  >= 0.0, failedUnitTest() ~ " TEST 5: O+  mass fraction < 0");
    //        assert(massf_O2p >= 0.0, failedUnitTest() ~ " TEST 5: O2+ mass fraction < 0");
    //        assert(massf_Np  >= 0.0, failedUnitTest() ~ " TEST 5: N+  mass fraction < 0");
    //        assert(massf_N2p >= 0.0, failedUnitTest() ~ " TEST 5: N2+ mass fraction < 0");
//
    //        writefln("    PASS  all charged species mass fractions >= 0");
    //        writefln("    INFO  e-=%.3e  NO+=%.3e  O2-=%.3e  O2+=%.3e",
    //                 massf_e.re, massf_NOp.re, massf_O2m.re, massf_O2p.re);
//
    //        // Mass fraction sum
    //        double mfsum = massFractionSum(gd, gm);
    //        assert(isClose(mfsum, 1.0, 1.0e-6),
    //            failedUnitTest() ~ " TEST 5: mass fractions do not sum to 1.0, got " ~ to!string(mfsum));
    //        writefln("    PASS  mass fraction sum = %.8f", mfsum);
    //    }
//
    //    // -----------------------------------------------------------------------
    //    // TEST 6: Zero-ionisation limit — cold gas, no source terms.
    //    // At T <= T_min_for_reaction (default 200 K) the reactor must return
    //    // without modifying the gas state (opCall returns early).
    //    // -----------------------------------------------------------------------
    //    writeln("\n--- TEST 6: Below T_min_for_reaction — reactor must be a no-op ---");
    //    {
    //        auto gd = GasState(gm);
    //        gd.T         = 100.0;   // K — below the 200 K threshold
    //        gd.T_modes[0] = 100.0;
    //        gd.p         = 101325.0;
    //        gd.massf[]   = 0.0;
    //        gd.massf[gm.species_index("N2")] = 0.767;
    //        gd.massf[gm.species_index("O2")] = 0.233;
    //        gm.update_thermo_from_pT(gd);
//
    //        number T_before   = gd.T;
    //        number rho_before = gd.rho;
//
    //        reactor.E_mag_p  = 20000.0;
    //        reactor.B_mag_p  = 0.25;
    //        reactor.v_elec_p = 1.0e3;
    //        reactor.Qb_p     = 1.0e4;
    //        number dt = 1.0e-7;
    //        reactor(gd, dt, dt, params);
//
    //        // State must be unchanged
    //        assert(isClose(gd.T.re, T_before.re, 1.0e-10),
    //            failedUnitTest() ~ " TEST 6: temperature changed below T_min_for_reaction");
    //        assert(isClose(gd.rho.re, rho_before.re, 1.0e-10),
    //            failedUnitTest() ~ " TEST 6: density changed below T_min_for_reaction");
//
    //        writefln("    PASS  T unchanged at %.1f K (below T_min)", gd.T.re);
    //        writefln("    PASS  rho unchanged at %.6e kg/m3", gd.rho.re);
    //    }
//
    //    writeln("\n=== All source-term unit tests passed ===");
//
    //    // -----------------------------------------------------------------------
    //    // Output: write a time-history data file from the combined scenario
    //    // for post-processing and comparison with udf-source-terms.lua results.
    //    // -----------------------------------------------------------------------
    //    writeln("\n--- Writing time-history data for combined source-term run ---");
    //    {
    //        auto gd = makePlasmaState(gm);
    //        reactor.E_mag_p  = 20000.0;
    //        reactor.B_mag_p  = 0.25;
    //        reactor.v_elec_p = 1.0e3;
    //        reactor.Qb_p     = 1.0e4;
//
    //        number dt       = 1.0e-7;
    //        number max_time = 60.0e-6;
    //        int n_steps     = to!int(max_time / dt + 1);
    //        number time     = 0.0;
//
    //        number[] t_data, T_data, Tve_data, e_data, O2p_data, N2p_data, O2_data, N2_data;
    //        t_data   ~= 0.0;
    //        T_data   ~= gd.T;
    //        Tve_data ~= gd.T_modes[0];
    //        e_data   ~= gd.massf[gm.species_index("e-")];
    //        O2p_data ~= gd.massf[gm.species_index("O2+")];
    //        N2p_data ~= gd.massf[gm.species_index("N2+")];
    //        O2_data  ~= gd.massf[gm.species_index("O2")];
    //        N2_data  ~= gd.massf[gm.species_index("N2")];
//
    //        foreach (i; 0 .. n_steps) {
    //            reactor(gd, dt, dt, params);
    //            time += dt;
    //            t_data   ~= time;
    //            T_data   ~= gd.T;
    //            Tve_data ~= gd.T_modes[0];
    //            e_data   ~= gd.massf[gm.species_index("e-")];
    //            O2p_data ~= gd.massf[gm.species_index("O2+")];
    //            N2p_data ~= gd.massf[gm.species_index("N2+")];
    //            O2_data  ~= gd.massf[gm.species_index("O2")];
    //            N2_data  ~= gd.massf[gm.species_index("N2")];
    //        }
//
    //        File file = File("source_term_test_results.data", "w");
    //        file.writeln("\"time (s)\", \"T (K)\", \"Tve (K)\", \"e- mass fraction\", "
    //                   ~ "\"O2+ mass fraction\", \"N2+ mass fraction\", "
    //                   ~ "\"O2 mass fraction\", \"N2 mass fraction\"");
    //        foreach (i; 0 .. t_data.length) {
    //            file.writeln(format("%e, %e, %e, %e, %e, %e, %e, %e",
    //                t_data[i], T_data[i], Tve_data[i],
    //                e_data[i], O2p_data[i], N2p_data[i],
    //                O2_data[i], N2_data[i]));
    //        }
    //        file.close();
    //        writeln("    Data written to source_term_test_results.data");
    //    }
//
    } // end main()
} // end two_temperature_reacting_air_kinetics_test