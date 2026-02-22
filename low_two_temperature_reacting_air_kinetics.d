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
    public double bulk_v_elec_p = 0.0; // Electron Velocity (m/s)
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
            _integration_method = getStringWithDefault(L, -1, "integration_method", "RK4");
            _T_min_for_reaction = getDoubleWithDefault(L, -1, "T_min_for_reaction", 200.0); // set to 200.0K from 3000.0K
            _Te_default = getDoubleWithDefault(L, -1, "Te_default", 10000.0); // User can set to 200.0K for low temperature applications
            _n_step_suggest = getIntWithDefault(L, -1, "n_step_suggest", 10);
            _max_iter_newton = getIntWithDefault(L, -1, "max_iter_newton", 30);
        }
        lua_close(L);
    }

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
            number bulk_v_elec = to!number(this.bulk_v_elec_p); // Electron Velocity (m/s)
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
            _u_min_heavy_particles = 3.0 / 2.0 * R * _T_min_for_reaction; // Just taken from Argon file but barely used/relevant
            number u_available = u_total - _u_min_heavy_particles; // Ditto

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

            // Calculate necessary characteristics of the plasma for the rate coefficients
            // E* = |E + Ve x B|/N
            number E_star = (E_mag + (bulk_v_elec * B_mag)) / N_plasma; // Reduced effective electric field in electron reference frame (Vm^2)
            // Qb* = Qb/N
            number Qb_star = Qb / N_plasma; // Ratio between electron beam power per unit volume and total number density of plasma (W)

            // Start with the suggested time step size
            int NumberSteps = cast(int) fmax(floor(tInterval/dtSuggest), 1.0);
            number[12] S;
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

                    switch (_integration_method) {
                    case "Forward_Euler":
                        foreach (n; 0 .. NumberSteps) {
                            // F calls update_reaction_rate which will set recent_eer (elastic energy rate)
                            number[12] myF = F(S, Q, E_star, Qb_star, N_plasma);
                            foreach (i; 0 .. 12) {
                                S[i] += _chem_dt * myF[i];
                            }

                            // Since the test cases are 0D, and the power deposition of the beam is constantly applied to this box,
                            // the beam energy should probably be influencing the total energy of the actual gas in this reactor.
                            // In one of the old pfe job scripts and source term file, it is said that there is energy addition because there is electron additions
                            // NOTE, REMOVE THE QB TERM WHEN DOING ACTUAL SIMULATIONS, THE SOURCE TERM FILE SHOULD ACCOUNT FOR THIS ITSELF
                            u_total += (Qb / Q.rho + recent_eer / Q.rho) * _chem_dt;

                            // Update Q to reflect the new energy transrotational energy
                            update_Q_from_state_vector(S, Q, E_star, N_plasma);

                            if (!state_vector_is_within_limits(S)) {
                                string msg = "State vector not within limits";
                                debug { msg ~= format("\n   y=[%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g]", S[0], S[1], S[2], S[3], S[4], S[5], S[6], S[7], S[8], S[9], S[10], S[11]); }
                                throw new ThermochemicalReactorUpdateException(msg);
                            }
                        } // end foreach n
                        break;
                    case "RK4": // Same as argon RK4 but with more care around the energy exchange
                        number [12] k1, k2_in, k2, k3_in, k3, k4_in, k4;
                        foreach (n; 0 .. NumberSteps) {
                            // Store the initial energy and density
                            number u_total_initial = u_total;
                            number rho_initial = Q.rho;

                            // Start RK4
                            number[12] myF = F(S, Q, E_star, Qb_star, N_plasma);

                            // Store the initial eer after updating the state vector with F()
                            number eer_initial = recent_eer;

                            foreach (i; 0 .. 12) { k1[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 11) { k2_in[i] = fmax(S[i] + k1[i]/2.0, 0.0); }
                            k2_in[11] = S[11] + k1[11]/2.0;
                            //foreach (i; 0 .. 13) { k2_in[i] = S[i] + k1[i]/2.0; }
                            myF = F(k2_in, Q, E_star, Qb_star, N_plasma);
                            foreach (i; 0 .. 12) { k2[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 11) { k3_in[i] = fmax(S[i] + k2[i]/2.0, 0.0); }
                            k3_in[11] = S[11] + k2[11]/2.0;
                            //foreach (i; 0 .. 13) { k3_in[i] = S[i] + k2[i]/2.0; }
                            myF = F(k3_in, Q, E_star, Qb_star, N_plasma);
                            foreach (i; 0 .. 12) { k3[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 11) { k4_in[i] = fmax(S[i] + k3[i], 0.0); }
                            k4_in[11] = S[11] + k3[11];
                            //foreach (i; 0 .. 13) { k4_in[i] = S[i] + k3[i]; }
                            myF = F(k4_in, Q, E_star, Qb_star, N_plasma);
                            foreach (i; 0 .. 12) { k4[i] = _chem_dt * myF[i]; }
                            foreach (i; 0 .. 11) { S[i] = fmax(S[i] + 1.0/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]), 0.0); }
                            S[11] += 1.0/6.0*(k1[11]+2.0*k2[11]+2.0*k3[11]+k4[11]);
                            //foreach (i; 0 .. 13) { S[i] += 1.0/6.0*(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i]); }

                            // Follow the same energy increase as in the forward euler case
                            // NOTE, REMOVE THE QB TERM WHEN DOING ACTUAL SIMULATIONS, THE SOURCE TERM FILE SHOULD ACCOUNT FOR THIS ITSELF
                            u_total = u_total_initial + (Qb / rho_initial + eer_initial / rho_initial) * _chem_dt;

                            // Update Q to reflect the new energy transrotational energy
                            update_Q_from_state_vector(S, Q, E_star, N_plasma);

                            if (!state_vector_is_within_limits(S)) {
                                string msg = "State vector not within limits";
                                debug { msg ~= format("\n    y=[%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g,%g]", S[0], S[1], S[2], S[3], S[4], S[5], S[6], S[7], S[8], S[9], S[10], S[11]); }
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
                    Q.u = S[11];
                    Q.u_modes[0] = fmax(u_total - Q.u, 0.0);

                    // Reconstruct the other parts of the flow state
                    // Utilise charge neutrality
                    n_NOp = n_e + n_O2m - n_Op - n_O2p - n_Np - n_N2p;
                    if (n_NOp < 0.0) { n_NOp = to!number(0.0); }
                
                    // Now to update the mass fractions, the electron temperature, thermodynamic behaviour of the gas, and the sound speed !
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
                    _gmodel.update_thermo_from_rhou(Q);
                    _gmodel.update_sound_speed(Q);

                    finished_integration = true;
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
        number bulk_v_elec = to!number(this.bulk_v_elec_p); // Electron Velocity (m/s)
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
        number E_star = (E_mag + (bulk_v_elec * B_mag)) / N_plasma;
        number Qb_star = Qb / N_plasma;

        // Determine the necessary rates to update the source terms
        number[12] myF = update_reaction_rates(Q, E_star, Qb_star);
        number energy_leaving_trans = myF[11];

        // Map the local enum indicies to Eilmer's dynamic indices and convert dn/dt
        // to the species production rates (kg/(s*m^3)) if the temperature is hot enough
        if (Q.T > _T_min_for_reaction) {
            // Species production rate source terms
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
            // Energy source terms
            source[12] = energy_leaving_trans * Q.rho; // Convert to W/m^3
            source[13] = (Qb + recent_eer) - source[12]; // Energy conservation for the vibroelectronic energy, hopefully interfaces okay with source term file?
            // Fairly certain the source term file will account for this for the kinetics file, so simply get rid of the Qb term in that line
        } else {
            foreach (i; 0 .. 14) { source[i] = 0.0; }
        }
    }

    private:
    @nogc
    number[12] F(ref const(number[12]) S, ref GasState Q, number E_star, number Qb_star, number N_plasma) {
        update_Q_from_state_vector(S, Q, E_star, N_plasma); // Update the gas state
        return update_reaction_rates(Q, E_star, Qb_star); // Update the rates based on this gas state and feed it back
    }

    @nogc
    void update_Q_from_state_vector(ref const(number[12]) S, ref GasState Q, number E_star, number N_plasma) {
        // Update the GasState from the state vector, that way we always work the rate calculations
        // from a physically realisable state
        number n_e = S[0]; number n_O = S[1]; number n_Op = S[2]; number n_O2 = S[3];
        number n_O2p = S[4]; number n_O2m = S[5]; number n_N = S[6]; number n_Np = S[7];
        number n_N2 = S[8]; number n_N2p = S[9]; number n_NO = S[10];

        // Utilise charge neutrality
        number n_NOp = n_e + n_O2m - n_Op - n_O2p - n_Np - n_N2p;
        if (n_NOp < 0.0) { n_NOp = to!number(0.0); }

        if (n_e < 0.0) {
            // Do not let the number of electrons go negative
            string msg = "Electron number density tried to go negative.";
            throw new Exception(msg);
        }
        Q.u = S[11];
        Q.u_modes[0] = fmax(u_total - Q.u, 0.0);

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
    number calc_ve(number Te) {
        // Average thermal velocity based on Maxwell-Boltzmann distirbution statistical mechanics
        return sqrt((8.0*k_b*Te)/(_m_e*to!double(PI)));
    }

    @nogc
    number[12] update_reaction_rates(in GasState Q, number E_star, number Qb_star) {
        // Compute the rate of change of the state vector for the reactions
        // Definition of this new rate of change state vector is:
        // [0] -> [10] dn/dt of e-, O, O+, ... , NO (excluding NO+ since it is a dependent variable)
        // [11] energy leaving translational mode
        number[12] S_dash; S_dash[0] = 0.0; S_dash[1] = 0.0; S_dash[2] = 0.0;
        S_dash[3] = 0.0; S_dash[4] = 0.0; S_dash[5] = 0.0; S_dash[6] = 0.0;
        S_dash[7] = 0.0; S_dash[8] = 0.0; S_dash[9] = 0.0; S_dash[10] = 0.0;
        S_dash[11] = 0.0;

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
        number n_ions = n_Op + n_O2p - n_O2m + n_Np + n_N2p + n_NOp; // Flag this for potential source of error,
        // Need to check definition of the n_ions parameter used in the energy exchange calculations later

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

        reaction_rates[0] = k1a * n_e * n_N2; // 1a: e- + N2 -> N2+ + e- + e-

        reaction_rates[1] = k1b * n_e * n_O2; // 1b: e- + O2 -> O2+ + e- + e-

        reaction_rates[2] = k1c * n_e * n_NO; // 1c: e- + NO -> NO+ + e- + e-

        reaction_rates[3] = k1d * n_e * n_O; // 1d: e- + O -> O+ + e- + e-

        reaction_rates[4] = k1e * n_e * n_N; //1e: e- + N -> N+ + e- + e-

        reaction_rates[5] = k2a * n_e * n_O2p; // 2a: e- + O2+ -> O + O

        reaction_rates[6] = k2b * n_e * n_N2p; // 2b: e- + N2+ -> N + N

        reaction_rates[7] = k2c * n_e * n_NOp; // 2c: e- + NO+ -> N + O

        reaction_rates[8] = k2d * n_e * n_e * n_Np; // 2d: e- + e- + N+ -> e- + N (3-Body)

        reaction_rates[9] = k2e * n_e * n_e * n_Op; // 2e: e- + e- + O+ -> e- + O (3-Body)

        reaction_rates[10] = k3a * n_O2m * n_N2p; // 3a: O2- + N2+ -> O2 + N2

        reaction_rates[11] = k3b * n_O2m * n_O2p; // 3b: O2- + O2+ -> O2 + O2

        reaction_rates[12] = k4a * n_O2m * n_N2p * n_N2; // 4a: O2- + N2+ + N2 -> O2 + N2 + N2 (3-Body)
        
        reaction_rates[13] = k4b * n_O2m * n_O2p * n_N2; // 4b: O2- + O2+ + N2 -> O2 + O2 + N2 (3-Body)

        reaction_rates[14] = k4c * n_O2m * n_N2p * n_O2; // 4c: O2- + N2+ + O2 -> O2 + N2 + O2 (3-Body)

        reaction_rates[15] = k4d * n_O2m * n_O2p * n_O2; // 4d: O2- + O2+ + O2 -> O2 + O2 + O2 (3-Body)

        reaction_rates[16] = k5a * n_e * n_O2 * n_O2; // 5a: e- + O2 + O2 -> O2- + O2

        reaction_rates[17] = k5b * n_e * n_O2 * n_N2; // 5b: e- + O2 + N2 -> O2- + N2

        reaction_rates[18] = k6 * n_O2m * n_O2; // 6: O2- + O2 -> e- + O2 + O2

        reaction_rates[19] = k7a * n_O2; // 7a: O2 -> e- + O2+

        reaction_rates[20] = k7b * n_N2; // 7b: N2 -> e- + N2+

        reaction_rates[21] = k8a * n_O2 * n_O2; // 8a: O2 + O2 -> 2O + O2

        reaction_rates[22] = k8b * n_O2 * n_N2; // 8b: O2 + N2 -> 2O + N2

        reaction_rates[23] = k8c * n_O2 * n_O; // 8c: O2 + O -> 3O

        reaction_rates[24] = k8d * n_O2 * n_N; // 8d: O2 + N → O + O + N

        reaction_rates[25] = k8e * n_O2 * n_NO; // 8e: O2 + NO → O + O + NO

        reaction_rates[26] = k8f * n_O2 * n_NOp; // 8f: O2 + NO+ → O + O + NO+

        reaction_rates[27] = k8g * n_N2 * n_O2; // 8g: N2 + O2 -> 2N + O2

        reaction_rates[28] = k8h * n_N2 * n_N2; // 8h: N2 + N2 -> 2N + N2

        reaction_rates[29] = k8i * n_N2 * n_O; // 8i: N2 + O -> 2N + O

        reaction_rates[30] = k8j * n_N2 * n_N; // 8j: N2 + N -> 3N

        reaction_rates[31] = k8k * n_N2 * n_NO; // 8k: N2 + NO → N + N + NO

        reaction_rates[32] = k8l * n_N2 * n_NOp; // 8l: N2 + NO+ → N + N + NO+

        reaction_rates[33] = k8m * n_NO * n_O; // 8m: NO + O → N + O + O

        reaction_rates[34] = k8n * n_NO * n_O2; // 8n: NO + O2 → N + O + O2

        reaction_rates[35] = k8o * n_NO * n_N; // 8o: NO + N → N + O + N

        reaction_rates[36] = k8p * n_NO * n_N2; // 8p: NO + N2 → N + O + N2

        reaction_rates[37] = k8q * n_NO * n_NOp; // 8q: NO + NO+ → N + O + NO+

        reaction_rates[38] = k8r * n_O * n_N2; // 8r: O + N2 -> N + NO

        reaction_rates[39] = k8s * n_O * n_NO; // 8s: O + NO -> N + O2

        reaction_rates[40] = k8t * n_N * n_NO; // 8t: N + NO -> O + N2

        reaction_rates[41] = k8u * n_N * n_O2; // 8u: N + O2 -> O + NO

        reaction_rates[42] = k9a * n_O * n_O * n_O2; // 9a: O + O + O2 -> 2O2

        reaction_rates[43] = k9b * n_O * n_O * n_N2; // 9b: O + O + N2 -> O2 + N2

        reaction_rates[44] = k9c * n_O * n_O * n_O; // 9c: O + O + O -> O2 + O

        reaction_rates[45] = k9d * n_N * n_N * n_O2; // 9d: N + N + O2 -> N2 + O2

        reaction_rates[46] = k9e * n_N * n_N * n_N2; // 9e: N + N + N2 -> 2N2

        reaction_rates[47] = k9f * n_N * n_N * n_O; // 9f: N + N + O -> N2 + O

        reaction_rates[48] = k9g * n_N * n_N * n_N; // 9g: N + N + N -> N2 + N

        reaction_rates[49] = k10 * n_N * n_O; // 10: N + O -> NO+ + e-

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
        // Energy moving from '-1' to the vibroelectronic mode due to reactions and collisions
        // '-1' Translational/Rotational
        // '0' Vibrational (electronic)
        // '1' Free Electron (electronic) -- This is a 2T model so we don't deal with this, I think it was recommended not to since I only had 6 weeks

        // Now we calculate the elastic electron-heavy particle energy exhange based on the Appleton-Bray expression, using Gnoffo 1989 and Imamura 2018 parameters.
        // This is mainlny based on the implementation from Nick Gibbons in two_temperature_air_kinetics.d
        number Q_e_O, Q_e_N, Q_e_O2, Q_e_N2, Q_e_NO;

        // Calculate collision cross sectional areas for each of the neutrals
        // Based on Gnoffo(1989) curve fits of electron-neutral energy exchange
        Q_e_O = 1.2e-20 + 1.7e-24*Te - 2.0e-28*Te^^2; // Higher order constant defined as -2.0e-29 in Gnoffo, but Roshan changed to -2.0e-28 ?
        Q_e_O2 = 2.0e-20 + 6.0e-24*Te;
        Q_e_N = 5.0e-20;
        Q_e_N2 = 7.5e-20 + 5.5e-24*Te - 1.0e-28*Te^^2;
        Q_e_NO = 1.0e-19;

        // Calculate the thermal electron velocity
        number v_elec = calc_ve(Te);

        // Calculate the collisions of each neutral with electrons, then sum them up for total electron-neutral collision, based on Imamura (2018)
        number nu_e_O = 4.0/3.0 * Q_e_O * n_O * v_elec;
        number nu_e_O2 = 4.0/3.0 * Q_e_O2 * n_O2 * v_elec;
        number nu_e_N = 4.0/3.0 * Q_e_N * n_N * v_elec;
        number nu_e_N2 = 4.0/3.0 * Q_e_N2 * n_N2 * v_elec;
        number nu_e_NO = 4.0/3.0 * Q_e_NO * n_NO * v_elec;
        number nu_e_n = nu_e_O + nu_e_O2 + nu_e_N + nu_e_N2 + nu_e_NO;

        // Calculate the collision of each ion with electrons, based on Imamura (2018) and two_temperature_air_kinetics.d
        number Le = 4.0 * to!number(PI) * eps_zero * k_b * Te / q_e / q_e;
        number Lambda = log(Le*Le*Le/to!number(PI) / fmax(n_e, 1.0));
        number nu_e_i = (6.0 * to!number(PI)) * pow(q_e*q_e / (12.0 * to!number(PI) * eps_zero * k_b * Te), 2) * Lambda * n_ions * v_elec;

        // Accumulate the elastic energy exchange contributions from neutrals and ions
        number elastic_energy_constants = 3.0 * n_e * _m_e * k_b * (Te - T); // Terms not related to the species, but it is clear that when Te > T, the energy will flow into transrotational and vice versa (good)

        // The neutral species have corresponding collisions, so we must manually tally up the elastic energy for the neutrals
        number een_neutrals = (nu_e_O / (_gmodel.mol_masses[i_O]/AvN))
                            + (nu_e_O2 / (_gmodel.mol_masses[i_O2]/AvN))
                            + (nu_e_N / (_gmodel.mol_masses[i_N]/AvN))
                            + (nu_e_N2 / (_gmodel.mol_masses[i_N2]/AvN))
                            + (nu_e_NO / (_gmodel.mol_masses[i_NO]/AvN));

        // The ionic species have one overarching collision, so we can just divide by the average mass of the ions instead of tallying
        number ions_avg_mass = ((_gmodel.mol_masses[i_Op] + _gmodel.mol_masses[i_O2p]
                            + _gmodel.mol_masses[i_O2m] + _gmodel.mol_masses[i_Np]
                            + _gmodel.mol_masses[i_N2p] + _gmodel.mol_masses[i_NOp]) / AvN) / 6.0;
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

        // Store this elastic energy for use in the ODE loops
        recent_eer = elastic_energy; // W/m3

        // Now to determine V-T relaxation
        number VT_relax_rate = 0.0;
        int[3] vib_species = [i_O2, i_N2, i_NO]; 
        number[3] theta_v_array = [theta_v_O2, theta_v_N2, theta_v_NO];

        foreach (idx; 0 .. 3) { // Only want O2, N2, and NO (the main heavy neutral species which have characteristic vibrational temperatures)
            int isp = vib_species[idx];
            number theta_v = theta_v_array[idx];

            // Add a guard against low vibroelectronic temp which would cause incorrect application
            if (Q.T_modes[0] < 50.0) continue;

            // Calculate relaxation time based on Millikan-White and Park high-T correction - could be changed to Parent, 2016 relaxation time
            number tau_vt = calculate_tau_vt(Q, isp);

            // 'Equilibrium' vibrational energy from transrotational temperature
            number e_v_star = (R_universal / _gmodel.mol_masses[isp]) * theta_v / (exp(theta_v / Q.T) - 1.0);

            // Current vibrational energy from vibroelectronic temperature
            number e_v = (R_universal / _gmodel.mol_masses[isp]) * theta_v / (exp(theta_v / Q.T_modes[0]) - 1.0);

            // Landau-Teller relaxation rate
            VT_relax_rate += Q.massf[isp] * (e_v_star - e_v) / tau_vt;
        }

        // Apply this relaxation to S_dash to feed the internal redistribution of energy. Qb and elastic_energy are external sources for u_total
        S_dash[11] = -VT_relax_rate; // J/kg/s
        return S_dash;
    } // end rates

    @nogc
    bool state_vector_is_within_limits(ref number[12] S)
    {
        bool result = true;
        if (S[0] < 0.0) { result = false; }
        if (S[0] > _n_e_max) { result = false; }
        //if (S[11] < _u_min_heavy_particles) { result = false; }
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
    // Set during update_reaction_rates for use in the ODE loops
    number recent_eer; // Elastic energy rate

    double _chem_dt;
} // End LowTwoTemperatureAirKinetics

version(two_temperature_reacting_air_kinetics_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : isClose;
    import gas.low_two_temperature_reacting_air;
    string unit_test_case = "sensitivity investigation";
    void main() {
        if ( unit_test_case == "relaxation") {
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

            // Mass fractions
            gd.massf[] = 0.0;
            gd.massf[gm.species_index("O2")] = 0.233;
            gd.massf[gm.species_index("N2")] = 0.767;

            gm.update_thermo_from_pT(gd);

            // Initialise time loop
            number time = 0.0;
            number dt = 5.0e-8; // mall timestep since the kinetics will be super stiff
            // If we do anything bigger than a nanosecond the stiffness of air chemical kinetics will push the test to failure
            number max_time = 60.0e-6; // 60 microseconds of simulation
            int n_steps = to!int(max_time/dt + 1);

            // Initialise storage arrays and values for the storage arrays
            number[] t_data, T_data, Tve_data, N_data, O_data, NO_data, N2_data, O2_data;
            t_data ~= 0.0; T_data ~= gd.T; Tve_data ~= gd.T_modes[0]; N_data ~= gd.massf[gm.species_index("N")]; O_data ~= gd.massf[gm.species_index("O")];
            NO_data ~= gd.massf[gm.species_index("NO")]; N2_data ~= gd.massf[gm.species_index("N2")]; O2_data ~= gd.massf[gm.species_index("O2")];
            // Set the constant flow field values (just random for now)
            reactor.E_mag_p = 0.0; // V/m
            reactor.bulk_v_elec_p = 0.0; // m/s
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
        } // End relaxation test case


        else if ( unit_test_case == "energy investigation") {
            // Special test to see how species behave in a low temperature but high energy 0D flow field chemical case
            auto L = init_lua_State();
            doLuaFile(L, "sample-input/two-temperature-reacting-air-model.lua");
            auto gm = new TwoTemperatureReactingAir(L);
            auto reactor = new LowTwoTemperatureAirKinetics("sample-input/two-temperature-reacting-air-model.lua", gm);
            auto gd = GasState(gm);

            number dt = 1.0e-9; // Super small timestep since the kinetics will be super stiff
            // If we do anything bigger than a nanosecond the stiffness of air chemical kinetics will (probably) push the test to failure?
            number max_time = 25.0e-6; // 45 microseconds of simulation
            int n_steps = to!int(max_time/dt + 1);

            // Set the constant flow field values for different scenarios to see how the plasma behaves
            number[4][4] test_cases; // Row 0-3 are the cases. Column 0 is E_mag_p (V/m), column 1 is bulk_v_elec_p (m/s), column 2 is B_mag_p (T), column 3 is Qb_p (W/m3)
            // Case 1
            test_cases[0][0] = 0.0; test_cases[0][1] = 2000.0; test_cases[0][2] = 1.0; test_cases[0][3] = 30.0e6; // Comparison case
            test_cases[1][0] = 15000.0; test_cases[1][1] = 2000.0; test_cases[1][2] = 1.0; test_cases[1][3] = 30.0e06;
            test_cases[2][0] = 17500.0; test_cases[2][1] = 2000.0; test_cases[2][2] = 1.0; test_cases[2][3] = 30.0e06;
            test_cases[3][0] = 20000.0; test_cases[3][1] = 2000.0; test_cases[3][2] = 1.0; test_cases[3][3] = 30.0e06;
            foreach (testcase; 0 .. 4) {
                writeln("Test case " ~ to!string(testcase) ~ " of 3");
                reactor.E_mag_p = test_cases[testcase][0]; // V/m
                reactor.bulk_v_elec_p = test_cases[testcase][1]; // m/s
                reactor.B_mag_p = test_cases[testcase][2]; // T
                reactor.Qb_p = test_cases[testcase][3]; // W/m3

                // Set gas state initial conditions for each test
                // We want slightly above room temp air and a massive amount of energy
                gd.T = 500.0; // K (Room temperature (i.e. transrotational temperature))
                gd.T_modes[0] = 500.0; // K (Vibroelectronic temperature)
                gd.p = 2000.0; // Pa (Not sure why, this was just in one of Roshan's test cases - I believe it is tied to a Parent, 2016 test)

                // Mass fractions
                gd.massf[] = 0.0;
                number Xi = 1.0e-4; // Ionization factor

                double W_e = 5.4858e-7;  // Electron molar mass, kg/mol
                // Compute mixture molar mass the same way pfe.lua does
                double M_mixt = 0.79*(1.0-2.0*Xi)*28.0e-3 + 0.21*(1.0-2.0*Xi)*32.0e-3
                            + 0.79*Xi*28.0e-3 + 0.21*Xi*32.0e-3 + Xi*W_e;
                // Compute mass fractions also the same way
                gd.massf[gm.species_index("N2")]  = 0.79*(1.0-2.0*Xi) * 28.0e-3 / M_mixt;
                gd.massf[gm.species_index("O2")]  = 0.21*(1.0-2.0*Xi) * 32.0e-3 / M_mixt;
                gd.massf[gm.species_index("N2+")] = 0.79*Xi * 28.0e-3 / M_mixt;
                gd.massf[gm.species_index("O2+")] = 0.21*Xi * 32.0e-3 / M_mixt;
                gd.massf[gm.species_index("e-")]  = Xi * W_e / M_mixt;

                gm.update_thermo_from_pT(gd);

                // Initialise time loop
                number time = 0.0;

                // Initialise storage arrays and values for the storage arrays
                number[] t_data, T_data, Tve_data, N_data, O_data, NO_data, N2_data, O2_data, Np_data, Op_data, O2p_data, O2m_data, N2p_data, NOp_data, e_data;
                t_data ~= 0.0; T_data ~= gd.T; Tve_data ~= gd.T_modes[0]; N_data ~= gd.massf[gm.species_index("N")]; O_data ~= gd.massf[gm.species_index("O")];
                NO_data ~= gd.massf[gm.species_index("NO")]; N2_data ~= gd.massf[gm.species_index("N2")]; O2_data ~= gd.massf[gm.species_index("O2")];
                Np_data ~= gd.massf[gm.species_index("N+")]; Op_data ~= gd.massf[gm.species_index("O+")]; O2p_data ~= gd.massf[gm.species_index("O2+")];
                O2m_data ~= gd.massf[gm.species_index("O2-")]; N2p_data ~= gd.massf[gm.species_index("N2+")]; NOp_data ~= gd.massf[gm.species_index("NO+")];
                e_data ~= gd.massf[gm.species_index("e-")];

                double[maxParams] params; // Ignore, is used in opCall but we use public variables instead
                    foreach (i; 0 .. n_steps) {
                        reactor(gd, dt, dt, params);
                        // Log the data
                        t_data ~= time;
                        T_data ~= gd.T;
                        Tve_data ~= gd.T_modes[0];
                        O_data ~= gd.massf[gm.species_index("O")];
                        N_data ~= gd.massf[gm.species_index("N")];
                        NO_data ~= gd.massf[gm.species_index("NO")];
                        N2_data ~= gd.massf[gm.species_index("N2")];
                        O2_data ~= gd.massf[gm.species_index("O2")];
                        Np_data ~= gd.massf[gm.species_index("N+")];
                        Op_data ~= gd.massf[gm.species_index("O+")];
                        O2p_data ~= gd.massf[gm.species_index("O2+")];
                        O2m_data ~= gd.massf[gm.species_index("O2-")];
                        N2p_data ~= gd.massf[gm.species_index("N2+")];
                        NOp_data ~= gd.massf[gm.species_index("NO+")];
                        e_data ~= gd.massf[gm.species_index("e-")];

                        time += dt;
                    } // End time loop for specific test case
                File file = File("two_temperature_reacting_air_kinetics_test_results_v" ~ to!string(testcase) ~".data", "w");
                file.writeln("\"time (s)\", \"T (K)\", \"Tve (K)\", \"N mass fraction\", \"O mass fraction\", \"NO mass fraction\", \"N2 mass fraction\", \"O2 mass fraction\", \"N+ mass fraction\", \"O+ mass fraction\", \"O2+ mass fraction\", \"O2- mass fraction\", \"N2+ mass fraction\", \"NO+ mass fraction\", \"e- mass fraction\"");
                foreach (i; 0 .. t_data.length) {
                    file.writeln(format("%e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e, %e", t_data[i], T_data[i], Tve_data[i], N_data[i], O_data[i], NO_data[i], N2_data[i], O2_data[i], Np_data[i], Op_data[i], O2p_data[i], O2m_data[i], N2p_data[i], NOp_data[i], e_data[i]));
                }
                file.close();
            }
        } // End energy investigation test case

        else if ( unit_test_case == "sensitivity investigation") {
            // Test to see how the electron mass fraction after 15 microseconds of 0D simulation behaves with different electric fields, initial temperatures,
            // and ratio of transrotational and vibroelectronic temperature

            auto L = init_lua_State();
            doLuaFile(L, "sample-input/two-temperature-reacting-air-model.lua");
            auto gm = new TwoTemperatureReactingAir(L);
            auto reactor = new LowTwoTemperatureAirKinetics("sample-input/two-temperature-reacting-air-model.lua", gm);
            auto gd = GasState(gm);

            // Setup parameter ranges
            number[4] T_ratios = [1.0, 10.0, 20.0, 30.0];
            
            // Initialise electric field from 0 to 25000 V/m every 1000 V/m
            number[] E_fields;
            number E_min = 0.0;
            number E_max = 25000.0;
            number E_n_steps = 250;
            for (number E = E_min; E <= E_max; E += (E_max - E_min) / E_n_steps) {
                E_fields ~= E;
            }

            // Fixed parameters
            number max_time = 25.0e-6; // Write in data at 25 microseconds
            number dt = 1.0e-9;
            int n_steps = to!int(max_time/dt + 1);

            foreach (testcase, ratio; T_ratios) {        
                writeln("Test case " ~ to!string(testcase) ~ " of 3");
                File file = File("two_temperature_reacting_air_kinetics_test_results_v" ~ to!string(testcase) ~".data", "w");
                file.writeln("\"E_field (V/m)\", \"e_mass_frac\", \"O2-_mass_fraction\"");        
                foreach (E_mag; E_fields) {
                    // Initial plasma parameters and gas state initial condition
                    reactor.bulk_v_elec_p = 2000.0; // m/s
                    reactor.B_mag_p = 1.0; // T
                    reactor.Qb_p = 30.0e6; // W/m3
                    gd.T_modes[0] = 300.0;
                    gd.T = gd.T_modes[0] * ratio; // Set Tve based on the ratio
                    gd.p = 2000.0; // Pa (Not sure why, this was just in one of Roshan's test cases - I believe it is tied to a Parent, 2016 test)
                    // Each electric field
                    reactor.E_mag_p = E_mag;

                    // Mass fractions
                    gd.massf[] = 0.0;
                    number Xi = 1.0e-4; // Ionization factor

                    double W_e = 5.4858e-7;  // Electron molar mass, kg/mol
                    // Compute mixture molar mass the same way pfe.lua does
                    double M_mixt = 0.79*(1.0-2.0*Xi)*28.0e-3 + 0.21*(1.0-2.0*Xi)*32.0e-3
                                + 0.79*Xi*28.0e-3 + 0.21*Xi*32.0e-3 + Xi*W_e;
                    // Compute mass fractions also the same way
                    gd.massf[gm.species_index("N2")]  = 0.79*(1.0-2.0*Xi) * 28.0e-3 / M_mixt;
                    gd.massf[gm.species_index("O2")]  = 0.21*(1.0-2.0*Xi) * 32.0e-3 / M_mixt;
                    gd.massf[gm.species_index("N2+")] = 0.79*Xi * 28.0e-3 / M_mixt;
                    gd.massf[gm.species_index("O2+")] = 0.21*Xi * 32.0e-3 / M_mixt;
                    gd.massf[gm.species_index("e-")]  = Xi * W_e / M_mixt;

                    gm.update_thermo_from_pT(gd);

                    double[maxParams] params; // Ignore, is used in opCall but we use public variables instead
                    // Run the time loop, we only need the last value
                    number time = 0.0;
                    foreach (i; 0 .. n_steps) {
                        reactor(gd, dt, dt, params);
                        time += dt;
                    }
                    file.writeln(format("%e, %e, %e", E_mag, gd.massf[gm.species_index("e-")], gd.massf[gm.species_index("O2-")]));
                }
                file.close();
            }
        }

        else { throw new Exception("Invalid low temperature air plasma test case"); }

    } // end main()
} // end two_temperature_reacting_air_kinetics_test