/**
 * low_two_temperature_reacting_air.d
 *
 * Low two-temperature reacting air based on the model of:
 * Bernard Parent and Mikhail N. Shneider and Sergey O. Macheret
 * "Detailed Modeling of Plasmas for Computational Aerodynamics"
 * AIAA Journal, Vol. 54, No. 3, 898 (2016); doi 10.2514/1.J054624
 *
 * Authors: Thomas Condon
 * Inspired By: Daniel Smith, Rory Kelly and Peter J.
 * Version: 03-February-2026: initial development
 */

module gas.low_two_temperature_reacting_air;

import std.math;
import std.stdio;
import std.string;
import std.file;
import std.json;
import std.conv;
import util.lua;
import util.lua_service;
import ntypes.complex;
import nm.number;
import std.algorithm : canFind;


import gas.gas_model;
import gas.gas_state;
import gas.physical_constants;
import gas.diffusion.viscosity;
import gas.diffusion.therm_cond;
import gas.thermo.perf_gas_mix_eos; // Added based on two_temperature_air
import gas.thermo.cea_thermo_curves; // Added based on two_temperature_air - not needed

immutable double T_REF = 298.15; // K, reference temperature for formation enthalpies
immutable double[12] formation_enthalpies = [
    0.0,            // e-
    249175.003,     // O
    1568787.228,    // O+
    0.0,            // O2
    1171828.436,    // O2+
    0.0,            // O2-, 0.0 for now, will try to find a reliable source for this at some point. But it is shallowly researched ?
    472680.0,       // N
    1882127.624,    // N+
    0.0,            // N2
    1509508.424,    // N2+
    91271.31,       // NO
    990809.704,     // NO+ 
]; // Taken from the Eilmer species database, all in J/mol, then converted to J/kg in _del_hf

// Now, our specific gas model.
enum Species {e_minus=0, O, O_plus, O2, O2_plus, O2_minus, N, N_plus, N2, N2_plus, NO, NO_plus}
static string[] molecularSpeciesNames = ["O", "O+", "O2", "O2+", "O2-", "N", "N+", "N2", "N2+", "NO", "NO+"];

@nogc Species getSpeciesId(string name)
{
    switch (name) {
    case "N": return Species.N;
    case "O": return Species.O;
    case "N2": return Species.N2;
    case "O2": return Species.O2;
    case "NO": return Species.NO;
    case "N+": return Species.N_plus;
    case "O+": return Species.O_plus;
    case "N2+": return Species.N2_plus;
    case "O2+": return Species.O2_plus;
    case "NO+": return Species.NO_plus;
    case "e-": return Species.e_minus;
    case "O2-": return Species.O2_minus;
    default: throw new Error("invalid species name.");
    } // end switch
} // end getSpeciesId

class TwoTemperatureReactingAir: GasModel {
public:
    int[] molecularSpecies;
    this(lua_State *L) {
        type_str = "TwoTemperatureReactingAir";
        // Some parameters are fixed and some come from the gas model file.
        _is_plasma = true;
        _n_species = 12;
        _n_modes = 1;
        _species_names.length = 12;
        _species_names[Species.e_minus] = "e-";
        _species_names[Species.O] = "O";
        _species_names[Species.O_plus] = "O+";
        _species_names[Species.O2] = "O2";
        _species_names[Species.O2_plus] = "O2+";
        _species_names[Species.O2_minus] = "O2-";
        _species_names[Species.N] = "N";
        _species_names[Species.N_plus] = "N+";
        _species_names[Species.N2] = "N2";
        _species_names[Species.N2_plus] = "N2+";
        _species_names[Species.NO] = "NO";
        _species_names[Species.NO_plus] = "NO+";
        
        // All molar masses taken from PubChem
        _mol_masses.length = 12;
        _mol_masses[Species.e_minus] = 5.485799e-7; // Units are kg/mol
        _mol_masses[Species.O] = 0.015999; // Units are kg/mol
        _mol_masses[Species.O_plus] = _mol_masses[Species.O] - _mol_masses[Species.e_minus]; // Units are kg/mol
        _mol_masses[Species.O2] = 0.031999; // Units are kg/mol
        _mol_masses[Species.O2_plus] = _mol_masses[Species.O2] - _mol_masses[Species.e_minus]; // Units are kg/mol
        _mol_masses[Species.O2_minus] = _mol_masses[Species.O2] + _mol_masses[Species.e_minus]; // Units are kg/mol
        _mol_masses[Species.N] = 0.014007; // Units are kg/mol
        _mol_masses[Species.N_plus] = _mol_masses[Species.N] - _mol_masses[Species.e_minus]; // Units are kg/mol
        _mol_masses[Species.N2] = 0.028014; // Units are kg/mol
        _mol_masses[Species.N2_plus] = _mol_masses[Species.N2] - _mol_masses[Species.e_minus]; // Units are kg/mol
        _mol_masses[Species.NO] = 0.03006; // Units are kg/mol
        _mol_masses[Species.NO_plus] = _mol_masses[Species.NO] - _mol_masses[Species.e_minus]; // Units are kg/mol
        _e_mass_over_ion_mass = _mol_masses[Species.e_minus] / (_mol_masses[Species.O_plus] + _mol_masses[Species.O2_plus] + _mol_masses[Species.O2_minus] + _mol_masses[Species.N_plus] + _mol_masses[Species.N2_plus] + _mol_masses[Species.NO_plus]);
        
        // Set the id of each species in a vector which we can call later
        _species_ids.length = _n_species;
        foreach (isp, sp; _species_names) { _species_ids[isp] = getSpeciesId(sp); }

        _R.length = _n_species;
        _molef.length = _n_species;
        _particleMass.length = _n_species;
        _del_hf.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            _R[isp] = R_universal / _mol_masses[isp];
            _particleMass[isp] = _mol_masses[isp]/Avogadro_number;
            _particleMass[isp] *= 1000.0; // Changed to grams
            // Formation enthalpies of each species in J/kg is the value of enthalpy at 298.15K
            _del_hf[isp] = formation_enthalpies[isp] / _mol_masses[isp];
        }        

        // Cp in translation and rotation is fixed, hence let's define it to make gas temperature easier to calculate later
        _Cp_tr_rot.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            if (_species_ids[isp] == Species.e_minus) {
                // The electron translation is governed by T_ve,
                // so it has no energy contribution in the T_tr mode.
                _Cp_tr_rot[isp] = 0.0;
                continue;
            }
            _Cp_tr_rot[isp] = (5./2.)*_R[isp];
            if (canFind(molecularSpeciesNames, _species_names[isp])) {
                _Cp_tr_rot[isp] += _R[isp];
                molecularSpecies ~= isp;
            }
        }

        // Compute the reference vibrational energies at T_REF, which ensures u_modes[0] = 0 at T_modes[0] = T_REF
        _reference_vib_energy.length = _n_species;
        foreach (isp; 0 .. _n_species) {
            if (_species_ids[isp] == Species.O2) {
                _reference_vib_energy[isp] = _R[isp] * theta_v_O2 / (exp(theta_v_O2 / T_REF) - 1.0);
            } else if (_species_ids[isp] == Species.N2) {
                _reference_vib_energy[isp] = _R[isp] * theta_v_N2 / (exp(theta_v_N2 / T_REF) - 1.0);
            } else if (_species_ids[isp] == Species.NO) {
                _reference_vib_energy[isp] = _R[isp] * theta_v_NO / (exp(theta_v_NO / T_REF) - 1.0);
            } else {
                _reference_vib_energy[isp] = 0.0;
            }
        }

        // Finish off constructor
        lua_getglobal(L, "TwoTemperatureReactingAir");
        // Now, pull out the remaining numeric value parameters.
        _ion_tol = getDoubleWithDefault(L, -1, "ion_tol", 0.0);
        _Te_default = getDoubleWithDefault(L, -1, "Te_default", 10000.0);
        lua_pop(L, 1); // dispose of the table
        create_species_reverse_lookup();
    } // end constructor

    override string toString() const
    {
        char[] repr;
        repr ~= "TwoTemperatureReactingAir =(";
        repr ~= "species=[\"e-\", \"O\", \"O+\", \"O2\", \"O2+\", \"O2-\", \"N\", \"N+\", \"N2\", \"N2+\", \"NO\", \"NO+\"], ";
        repr ~= ", Mmass=[" ~ to!string(_mol_masses[Species.e_minus]);
        repr ~= "," ~ to!string(_mol_masses[Species.O]);
        repr ~= "," ~ to!string(_mol_masses[Species.O_plus]);
        repr ~= "," ~ to!string(_mol_masses[Species.O2]);
        repr ~= "," ~ to!string(_mol_masses[Species.O2_plus]);
        repr ~= "," ~ to!string(_mol_masses[Species.O2_minus]);
        repr ~= "," ~ to!string(_mol_masses[Species.N]);
        repr ~= "," ~ to!string(_mol_masses[Species.N_plus]);
        repr ~= "," ~ to!string(_mol_masses[Species.N2]);
        repr ~= "," ~ to!string(_mol_masses[Species.N2_plus]);
        repr ~= "," ~ to!string(_mol_masses[Species.NO]);
        repr ~= "," ~ to!string(_mol_masses[Species.NO_plus]) ~ "]";
        repr ~= ")";
        return to!string(repr);
    }

    override void update_thermo_from_pT(ref GasState Q)
    {
        number alpha = ionisation_fraction_from_mass_fractions(Q);
        if (Q.T <= 0.0 || Q.p <= 0.0 || Q.T_modes[0] <= 0.0) {
            string msg = "Temperature and/or pressure was negative for update_thermo_from_pT.";
            debug { msg ~= format("\nQ=%s\n", Q); }
            throw new GasModelException(msg);
        }
        Q.rho = update_density(Q, Q.T_modes[0]);
        Q.u = transRotEnergy(Q);
        Q.u_modes[0] = vibElecEnergy(Q, Q.T_modes[0]);
    }
    override void update_thermo_from_rhou(ref GasState Q)
    {
        number alpha = ionisation_fraction_from_mass_fractions(Q);
        if (Q.rho <= 0.0) {
            string msg = "Density was negative for update_thermo_from_rhou.";
            debug { msg ~= format("\nQ=%s\n", Q); }
            throw new GasModelException(msg);
        }
        Q.T = update_gas_temperature(Q);
        number Te = _Te_default; // in case alpha == 0.0
        if (alpha > _ion_tol) {
            Te = Q.T_modes[0];
        }
        if (Te > 500.0e3) { Te = 500.0e3; }
        //if (Te < 200.0) { Te = 200.0; }
        Q.T_modes[0] = vibElecTemperature(Q);
        Q.p = update_gas_pressure(Q, Te);
    }
    override void update_thermo_from_rhoT(ref GasState Q)
    {
        number alpha = ionisation_fraction_from_mass_fractions(Q);
        if (Q.T <= 0.0 || Q.rho <= 0.0 || Q.T_modes[0] <= 0.0) {
            string msg = "Temperature and/or density was negative for update_thermo_from_rhoT.";
            debug { msg ~= format("\nQ=%s\n", Q); }
            throw new GasModelException(msg);
        }
        number Te = Q.T_modes[0];
        Q.p = update_gas_pressure(Q, Te);

        // Calculate energy
        Q.u = transRotEnergy(Q);
        Q.u_modes[0] = vibElecEnergy(Q, Q.T_modes[0]);
    }
    override void update_thermo_from_rhop(ref GasState Q)
    {
        // Assume Q.T_modes[0] remains fixed. <- this was here in the argon file. Let's not do this? I am unsure why they did this
        number alpha = ionisation_fraction_from_mass_fractions(Q);
        if (Q.p <= 0.0 || Q.rho <= 0.0) {
            string msg = "Pressure and/or density was negative for update_thermo_from_rhop.";
            debug { msg ~= format("\nQ=%s\n", Q); }
            throw new GasModelException(msg);
        }
        Q.T = update_gas_temperature(Q);
        number Te = _Te_default; // in case alpha == 0.0
        if (alpha > _ion_tol) {
            Te = Q.T_modes[0];
        }
        if (Te > 500.0e3) { Te = 500.0e3; }
        //if (Te < 200.0) { Te = 200.0; }
        Q.T_modes[0] = vibElecTemperature(Q);
        
        // Calculate energy
        Q.u = transRotEnergy(Q);
        Q.u_modes[0] = vibElecEnergy(Q, Q.T_modes[0]);
    }
    override void update_thermo_from_ps(ref GasState Q, number s)
    {
        throw new GasModelException("update_thermo_from_ps not implemented in TwoTemperatureReactingAir.");
    }
    override void update_thermo_from_hs(ref GasState Q, number h, number s)
    {
        throw new GasModelException("update_thermo_from_hs not implemented in TwoTemperatureReactingAir.");
    }
    override void update_sound_speed(ref GasState Q)
    {
        if (Q.T <= 0.0 || Q.T_modes[0] < 0.0) {
            string msg = "Temperature was negative for update_sound_speed.";
            debug { msg ~= format("\nQ=%s\n", Q); }
            throw new GasModelException(msg);
        }
        // Let's assume a frozen speed of sound, as it could be expected that the sound wave passes fast enough
        // that the chemistry and vibrational energy does not have time to change (i.e. frozen)
        // We will calculate it assuming a boltzmann distribution and hence use the equipartition theorem to calculate gamma
        // Although in hindsight it was suggested I could just use 1.4 for gamma, as stated in the air species-database in Eilmer
        number Cp_mix = 0.0;
        number Cv_mix = 0.0;

        foreach (i; 0 .. _n_species) {
            // Same logic as my original hard coded transrotational energy mon/di atomic calculations before switching to
            // the logic in the two_temperature_air.d model
            number Cv_i = 0.0;
            if (_species_ids[i] == Species.O2 || _species_ids[i] == Species.N2 || _species_ids[i] == Species.NO || _species_ids[i] == Species.O2_plus ||
            _species_ids[i] == Species.O2_minus || _species_ids[i] == Species.N2_plus || _species_ids[i] == Species.NO_plus) {
                Cv_i = 5.0/2.0 * _R[i]; // Diatomic molecules
            } else {
                Cv_i = 3.0/2.0 * _R[i]; // Monatomic molecules including electrons
            }
            number Cp_i = Cv_i + _R[i];
            Cv_mix += Q.massf[i] * Cv_i;
            Cp_mix += Q.massf[i] * Cp_i;
        }
        number gamma_frozen = Cp_mix / Cv_mix;
        Q.a = sqrt(gamma_frozen * Q.p / Q.rho);
    }
    override void update_trans_coeffs(ref GasState Q)
    {
        // We will assume the properties of air molecules, technically ignoring the ionisation. Taken from the Eilmer
        // species-database, will we use the Sutherland law from "Table 1-2, White (2006)" and "Table 1-3, White (2006)"
        double mu_ref = 1.716e-5;
        double T_ref = 273.0;
        double S_mu = 111.0;
        double k_ref = 0.0241;
        double S_k = 194.0;
        Q.mu = (mu_ref * pow(Q.T/T_ref, 3.0/2.0) * (T_ref + S_mu)) / (Q.T + S_mu);
        Q.k = (k_ref * pow(Q.T/T_ref, 3.0/2.0) * (T_ref + S_k)) / (Q.T + S_k);
        Q.k_modes[0] = 0.0;
    }
    override number dudT_const_v(in GasState Q)
    {
        number Cv_mix = 0.0;

        foreach (i; 0 .. _n_species) {
            // Same logic as my original hard coded transrotational energy mon/di atomic calculations before switching to
            // the logic in the two_temperature_air.d model
            number Cv_i = 0.0;
            if (_species_ids[i] == Species.O2 || _species_ids[i] == Species.N2 || _species_ids[i] == Species.NO || _species_ids[i] == Species.O2_plus ||
            _species_ids[i] == Species.O2_minus || _species_ids[i] == Species.N2_plus || _species_ids[i] == Species.NO_plus) {
                Cv_i = 5.0/2.0 * _R[i]; // Diatomic molecules
            } else {
                Cv_i = 3.0/2.0 * _R[i]; // Monatomic molecules including electrons
            }
            Cv_mix += Q.massf[i] * Cv_i;
        }
        return Cv_mix;
    }
    override number dhdT_const_p(in GasState Q)
    {
        number Cp_mix = 0.0;

        foreach (i; 0 .. _n_species) {
            // Same logic as my original hard coded transrotational energy mon/di atomic calculations before switching to
            // the logic in the two_temperature_air.d model
            number Cv_i = 0.0;
            if (_species_ids[i] == Species.O2 || _species_ids[i] == Species.N2 || _species_ids[i] == Species.NO || _species_ids[i] == Species.O2_plus ||
            _species_ids[i] == Species.O2_minus || _species_ids[i] == Species.N2_plus || _species_ids[i] == Species.NO_plus) {
                Cv_i = 5.0/2.0 * _R[i]; // Diatomic molecules
            } else {
                Cv_i = 3.0/2.0 * _R[i]; // Monatomic molecules including electrons
            }
            number Cp_i = Cv_i + _R[i];
            Cp_mix += Q.massf[i] * Cp_i;
        }
        return Cp_mix;
    }
    override number dpdrho_const_T(in GasState Q)
    {
        return gas_constant(Q)*Q.T;
    }
    override number gas_constant(in GasState Q)
    {
        return mass_average(Q, _R);
    }
    override number internal_energy(in GasState Q)
    {
        return Q.u + Q.u_modes[0];
    }
    override number enthalpy(in GasState Q)
    {
        return Q.u + Q.u_modes[0] + Q.p/Q.rho;
    }
    override number entropy(in GasState Q)
    {
        throw new GasModelException("entropy not implemented in TwoTemperatureReactingAir.");
    }
    override void balance_charge(ref GasState Q)
    {
        // Let's sum all of the ion mass fractions
        number ion_massf = Q.massf[Species.O_plus] + Q.massf[Species.O2_plus] - Q.massf[Species.O2_minus] +
                           Q.massf[Species.N_plus] + Q.massf[Species.N2_plus] + Q.massf[Species.NO_plus];
        
        // Now ensure that this ion mass fraction multiplied by electron mass over the ion mass gives us the electron mass fraction
        Q.massf[Species.e_minus] = ion_massf * _e_mass_over_ion_mass;
    }
    @nogc number ionisation_fraction_from_mass_fractions(ref const(GasState) Q)
    {
        number ions = (Q.massf[Species.O_plus] / _mol_masses[Species.O_plus]) +
                      (Q.massf[Species.O2_plus] / _mol_masses[Species.O2_plus]) -
                      (Q.massf[Species.O2_minus] / _mol_masses[Species.O2_minus]) +
                      (Q.massf[Species.N_plus] / _mol_masses[Species.N_plus]) +
                      (Q.massf[Species.N2_plus] / _mol_masses[Species.N2_plus]) +
                      (Q.massf[Species.NO_plus] / _mol_masses[Species.NO_plus]);
        number atoms = (Q.massf[Species.O] / _mol_masses[Species.O]) +
                       (Q.massf[Species.O2] / _mol_masses[Species.O2]) +
                       (Q.massf[Species.N] / _mol_masses[Species.N]) +
                       (Q.massf[Species.N2] / _mol_masses[Species.N2]) +
                       (Q.massf[Species.NO] / _mol_masses[Species.NO]);
        return ions/(ions+atoms);
    }
    @nogc number vibElecEnergy(in GasState Q, number Te)
    {
        // Calculate electron/vibrational energy
        number e_ve = 0.0;
        number theta_v = 0.0;
        foreach(i; 0 .. _n_species) {
            // Electron transitional energy
            if (_species_ids[i] == Species.e_minus) {
                // Referenced to T_REF
                e_ve += Q.massf[i] * ((3.0/2.0) * _R[i] * (Te - T_REF) + _del_hf[i] - (_R[i] * T_REF));
                continue;
            }
            // Vibrational energy (Parent, 2025) -> e_v = (R_i*theta_i)/(exp(theta_i/T_v) - 1) (which is just the classic definition)
            // Then this needs to be in reference to T_REF vibrational energies
            if (_species_ids[i] == Species.O2) theta_v = theta_v_O2;
            else if (_species_ids[i] == Species.N2) theta_v = theta_v_N2;
            else if (_species_ids[i] == Species.NO) theta_v = theta_v_NO;

            if (theta_v > 0.0) {
                // Referenced to T_REF vibrational energy
                number e_v = _R[i] * theta_v / (exp(theta_v / Te) - 1.0) - _reference_vib_energy[i];
                e_ve += Q.massf[i] * e_v;
            }            
        }
        return e_ve;
    }
    @nogc number vibElecSpecHeatConstV(in GasState Q, number Te)
    {
        // Calculate Cv for vibroelectronic modes
        // Cv = dE/dT
        number Cv_ve = 0.0;
        number theta_v = 0.0;
        foreach(i; 0 .. _n_species) {
            // Electron specific heat (constant)
            if (_species_ids[i] == Species.e_minus) {
                Cv_ve += Q.massf[i] * (3.0/2.0) * _R[i];
                continue;
            }
            // Vibrational specific heat (harmonic oscillator)
            // Cv_vib = R * (theta_v/T)^2 * exp(theta_v/T) / [exp(theta_v/T) - 1]^2
            if (_species_ids[i] == Species.O2) theta_v = theta_v_O2;
            else if (_species_ids[i] == Species.N2) theta_v = theta_v_N2;
            else if (_species_ids[i] == Species.NO) theta_v = theta_v_NO;

            if (theta_v > 0.0) {
                // Implementation of harmonic oscillator
                number x = theta_v / Te;
                number exp_x = exp(x);
                number Cv_v = _R[i] * x * x * exp_x / ((exp_x - 1.0) * (exp_x - 1.0));
                Cv_ve += Q.massf[i] * Cv_v;
            }
        }
        return Cv_ve;
    }

    @nogc number vibElecTemperature(in GasState Q)
    {
        int MAX_ITERATIONS = 20;
        double TOL = 1.0e-6;
        double E_TOL = 0.01; // J

        // Use current T_modes[0] as initial guess
        number T_guess = Q.T_modes[0];
        number f_guess = vibElecEnergy(Q, T_guess) - Q.u_modes[0];
        
        // Check if initial guess is good enough
        if (fabs(f_guess) < E_TOL) {
            return Q.T_modes[0];
        }

        // Newton-Raphson iteration
        int count = 0;
        number Cv, dT;
        foreach (iter; 0 .. MAX_ITERATIONS) {
            Cv = vibElecSpecHeatConstV(Q, T_guess);
            dT = -f_guess / Cv;
            T_guess += dT;
            
            if (fabs(dT) < TOL) {
                break;
            }
            
            f_guess = vibElecEnergy(Q, T_guess) - Q.u_modes[0];
            count++;
        }

        if (count == MAX_ITERATIONS) {
            string msg = "The 'vibElecTemperature' function failed to converge.\n";
            debug {
                msg ~= format("The final value for T_modes[0] was: %12.6f\n", T_guess);
                msg ~= "The supplied GasState was:\n";
                msg ~= Q.toString() ~ "\n";
            }
            throw new GasModelException(msg);
        }

        return T_guess;
    }


private:
    // Thermodynamic constants
    // Added new for reacting air
    double R_uni = 8.31451; // J/mol/K
    double _ion_tol = 0.0; // may be over-ridden from value in Lua file
    double _Te_default = 10000.0; // may be over-ridden for a specific situation
    double _e_mass_over_ion_mass; // set once mol masses are set in constructor
    double theta_v_N2 = 3353; // K
    double theta_v_NO = 2719; // K
    double theta_v_O2 = 2240; // K
    number[] _del_hf;
    number[] _R;
    number[] _Cp_tr_rot;
    Species[] _species_ids;
    number[] _reference_vib_energy; // Vibrational energy at T_REF, for consistent referencing

    // From air gas model
    PerfectGasMixEOS _pgMixEOS;
    double _R_U_cal = 1.987; // cal/(mole.K)
    number[] _molef; // will be getting mole-fractions from outside, so may be complex
    number[] _s;
    double[] _particleMass;
    //double[] _R;
    CEAThermoCurve[] _curves;
    //number[] _del_hf;
    //double[] _Cp_tr_rot;
    number[][] _A_11, _B_11, _C_11, _D_11, _Delta_11, _alpha;
    number[][] _A_22, _B_22, _C_22, _D_22, _Delta_22, _mu;
    //Species[] _species_ids;

    @nogc
    number update_density(ref GasState Q, number Te)
    {
        // Unfortunately, as defined in Parent, 2025. We cannot just simply use the universal gas constant
        number R_heavy = 0.0; number R_elec = 0.0;
        // Iterate through the species to sum up their partial pressures and find their gas constants
        foreach (i; 0 .. _n_species) {
            if (_species_ids[i] == Species.e_minus) {
                R_elec += Q.massf[i] * _R[i];
            } else {
                R_heavy += Q.massf[i] * _R[i];
            }
        }
        return Q.p / (R_heavy*Q.T + R_elec*Te);
    }

    @nogc
    number update_gas_temperature(ref GasState Q)
    {
        // We can compute T by direct inversion since the Cp in
        // translation and rotation are fully excited,
        // and, as such, constant.
        number sumA = 0.0;
        number sumB = 0.0;
        number updated_temp = 0.0;
        foreach (isp; 0 .. _n_species) {
            if (_species_ids[isp] == Species.e_minus) continue;
            sumA += Q.massf[isp]*(_Cp_tr_rot[isp]*T_REF - _del_hf[isp]);
            sumB += Q.massf[isp]*(_Cp_tr_rot[isp] - _R[isp]);
        }
        updated_temp = (Q.u + sumA)/sumB;
        return updated_temp;
    }

    @nogc
    number heavyParticleGasConstant(ref GasState Q)
    {
        number Rmix = 0.0;
        foreach (isp; 0 .. _n_species) {
            if (_species_ids[isp] == Species.e_minus) continue;
            Rmix += Q.massf[isp]*_R[isp];
        }
        return Rmix;
    }

    @nogc
    number update_gas_pressure(ref GasState Q, number Te)
    {
        // Now we calculate pressure, definitions converted to a hard coded format from two_temperature_air
        number Rmix = heavyParticleGasConstant(Q);
        number R_elec = 0.0;
        number updated_pressure = 0.0;
        foreach(i; 0 .. _n_species) {
            if (_species_ids[i] == Species.e_minus) {
                R_elec = _R[i];
            }
        }
        updated_pressure = Q.rho*((Rmix*Q.T) + (Q.massf[Species.e_minus]*R_elec*Te));
        return updated_pressure;
    }

    @nogc number transRotEnergy(in GasState Q)
    {
        // Calculate thermal energy of the heavy particles
        number e_tr_rot = 0.0;
        foreach (isp; 0 .. _n_species) {
            if (_species_ids[isp] == Species.e_minus) continue;
            number h_tr_rot = _Cp_tr_rot[isp]*(Q.T - T_REF) + _del_hf[isp];
            e_tr_rot += Q.massf[isp]*(h_tr_rot - _R[isp]*Q.T);
        }
        return e_tr_rot;
    }
} // end class

// Unit test of the basic gas model

version(two_temperature_reacting_air_test) {
    import std.stdio;
    import util.msg_service;
    import std.math : isClose;
    int main() {
        lua_State* L = init_lua_State();
        doLuaFile(L, "sample-data/two-temperature-reacting-air-model.lua");
        auto gm = new TwoTemperatureReactingAir(L);
        lua_close(L);
        auto gd = GasState(gm);

        // Testing standard atmosphere on new gas model (300K, 1atm, 77% N2, 23% O2)
        gd.p = 101325.0;
        gd.T = 300.0;
        gd.T_modes[0] = 300.0;
        gd.massf[] = 0.0;
        gd.massf[Species.N2] = 0.77;
        gd.massf[Species.O2] = 0.23;

        assert(gm.n_modes == 1, failedUnitTest() ~ " : Incorrect number of modes");
        assert(gm.n_species == 12, failedUnitTest() ~ " : Incorrect number of species");
        assert(isClose(gd.p, 101325.0, 1.0e-6), failedUnitTest() ~ " : Pressure initialisation failed");
        assert(isClose(gd.T, 300.0, 1.0e-6), failedUnitTest() ~ " : Temperature initialisation failed");
        assert(isClose(gd.massf[Species.N2], 0.77, 1.0e-6), failedUnitTest() ~ " : N2 mass fraction initialisation failed");
        assert(isClose(gd.massf[Species.O2], 0.23, 1.0e-6), failedUnitTest() ~ " : O2 mass fraction initialisation failed");

        gm.update_thermo_from_pT(gd);
        gm.update_sound_speed(gd);

        assert(isClose(gm.R(gd), 287.0, 1.0), failedUnitTest() ~ " : Gas constant mismatch " ~ to!string(gm.R(gd)) ~ " J/kg/K");
        number my_rho = 101325.0 / (287.0 * 300.0); // Ideal gas law for air, where R =~ 287 J/kg/K
        assert(isClose(gd.rho, my_rho, 0.01), failedUnitTest() ~ " : Air density mismatch " ~ to!string(gd.rho) ~ " kg/m^3");
        assert(isClose(gd.a, 347.0, 5.0), failedUnitTest() ~ "Sound speed mismatch" ~ to!string(gd.a) ~ " m/s");

        // Testing high temperature non-equilibrium case on new gas model (T = 10000K, Tve = 4000K, 10kPa, 50% O, 50% N)
        gd.p = 10000.0;
        gd.T = 10000.0;
        gd.T_modes[0] = 4000.0;

        gd.massf[] = 0.0;
        gd.massf[Species.N] = 0.5;
        gd.massf[Species.O] = 0.5;

        gm.update_thermo_from_pT(gd);
        number h_test = gm.enthalpy(gd);

        assert(gd.u > 0.0, failedUnitTest() ~ " : Internal energy" ~ to!string(gd.u) ~ " J should be positive");
        assert(h_test > gd.u, failedUnitTest() ~ " : Enthalpy " ~ to!string(h_test) ~ " J/kg should be greater than internal energy");
        
        return 0; // The Argon gas file manually returned 0 because it was failing I assume?
    }
}