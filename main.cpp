//#!/bin/sh -c c++ main.cpp -lnekRK -locca -o /tmp/nekrk && gpu-on && /tmp/nekrk; gpu-off; sensors | rg GPU
//watchexec -e cpp,okl 'sudo make -j9 -l8 install -C ~/nekRK/build && c++ main.cpp -lnekrk -locca -g -o /tmp/nekrk -ferror-limit=1 && gpu-on && rm ~/occa -r && /tmp/nekrk; gpu-off; sensors | rg GPU'
#include <cstddef>
using usize = size_t ;
using f64 = double;
#include <string>
using namespace std;
#include <occa.h>
#include <nekrk.h>
#include <cassert>

#define AMREX_GPU_HOST_DEVICE
#define AMREX_FORCE_INLINE
namespace amrex {
	using Real = double;
	//using Vector<typename T> = std::vector<T>;
}
void CKRHOY(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  rho);
void CKCPBS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  cpbs);
void CKMMWY(amrex::Real *  y,  amrex::Real *  wtm);
/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTY(amrex::Real *  x,  amrex::Real *  y);
void get_imw(amrex::Real *imw_new);
/* Returns the vector of strings of species names */
void CKSYMS_STR(std::vector<std::string>& kname);

auto main(int argc, const char **argv) -> int {
	std::vector<std::string> species_names; /*=*/ CKSYMS_STR(/*&mut*/ species_names);

	using namespace nekRK;
	occa::device device((string)"{mode: 'Serial', 'kernel/compiler': 'clang', 'kernel/compiler_language': 'C++'}");
	f64 pressure = 101325.;
	f64 temperature = 1000.;
	const usize species_len = 53;
	f64 mole_fractions[species_len] = {};
	/// !!! /!\ WARNING /!\ !!! C introduces undefined behavior by allowing to use uninitialized arrays and not zero initializing by default !!!
	for(int i=0; i<species_len; i+=1) { mole_fractions[i] = 0.; }
	mole_fractions[/*O2*/3] = 2./5.;
	mole_fractions[/*CH4*/13] = 1./5.;
	mole_fractions[/*Ar*/48] = 2./5.;
	f64 mass_fractions[species_len] = {};
	for(int i=0; i<species_len; i+=1) { mass_fractions[i] = 0.; }
	// mass_fractions ( mole_fractions )
	CKXTY(mole_fractions, mass_fractions);
	for(usize i=0; i<species_len-1; i++) printf("%f ", mass_fractions[i]);
	printf("\n");

	const usize states_len = 1;
	//const usize states_len = 100000;

	//auto reaction_kinetics = compile("Fuego", "grimech30", device, {}, /*Reference{*/1., 1., 1., 1., 1., 1./*}*/, states_len);
	f64 reference_temperature = temperature;
	f64 reference_pressure = pressure;
	f64 reference_length = 1.;
	f64 reference_velocity = 1.;
	init("Fuego", "grimech30", device, {},
			reference_temperature,
			reference_pressure,
			reference_length,
			reference_velocity,
			mass_fractions/*_for_reference_heat_capacity*/, // mass_fractions[species_len]
			states_len); // Sets global OCCA reaction kinetics settings

	f64 temperatures[states_len];
	for(usize i=0; i<states_len; i+=1) temperatures[i] = temperature / reference_temperature;

	f64 mass_fractions_arrays[species_len][states_len];
	for(usize i=0; i<species_len; i++) for(usize id=0; id<states_len; id+=1) mass_fractions_arrays[i][id] = mass_fractions[i];

  auto device_temperature = device.malloc<f64>(states_len);
	device_temperature.copyFrom(temperatures);
	auto device_contiguous_mass_fractions_arrays = device.malloc<f64>(species_len*states_len);
	{
		f64 contiguous_mass_fractions_arrays[species_len*states_len];
		for(usize i=0; i<species_len; i+=1) for(int id=0; id<states_len; id+=1) contiguous_mass_fractions_arrays[i*states_len+id] = mass_fractions_arrays[i][id];
		device_contiguous_mass_fractions_arrays.copyFrom(contiguous_mass_fractions_arrays);
	}
	auto device_contiguous_specie_mass_production_rate_arrays = device.malloc<f64>(species_len*states_len);
	auto device_heat_release_rate = device.malloc<f64>(states_len);

	// Uses global OCCA reaction kinetics settings
	assert(device_contiguous_specie_mass_production_rate_arrays.length() == species_len*states_len);
	production_rates(pressure / reference_pressure,
									 device_temperature, device_contiguous_mass_fractions_arrays,
									/*&mut*/ device_contiguous_specie_mass_production_rate_arrays, /*&mut*/ device_heat_release_rate);

	f64 scaled_specie_mass_production_rate[species_len][states_len]; // V*dtω*W
	f64 contiguous_specie_mass_production_rate_arrays[species_len*states_len];
	device_contiguous_specie_mass_production_rate_arrays.copyTo(contiguous_specie_mass_production_rate_arrays);
	for(usize i=0; i<species_len; i+=1) for(usize id=0; id<states_len; id+=1) scaled_specie_mass_production_rate[i][id] = contiguous_specie_mass_production_rate_arrays[i*states_len+id];

	f64 heat_release_rate[states_len];
	device_heat_release_rate.copyTo(heat_release_rate);

	double P_ref = reference_pressure;
	double T_ref = reference_temperature;
	double* Y_ref = mass_fractions;
	double L_ref = reference_length;
	double v_ref = reference_velocity;
	double t_ref = L_ref / v_ref; // time = length / velocity
	// density (pressure, temperature, mass fractions)
	double rho_ref=0.; CKRHOY(&P_ref, &T_ref, Y_ref, &rho_ref);
	// specific heat capacity at constant pressure in mass units (temperature, mass fractions)
	double cpbs_ref=0.; CKCPBS(&T_ref, Y_ref, &cpbs_ref);
	// mean molecular weight (mass fractions)
	double mmw_ref=0.; CKMMWY(Y_ref, /*wtm?:*/ &mmw_ref);
	// 1/ molecular weight
	double rcp_molar_mass[species_len]; get_imw(rcp_molar_mass);

	f64 specie_production_rate[species_len][states_len]; // V*dtω*W
	for(usize i=0; i<species_len; i+=1) for(usize id=0; id<states_len; id+=1) {
		specie_production_rate[i][id] = scaled_specie_mass_production_rate[i][id] / ( /*reference_time*/t_ref * /*reference_mean_molecular_weight*/mmw_ref / /*reference_density*/rho_ref ) * 1e8 /*/W*/ * rcp_molar_mass[i] /* mass -> quantity */;
	}

	for(usize i=0; i<species_len; i++) if (specie_production_rate[i][0] != 0.) { printf("%s %e ", species_names[i].c_str(), specie_production_rate[i][0]); }
	return 0;
}
