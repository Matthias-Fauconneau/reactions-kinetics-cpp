//#!/bin/sh -c c++ main.cpp -lnekRK -locca -o /tmp/nekrk && gpu-on && /tmp/nekrk; gpu-off; sensors | rg GPU
//watchexec -e cpp,okl 'sudo make -j9 -l8 install -C ~/nekRK/build && c++ main.cpp -lnekrk -locca -g -o /tmp/nekrk gpu-on && /tmp/nekrk; gpu-off; sensors | rg GPU'
#include <cstddef>
using usize = size_t ;
using f64 = double;
#include <string>
using namespace std;
#include <occa.h>
#include <nekrk.h>

auto main(int argc, const char **argv) -> int {
	using namespace nekRK;
	occa::device device((string)"{mode: 'Serial', 'kernel/compiler': 'clang', 'kernel/compiler_language': 'C++'}");
	f64 pressure = 101325.;
	f64 temperature = 1000.;
	const usize species_len = 53;
	f64 mass_fractions[species_len] = {};
	/// !!! /!\ WARNING /!\ !!! C introduces undefined behavior by allowing to use uninitialized arrays and not zero initializing by default !!!
	for(int i=0; i<species_len; i+=1) { mass_fractions[i] = 0.; }
	mass_fractions[/*O2*/3] = 2./5.;
	mass_fractions[/*CH4*/13] = 1./5.;
	mass_fractions[/*Ar*/48] = 2./5.;

	const usize len = 1;
	//const usize len = 100000;

	//auto reaction_kinetics = compile("Fuego", "grimech30", device, {}, /*Reference{*/1., 1., 1., 1., 1., 1./*}*/, len);
	init("Fuego", "grimech30", device, {},
			1., // temperature
			1., // pressure
			1., // length
			1., // velocity
			mass_fractions/*_for_reference_heat_capacity*/, // mass_fractions[species_len]
			len); // Sets global OCCA reaction kinetics settings

	f64 temperatures[len];
	for(usize i=0; i<len; i+=1) temperatures[i] = temperature;

	f64 mass_fractions_arrays[species_len][len];
	for(usize i=0; i<species_len; i++) for(usize id=0; id<len; id+=1) mass_fractions_arrays[i][id] = mass_fractions[i];

  auto device_temperature = device.malloc<f64>(len);
	device_temperature.copyFrom(temperatures);
	auto device_contiguous_mass_fractions_arrays = device.malloc<f64>(species_len*len);
	{
		f64 contiguous_mass_fractions_arrays[species_len*len];
		for(usize i=0; i<species_len; i+=1) for(int id=0; id<len; id+=1) contiguous_mass_fractions_arrays[i*len+id] = mass_fractions_arrays[i][id];
		device_contiguous_mass_fractions_arrays.copyFrom(contiguous_mass_fractions_arrays);
	}
	auto device_contiguous_specie_mass_production_rate_arrays = device.malloc<f64>(species_len*len);
	auto device_heat_release_rate = device.malloc<f64>(len);

	// Uses global OCCA reaction kinetics settings
	production_rates(pressure,
									 device_temperature, device_contiguous_mass_fractions_arrays,
									/*&mut*/ device_contiguous_specie_mass_production_rate_arrays, /*&mut*/ device_heat_release_rate);

	f64 specie_mass_production_rate[species_len][len]; // V*dtω*W
	f64 contiguous_specie_mass_production_rate_arrays[species_len*len];
	device_contiguous_specie_mass_production_rate_arrays.copyTo(contiguous_specie_mass_production_rate_arrays);
	for(usize i=0; i<species_len; i+=1) for(usize id=0; id<len; id+=1) specie_mass_production_rate[i][id] = contiguous_specie_mass_production_rate_arrays[i*len+id];

	f64 heat_release_rate[len];
	device_heat_release_rate.copyTo(heat_release_rate);

	for(usize i=0; i<species_len-1; i++) printf("%e ", specie_mass_production_rate[i][0]);
	return 0;
}
