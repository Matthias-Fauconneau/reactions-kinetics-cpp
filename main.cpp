//#!/bin/sh -c c++ main.cpp -lnekRK -locca -o /tmp/nekrk && gpu-on && /tmp/nekrk; gpu-off; sensors | rg GPU
//watchexec -e cpp,okl 'sudo make -j9 -l8 install -C ~/nekRK/build && c++ main.cpp -lnekrk -locca -g -o /tmp/nekrk gpu-on && /tmp/nekrk; gpu-off; sensors | rg GPU'
#include <cstddef>
using usize = size_t ;
using f64 = double;
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

	auto reaction_kinetics = setup("Fuego", "grimech30", device, {}, /*Reference{*/1., 1., 1., 1., 1., 1./*}*/, len);

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
	auto device_contiguous_wdot_arrays = device.malloc<f64>(species_len*len);
	auto device_heat_release_rate = device.malloc<f64>(len);

	production_rates(reaction_kinetics, pressure,
									 device_temperature, device_contiguous_mass_fractions_arrays,
									/*&mut*/ device_contiguous_wdot_arrays, /*&mut*/ device_heat_release_rate);

	f64 /*dt_mass*/wdot[species_len][len]; // V*dtω*W: mass rate
	f64 contiguous_wdot_arrays[species_len*len];
	device_contiguous_wdot_arrays.copyTo(contiguous_wdot_arrays);
	for(usize i=0; i<species_len; i+=1) for(usize id=0; id<len; id+=1) wdot[i][id] = contiguous_wdot_arrays[i*len+id];

	f64 heat_release_rate[len];
	device_heat_release_rate.copyTo(heat_release_rate);

	for(usize i=0; i<species_len-1; i++) printf("%e ", wdot[i][0]);
	return 0;
}
