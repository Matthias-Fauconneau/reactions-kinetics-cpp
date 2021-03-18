//#!/bin/sh -c c++ main.cpp -lnekRK -locca -o /tmp/nekrk && gpu-on && /tmp/nekrk; gpu-off; sensors | rg GPU
#include <cstddef>
using usize = size_t ;
using f64 = double;
#include <nekrk.h>

auto main(int argc, const char **argv) -> int {
	using namespace nekRK;
	occa::device device((string)"{mode: 'Serial', 'kernel/compiler': 'clang', 'kernel/compiler_language': 'C++'}");
	auto reaction_kinetics = setup("Fuego", "grimech30", device, {}, /*Reference{*/1., 1., 1., 1., 1., 1./*}*/);
	f64 pressure = 101325.;
	f64 temperature = 1000.;
	const usize species_len = 53;
	f64 mass_fractions[species_len] = {};
	/// !!! /!\ WARNING /!\ !!! C introduces undefined behavior by allowing to use uninitialized arrays and not zero initializing by default !!!
	for(int i=0; i<species_len; i+=1) { mass_fractions[i] = 0.; }
	mass_fractions[/*O2*/3] = 2./5.;
	mass_fractions[/*CH4*/13] = 1./5.;
	mass_fractions[/*Ar*/48] = 2./5.;

	//const usize N = 1;
	const usize N = 100000;

	f64 pressures[N];
	for(usize i=0; i<N; i+=1) pressures[i] = pressure;

	f64 temperatures[N];
	for(usize i=0; i<N; i+=1) temperatures[i] = temperature;

	f64 mass_fractions_arrays[species_len][N];
	for(usize i=0; i<species_len; i++) for(usize id=0; id<N; id+=1) mass_fractions_arrays[i][id] = mass_fractions[i];

	const usize len = N;
  auto device_pressure = device.malloc<f64>(len);
	device_pressure.copyFrom(pressures);
	auto device_temperature = device.malloc<f64>(len);
	device_temperature.copyFrom(temperatures);
	auto device_contiguous_mass_fractions_arrays = device.malloc<f64>(species_len*len);
	{
		f64 contiguous_mass_fractions_arrays[species_len*len];
		for(usize i=0; i<species_len; i+=1) for(int id=0; id<len; id+=1) contiguous_mass_fractions_arrays[i*len+id] = mass_fractions_arrays[i][id];
		device_contiguous_mass_fractions_arrays.copyFrom(contiguous_mass_fractions_arrays);
	}
	auto device_dt_heat = device.malloc<f64>(len);
	auto device_dt_contiguous_mass_fractions_arrays = device.malloc<f64>(species_len*len);

	production_rates(reaction_kinetics, len,
																device_pressure, device_temperature, device_contiguous_mass_fractions_arrays,
																/*&mut*/ device_dt_heat, /*&mut*/ device_dt_contiguous_mass_fractions_arrays);

	f64 dt_heat[len];
	device_dt_heat.copyTo(dt_heat);

	f64 dt_mass_fractions[species_len][len];
	f64 dt_contiguous_mass_fractions_arrays[species_len*len];
	device_dt_contiguous_mass_fractions_arrays.copyTo(dt_contiguous_mass_fractions_arrays);
	for(usize i=0; i<species_len; i+=1) for(usize id=0; id<len; id+=1) dt_mass_fractions[i][id] = dt_contiguous_mass_fractions_arrays[i*len+id];

	for(usize i=0; i<species_len-1; i++) printf("%e ", dt_mass_fractions[i][0]);
	return 0;
}
