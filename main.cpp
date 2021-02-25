//#!/bin/sh -c c++ main.cpp -lnekrk -locca -o main && ./main
#include <cstddef>
using usize = size_t ;
using f64 = double;
#define Vec vector
#include <nekrk.h>

auto main(int argc, const char **argv) -> int {
	using namespace ReactionKinetics;
	occa::device device((string)"{mode: 'Serial', 'kernel/compiler': 'clang', 'kernel/compiler_language': 'C++'}");
	auto reaction_kinetics = compile("Fuego", "grimech30", device, "", Reference{1., 1., 1., 1., 1.});
	f64 pressure = 101325.;
	f64 temperature = 1000.;
	const usize S = 53;
	f64 amounts[S] = {};
	/// !!! WARNING !!! C introduces undefined behavior by allowing to use uninitialized arrays and not zero initializing by default !!!
	{unsigned int i=0; while(i<S) { amounts[i] = 0.; i+=1; }}
	amounts[/*O2*/3] = 2./5.;
	amounts[/*CH4*/13] = 1./5.;
	amounts[/*Ar*/48] = 2./5.;

	//const usize N = 1;
	const usize N = 100000;

	Vec<f64> pressures(N);
	for(usize i=0; i<N; i+=1) pressures[i] = pressure;

	Vec<f64> temperatures(N);
	for(usize i=0; i<N; i+=1) temperatures[i] = temperature;

	Vec<f64> amounts_arrays[S]; for(usize i=0; i<S; i++) amounts_arrays[i] = Vec<f64>(N);
	for(usize i=0; i<S; i++) for(usize id=0; id<N; id+=1) amounts_arrays[i][id] = amounts[i];

	const usize len = N;
  auto device_pressure = device.malloc<f64>(len);
	device_pressure.copyFrom(pressures.data());
	auto device_temperature = device.malloc<f64>(len);
	device_temperature.copyFrom(temperatures.data());
	auto device_contiguous_amounts_arrays = device.malloc<f64>(S*len);
	{
		auto contiguous_amounts_arrays = Vec<f64>(S*len);
		for(usize i=0; i<S; i+=1) for(int id=0; id<len; id+=1) contiguous_amounts_arrays[i*len+id] = amounts_arrays[i][id];
		device_contiguous_amounts_arrays.copyFrom(contiguous_amounts_arrays.data());
	}
	auto device_dt_temperature = device.malloc<f64>(len);
	auto device_dt_contiguous_amounts_arrays = device.malloc<f64>(S*len);

	production_rates(reaction_kinetics, len,
																device_pressure, device_temperature, device_contiguous_amounts_arrays,
																/*&mut*/ device_dt_temperature, /*&mut*/ device_dt_contiguous_amounts_arrays);

	auto dt_temperature = Vec<f64>(len);
	device_dt_temperature.copyTo(dt_temperature.data());

	Vec<f64> dt_amounts[S]; for(usize i=0; i<S; i++) dt_amounts[i] = Vec<f64>(len);
	auto dt_contiguous_amounts_arrays = Vec<f64>(S*len);
	device_dt_contiguous_amounts_arrays.copyTo(dt_contiguous_amounts_arrays.data());
	for(usize i=0; i<S; i+=1) for(usize id=0; id<len; id+=1) dt_amounts[i][id] = dt_contiguous_amounts_arrays[i*len+id];

	for(usize i=0; i<S-1; i++) printf("%e ", dt_amounts[i][0]);
	return 0;
}
