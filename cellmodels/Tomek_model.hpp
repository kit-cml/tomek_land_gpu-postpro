#ifndef TOMEK_MODEL_ENDO_HPP
#define TOMEK_MODEL_ENDO_HPP

#include "cellmodel.hpp"
#include "enums/enum_Tomek_model.hpp"

	__device__ void initConsts(double *CONSTANTS, double *STATES, double type, double conc, double *ic50, double *cvar,  bool is_cvar, double bcl, double epsilon, int sample_id);
	__device__ void computeRates(double TIME, double* CONSTANTS, double* RATES, double* STATES, double* ALGEBRAIC, int sample_id,  double land_trpn);
	__device__ void solveAnalytical(double *CONSTANTS, double *STATES, double *ALGEBRAIC, double *RATES, double dt, int sample_id);
	__device__ double set_time_step(double TIME, double time_point, double max_time_step, double* CONSTANTS, double* RATES, int sample_id);
    __device__ void applyDrugEffect(double *CONSTANTS, double conc, double *hill, int sample_id);
	// __device__ void ___applyCvar(double *CONSTANTS, double *cvar, int sample_id);
	__device__ void ___initConsts(double *CONSTANTS, double *STATES, double type, double bcl, int sample_id);
	__device__ void ___gaussElimination(double *A, double *b, double *x, int N);

#endif