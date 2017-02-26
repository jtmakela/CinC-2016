/*
 * Iir.h
 *
 *  Created on: Apr 7, 2016
 *      Author: jtmakela
 */

#ifndef SRC_MYDSP_IIR_H_
#define SRC_MYDSP_IIR_H_

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include <vector>
#include "../utils/memory_manager.h"

#define M_PI		3.14159265358979323846

namespace DSP {

enum iir_mode_e {
	LOW_PASS_2_POLE_0_01F,
	LOW_PASS_2_POLE_0_025F,
	LOW_PASS_2_POLE_0_05F,
	LOW_PASS_2_POLE_0_075F,
	LOW_PASS_2_POLE_0_1F,
	LOW_PASS_2_POLE_0_15F,
	LOW_PASS_2_POLE_0_2F,
	LOW_PASS_2_POLE_0_25F,
	LOW_PASS_2_POLE_0_3F,
	LOW_PASS_2_POLE_0_35F,
	LOW_PASS_2_POLE_0_40F,
	LOW_PASS_2_POLE_0_45F,

	LOW_PASS_4_POLE_0_01F,
	LOW_PASS_4_POLE_0_025F,
	LOW_PASS_4_POLE_0_05F,
	LOW_PASS_4_POLE_0_075F,
	LOW_PASS_4_POLE_0_1F,
	LOW_PASS_4_POLE_0_15F,
	LOW_PASS_4_POLE_0_2F,
	LOW_PASS_4_POLE_0_25F,
	LOW_PASS_4_POLE_0_3F,
	LOW_PASS_4_POLE_0_35F,
	LOW_PASS_4_POLE_0_40F,
	LOW_PASS_4_POLE_0_45F,

	HIGH_PASS_2_POLE_0_01F,
	HIGH_PASS_2_POLE_0_025F,
	HIGH_PASS_2_POLE_0_05F,
	HIGH_PASS_2_POLE_0_075F,
	HIGH_PASS_2_POLE_0_1F,
	HIGH_PASS_2_POLE_0_15F,
	HIGH_PASS_2_POLE_0_2F,
	HIGH_PASS_2_POLE_0_25F,
	HIGH_PASS_2_POLE_0_3F,
	HIGH_PASS_2_POLE_0_35F,
	HIGH_PASS_2_POLE_0_40F,
	HIGH_PASS_2_POLE_0_45F,

	HIGH_PASS_4_POLE_0_01F,
	HIGH_PASS_4_POLE_0_025F,
	HIGH_PASS_4_POLE_0_05F,
	HIGH_PASS_4_POLE_0_075F,
	HIGH_PASS_4_POLE_0_1F,
	HIGH_PASS_4_POLE_0_15F,
	HIGH_PASS_4_POLE_0_2F,
	HIGH_PASS_4_POLE_0_25F,
	HIGH_PASS_4_POLE_0_3F,
	HIGH_PASS_4_POLE_0_35F,
	HIGH_PASS_4_POLE_0_40F,
	HIGH_PASS_4_POLE_0_45F,
};

struct coefficients {
	std::vector<float> a;
	std::vector<float> b;
};

class Iir {
public:
	Iir();
	virtual ~Iir();

	void init_coefficients(enum iir_mode_e const mode);
	int calc(float **out, float const *in, size_t const &len);

	int filter(float **out, float const *in, size_t const &len,
			enum iir_mode_e const &mode);

	int low_pass(float **out, const float *in, const size_t len,
			const float cutoff_freq, const float sample_freq);

	int low_pass(float **out, const float *in, const size_t len,
			const float cutoff_freq, float const ripple_percent,
			const size_t number_of_poles, const float sample_freq);

	int high_pass(float **out, const float *in, const size_t len,
			const float cutoff_freq, const float sample_freq);

	int high_pass(float **out, const float *in, const size_t len,
			const float cutoff_freq, float const ripple_percent,
			const size_t number_of_poles, const float sample_freq);

	int narrow_pass(float **out, const float *in, const size_t len,
			const float center_freq, const float bandwidth,
			const float sample_freq);

	int bandpass(float **out, const float *in, const size_t len,
			const float low_freq, const float high_freq,
			const float sample_freq);

	int bandpass(float **out, const float *in, const size_t len,
			const float low_freq, const float high_freq,
			float const ripple_percent, const size_t number_of_poles,
			const float sample_freq);

	int bandpass(double **out, const double *in, const size_t len,
			const float low_freq, const float high_freq,
			float const ripple_percent, const size_t number_of_poles,
			const float sample_freq);

private:
	struct coefficients coeff;

	void calc_narrow_pass_coefficients(const float bandwidth,
			const float center_freq);
	void calc_simple_coefficients(float const cutoff_freq,
			float const sample_freq);
	void calc_chebyshev_coefficients(float const cutoff_freq,
			bool const is_high_pass, float const ripple_percent,
			const int number_of_poles, const float sample_freq);

	void chebyshev_coefficient_iterator(float const cutoff_freq,
			bool const is_high_pass, float const ripple_percent,
			int const number_of_poles, int ii, struct coefficients &out);
};

} /* namespace DSP */

#endif /* SRC_MYDSP_IIR_H_ */
