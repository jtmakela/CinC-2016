/**
 *  tftrig_final - a heart sound classifier
 *  Copyright (C) 2016 Jarno Mäkelä and Heikki Väänänen, RemoteA Ltd
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Iir.cpp
 *
 *  Created on: Apr 7, 2016
 *      Author: jtmakela
 */

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <math.h>
#include "macro.h"
#include "../utils/memory_manager.h"

#include "Iir.h"

#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif

namespace DSP {

Iir::Iir() {
}

Iir::~Iir() {
}

/**
 * Initialization of precalculated coefficients.
 */
void Iir::init_coefficients(enum iir_mode_e const mode) {
	coeff.a.clear();
	coeff.b.clear();
	if (mode <= LOW_PASS_2_POLE_0_45F
			|| (mode >= HIGH_PASS_2_POLE_0_01F && mode <= HIGH_PASS_2_POLE_0_45F)) {
		coeff.a.resize(3);
		coeff.b.resize(3);
	} else {
		coeff.a.resize(5);
		coeff.b.resize(5);
	}

	switch (mode) {
	case LOW_PASS_2_POLE_0_01F:
		coeff.a[0] = 8.663387E-04;
		coeff.a[1] = 1.732678E-03;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 1.919129E+00;
		coeff.b[2] = -9.225943E-01;
		break;
	case LOW_PASS_2_POLE_0_025F:
		coeff.a[0] = 5.112374E-03;
		coeff.a[1] = 1.022475E-02;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 1.797154E+00;
		coeff.b[2] = -8.176033E-01;
		break;
	case LOW_PASS_2_POLE_0_05F:
		coeff.a[0] = 1.868823E-02;
		coeff.a[1] = 3.737647E-02;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 1.593937E+00;
		coeff.b[2] = -6.686903E-01;
		break;
	case LOW_PASS_2_POLE_0_075F:
		coeff.a[0] = 3.869430E-02;
		coeff.a[1] = 7.738860E-01;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 1.392667E+00;
		coeff.b[2] = -5.474446E-01;
		break;
	case LOW_PASS_2_POLE_0_1F:
		coeff.a[0] = 6.372802E-02;
		coeff.a[1] = 1.274560E-01;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 1.194365E+00;
		coeff.b[2] = -4.492774E-01;
		break;
	case LOW_PASS_2_POLE_0_15F:
		coeff.a[0] = 1.254285E-01;
		coeff.a[1] = 2.508570E-01;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 8.070778E-01;
		coeff.b[2] = -3.087918E-01;
		break;
	case LOW_PASS_2_POLE_0_2F:
		coeff.a[0] = 1.997396E-01;
		coeff.a[1] = 3.994792E-01;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 4.291048E-01;
		coeff.b[2] = -2.280633E-01;
		break;
	case LOW_PASS_2_POLE_0_25F:
		coeff.a[0] = 2.858110E-01;
		coeff.a[1] = 5.716221E-01;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 5.423258E-02;
		coeff.b[2] = -1.974768E-01;
		break;
	case LOW_PASS_2_POLE_0_3F:
		coeff.a[0] = 3.849163E-01;
		coeff.a[1] = 7.698326E-01;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = -3.249116E-01;
		coeff.b[2] = -2.147536E-01;
		break;
	case LOW_PASS_2_POLE_0_35F:
		coeff.a[0] = 5.001024E-01;
		coeff.a[1] = 1.000205E+00;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = -7.158993E-01;
		coeff.b[2] = -2.845103E-01;
		break;
	case LOW_PASS_2_POLE_0_40F:
		coeff.a[0] = 6.362308E-01;
		coeff.a[1] = 1.272462E+00;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = -1.125379E+00;
		coeff.b[2] = -4.195441E-01;
		break;
	case LOW_PASS_2_POLE_0_45F:
		coeff.a[0] = 8.001101E-01;
		coeff.a[1] = 1.600220E+00;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = -1.556269E+00;
		coeff.b[2] = -6.441713E-01;
		break;

	case LOW_PASS_4_POLE_0_01F:
#ifdef VERBOSE
		fprintf(stdout, "*** unstable filter ***\n");
#endif
		coeff.a[0] = 4.149425E-07;
		coeff.a[1] = 1.659770E-06;
		coeff.a[2] = 2.489655E-06;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 3.893453E+00;
		coeff.b[2] = -5.688233E+00;
		coeff.b[3] = 3.695783E+00;
		coeff.b[4] = -9.010106E-01;
		break;
	case LOW_PASS_4_POLE_0_025F:
		coeff.a[0] = 1.504626E-05;
		coeff.a[1] = 6.018503E-05;
		coeff.a[2] = 9.027754E-05;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 3.725385E+00;
		coeff.b[2] = -5.226004E+00;
		coeff.b[3] = 3.270902E+00;
		coeff.b[4] = -7.705239E-01;
		break;
	case LOW_PASS_4_POLE_0_05F:
		coeff.a[0] = 2.141509E-04;
		coeff.a[1] = 8.566037E-04;
		coeff.a[2] = 1.284906E-03;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 3.425455E+00;
		coeff.b[2] = -4.479272E+00;
		coeff.b[3] = 2.643718E+00;
		coeff.b[4] = -5.933269E-01;
		break;
	case LOW_PASS_4_POLE_0_075F:
		coeff.a[0] = 9.726342E-04;
		coeff.a[1] = 3.890530E-03;
		coeff.a[2] = 5.835806E-03;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 3.103944E+00;
		coeff.b[2] = -3.774453E+00;
		coeff.b[3] = 2.111238E+00;
		coeff.b[4] = -4.562908E-01;
		break;
	case LOW_PASS_4_POLE_0_1F:
		coeff.a[0] = 2.780755E-03;
		coeff.a[1] = 1.112302E-02;
		coeff.a[2] = 1.668453E-02;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 2.764031E+00;
		coeff.b[2] = -3.122854E+00;
		coeff.b[3] = 1.664554E+00;
		coeff.b[4] = -3.502232E-01;
		break;
	case LOW_PASS_4_POLE_0_15F:
		coeff.a[0] = 1.180009E-02;
		coeff.a[1] = 4.720034E-02;
		coeff.a[2] = 7.080051E-02;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 2.039039E+00;
		coeff.b[2] = -2.012961E+00;
		coeff.b[3] = 9.897915E-01;
		coeff.b[4] = -2.046700E-01;
		break;
	case LOW_PASS_4_POLE_0_2F:
		coeff.a[0] = 3.224554E-02;
		coeff.a[1] = 1.289821E-01;
		coeff.a[2] = 1.934732E-01;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 1.265912E+00;
		coeff.b[2] = -1.203878E+00;
		coeff.b[3] = 5.405908E-01;
		coeff.b[4] = -1.185538E-01;
		break;
	case LOW_PASS_4_POLE_0_25F:
		coeff.a[0] = 7.015301E-02;
		coeff.a[1] = 2.806120E-01;
		coeff.a[2] = 4.209180E-01;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 4.541481E-01;
		coeff.b[2] = -7.417536E-01;
		coeff.b[3] = 2.361222E-01;
		coeff.b[4] = -7.096476E-02;
		break;
	case LOW_PASS_4_POLE_0_3F:
		coeff.a[0] = 1.335566E-01;
		coeff.a[1] = 5.342263E-01;
		coeff.a[2] = 8.013394E-01;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = -3.904486E-01;
		coeff.b[2] = -6.784138E-01;
		coeff.b[3] = -1.412021E-02;
		coeff.b[4] = -5.392238E-02;
		break;
	case LOW_PASS_4_POLE_0_35F:
		coeff.a[0] = 2.340973E-01;
		coeff.a[1] = 9.363892E-01;
		coeff.a[2] = 1.404584E+00;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = -1.263672E+00;
		coeff.b[2] = -1.080487E+00;
		coeff.b[3] = -3.276296E-01;
		coeff.b[4] = -7.376791E-02;
		break;
	case LOW_PASS_4_POLE_0_40F:
		coeff.a[0] = 3.896966E-01;
		coeff.a[1] = 1.558787E+00;
		coeff.a[2] = 2.338180E+00;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = -2.161179E+00;
		coeff.b[2] = -2.033992E+00;
		coeff.b[3] = -8.789098E-01;
		coeff.b[4] = -1.610655E-01;
		break;
	case LOW_PASS_4_POLE_0_45F:
		coeff.a[0] = 6.291693E-01;
		coeff.a[1] = 2.516677E+00;
		coeff.a[2] = 3.775016E+00;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = -3.077062E+00;
		coeff.b[2] = -3.641323E+00;
		coeff.b[3] = -1.949229E+00;
		coeff.b[4] = -3.990945E-01;
		break;

	case HIGH_PASS_2_POLE_0_01F:
		coeff.a[0] = 9.567529E-01;
		coeff.a[1] = -1.913506E+00;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 1.911437E+00;
		coeff.b[2] = -9.115749E-01;
		break;
	case HIGH_PASS_2_POLE_0_025F:
		coeff.a[0] = 8.950355E-01;
		coeff.a[1] = -1.790071E+00;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 1.777932E+00;
		coeff.b[2] = -8.022106E-01;
		break;
	case HIGH_PASS_2_POLE_0_05F:
		coeff.a[0] = 8.001102E-01;
		coeff.a[1] = -1.600220E+00;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 1.556269E+00;
		coeff.b[2] = -6.441715E-01;
		break;
	case HIGH_PASS_2_POLE_0_075F:
		coeff.a[0] = 7.142028E-01;
		coeff.a[1] = -1.428406E+00;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 1.338264E+00;
		coeff.b[2] = -5.185469E-01;
		break;
	case HIGH_PASS_2_POLE_0_1F:
		coeff.a[0] = 6.362307E-01;
		coeff.a[1] = -1.272461E+00;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 1.125379E+00;
		coeff.b[2] = -4.195440E-01;
		break;
	case HIGH_PASS_2_POLE_0_15F:
		coeff.a[0] = 5.001024E-01;
		coeff.a[1] = -1.000205E+00;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 7.158993E-01;
		coeff.b[2] = -2.845103E-01;
		break;
	case HIGH_PASS_2_POLE_0_2F:
		coeff.a[0] = 3.849163E-01;
		coeff.a[1] = -7.698326E-01;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = 3.249116E-01;
		coeff.b[2] = -2.147536E-01;
		break;
	case HIGH_PASS_2_POLE_0_25F:
		coeff.a[0] = 2.858111E-01;
		coeff.a[1] = -5.716222E-01;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = -5.423243E-02;
		coeff.b[2] = -1.974768E-01;
		break;
	case HIGH_PASS_2_POLE_0_3F:
		coeff.a[0] = 1.997396E-01;
		coeff.a[1] = -3.994792E-01;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = -4.291049E-01;
		coeff.b[2] = -2.280633E-01;
		break;
	case HIGH_PASS_2_POLE_0_35F:
		coeff.a[0] = 1.254285E-01;
		coeff.a[1] = -2.508570E-01;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = -8.070777E-01;
		coeff.b[2] = -3.087918E-01;
		break;
	case HIGH_PASS_2_POLE_0_40F:
		coeff.a[0] = 6.372801E-02;
		coeff.a[1] = -1.274560E-01;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = -1.194365E+00;
		coeff.b[2] = -4.492774E-01;
		break;
	case HIGH_PASS_2_POLE_0_45F:
		coeff.a[0] = 1.868823E-02;
		coeff.a[1] = -3.737647E-02;
		coeff.a[2] = coeff.a[0];
		coeff.b[1] = -1.593937E+00;
		coeff.b[2] = -6.686903E-01;
		break;

	case HIGH_PASS_4_POLE_0_01F:
#ifdef VERBOSE
		fprintf(stdout, "*** unstable filter ***\n");
#endif
		coeff.a[0] = 9.121579E-01;
		coeff.a[1] = -3.648632E+00;
		coeff.a[2] = 5.472947E+00;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 3.815952E+00;
		coeff.b[2] = -5.465026E+00;
		coeff.b[3] = 3.481295E+00;
		coeff.b[4] = -8.322529E-01;
		break;
	case HIGH_PASS_4_POLE_0_025F:
		coeff.a[0] = 7.941874E-01;
		coeff.a[1] = -3.176750E+00;
		coeff.a[2] = 4.765125E+00;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 3.538919E+00;
		coeff.b[2] = -4.722213E+00;
		coeff.b[3] = 2.814036E+00;
		coeff.b[4] = -6.318300E-01;
		break;
	case HIGH_PASS_4_POLE_0_05F:
		coeff.a[0] = 6.291694E-01;
		coeff.a[1] = -2.516678E+00;
		coeff.a[2] = 3.775016E+00;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 3.077062E+00;
		coeff.b[2] = -3.641324E+00;
		coeff.b[3] = 1.949230E+00;
		coeff.b[4] = -3.990947E-01;
		break;
	case HIGH_PASS_4_POLE_0_075F:
		coeff.a[0] = 4.965350E-01;
		coeff.a[1] = -1.986140E+00;
		coeff.a[2] = 2.979210E+00;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 2.617304E+00;
		coeff.b[2] = -2.749252E+00;
		coeff.b[3] = 1.325548E+00;
		coeff.b[4] = -2.524546E-01;
		break;
	case HIGH_PASS_4_POLE_0_1F:
		coeff.a[0] = 3.896966E-01;
		coeff.a[1] = -1.558786E+00;
		coeff.a[2] = 2.338179E+00;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 2.161179E+00;
		coeff.b[2] = -2.033991E+00;
		coeff.b[3] = 8.789094E-01;
		coeff.b[4] = -1.610655E-01;
		break;
	case HIGH_PASS_4_POLE_0_15F:
		coeff.a[0] = 2.340973E-01;
		coeff.a[1] = -9.363892E-01;
		coeff.a[2] = 1.404584E+00;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 1.263672E+00;
		coeff.b[2] = -1.080487E+00;
		coeff.b[3] = 3.276296E-01;
		coeff.b[4] = -7.376791E-02;
		break;
	case HIGH_PASS_4_POLE_0_2F:
		coeff.a[0] = 1.335566E-01;
		coeff.a[1] = -5.342262E-01;
		coeff.a[2] = 8.013393E-01;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = 3.904484E-01;
		coeff.b[2] = -6.784138E-01;
		coeff.b[3] = 1.412016E-02;
		coeff.b[4] = -5.392238E-02;
		break;
	case HIGH_PASS_4_POLE_0_25F:
		coeff.a[0] = 7.015302E-02;
		coeff.a[1] = -2.806121E-01;
		coeff.a[2] = 4.209182E-01;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = -4.541478E-01;
		coeff.b[2] = -7.417535E-01;
		coeff.b[3] = -2.361221E-01;
		coeff.b[4] = -7.096475E-02;
		break;
	case HIGH_PASS_4_POLE_0_3F:
		coeff.a[0] = 3.224553E-02;
		coeff.a[1] = -1.289821E-01;
		coeff.a[2] = 1.934732E-01;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = -1.265912E+00;
		coeff.b[2] = -1.203878E+00;
		coeff.b[3] = -5.405908E-01;
		coeff.b[4] = -1.185538E-01;
		break;
	case HIGH_PASS_4_POLE_0_35F:
		coeff.a[0] = 1.180009E-02;
		coeff.a[1] = -4.720035E-02;
		coeff.a[2] = 7.080051E-02;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = -2.039039E+00;
		coeff.b[2] = -2.012961E+00;
		coeff.b[3] = -9.897915E-01;
		coeff.b[4] = -2.046700E-01;
		break;
	case HIGH_PASS_4_POLE_0_40F:
		coeff.a[0] = 2.780754E-03;
		coeff.a[1] = -1.112302E-02;
		coeff.a[2] = 1.668453E-02;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = -2.764031E+00;
		coeff.b[2] = -3.122854E+00;
		coeff.b[3] = -1.664554E+00;
		coeff.b[4] = -3.502233E-01;
		break;
	case HIGH_PASS_4_POLE_0_45F:
		coeff.a[0] = 2.141509E-04;
		coeff.a[1] = -8.566037E-04;
		coeff.a[2] = 1.284906E-03;
		coeff.a[3] = coeff.a[1];
		coeff.a[4] = coeff.a[0];
		coeff.b[1] = -3.425455E+00;
		coeff.b[2] = -4.479272E+00;
		coeff.b[3] = -2.643718E+00;
		coeff.b[4] = -5.933269E-01;
		break;
	}
}

/**
 * The actual IIR calculator
 *
 * @param Pointer to pointer to output memory. The memory is allocated internally. Preallocated memory is freed.
 * @param Pointer to input data
 * @param length of input in samples
 * @return libc errno
 */
int Iir::calc(float **out, float const *in, size_t const &len) {
	off_t const padding = coeff.a.size();

	float *d1 = static_cast<float *>(mm_malloc(len * sizeof(float)));
	if (d1 == NULL) {
		PERROR("malloc");
		return errno;
	}

#ifdef VERBOSE
	fprintf(stdout, "  Calculating frequency response..\n");
#endif

	// first pass
	for (off_t i = 0; i < padding; ++i) {
		register float _d = coeff.a[0] * in[i];
		for (off_t j = 1, ___len = coeff.a.size(); j < ___len; ++j) {
			if (i < j) {
				_d += coeff.a[j] * in[0];
			} else {
				_d += coeff.a[j] * in[i - j];
				_d += coeff.b[j] * d1[i - j];
			}
		}

		d1[i] = _d;
	}

	for (off_t i = padding; i < len; ++i) {
		register float _d = coeff.a[0] * in[i];
		for (off_t j = 1, ___len = coeff.a.size(); j < ___len; ++j) {
			_d += coeff.a[j] * in[i - j];
			_d += coeff.b[j] * d1[i - j];
		}
		d1[i] = _d;
	}

	// second pass
	float *d2 = static_cast<float *>(mm_malloc(len * sizeof(float)));
	if (d2 == NULL) {
		PERROR("malloc");
		return errno;
	}

	for (off_t i = len - 1; i >= len - 1 - padding; --i) {
		register float _d = coeff.a[0] * d1[i];
		for (off_t j = 1, ___len = coeff.a.size(); j < ___len; ++j) {
			if (i + j >= len) {
				_d += coeff.a[j] * d1[len - 1];
			} else {
				_d += coeff.a[j] * d1[i + j];
				_d += coeff.b[j] * d2[i + j];
			}
		}

		d2[i] = _d;
	}

	for (off_t i = len - 1 - padding; i >= 0; --i) {
		register float _d = coeff.a[0] * d1[i];
		for (off_t j = 1, ___len = coeff.a.size(); j < ___len; ++j) {
			_d += coeff.a[j] * d1[i + j];
			_d += coeff.b[j] * d2[i + j];
		}
		d2[i] = _d;
	}
	mm_free(d1);

	if (*out) {
		mm_free(*out);
	}

	*out = d2;

	return 0;
}

/**
 * Run predefined filter
 *
 * @param Pointer to pointer to output memory. The memory is allocated internally. Preallocated memory is freed.
 * @param Pointer to input data
 * @param length of input in samples
 * @param filter to use
 * @return libc errno
 */
int Iir::filter(float **out, float const *in, size_t const &len,
		enum iir_mode_e const &mode) {
#ifdef VERBOSE
	fprintf(stdout, "  Initializing coefficients..\n");
#endif
	init_coefficients(mode);
	return calc(out, in, len);
}

/**
 * First order low-pass filter
 *
 * @param Pointer to pointer to output memory. The memory is allocated internally. Preallocated memory is freed.
 * @param Pointer to input data
 * @param length of input in samples
 * @param Cutoff frequency
 * @param Sample frequency
 * @return libc errno
 */
int Iir::low_pass(float **out, const float *in, const size_t len,
		const float cutoff_freq, const float sample_freq) {
	calc_simple_coefficients(cutoff_freq, sample_freq);
	return calc(out, in, len);
}

/**
 * Nth order Chebyshev low-pass filter
 *
 * @param Pointer to pointer to output memory. The memory is allocated internally. Preallocated memory is freed.
 * @param Pointer to input data
 * @param length of input in samples
 * @param Cutoff frequency
 * @param Allowed ripple percentage
 * @parma Number of poles
 * @param Sample frequency
 * @return libc errno
 */
int Iir::low_pass(float **out, const float *in, const size_t len,
		const float cutoff_freq, float const ripple_percent,
		const size_t number_of_poles, const float sample_freq) {
	calc_chebyshev_coefficients(cutoff_freq, false, ripple_percent,
			number_of_poles, sample_freq);
	return calc(out, in, len);
}

/**
 * First order high-pass filter
 *
 * @param Pointer to pointer to output memory. The memory is allocated internally. Preallocated memory is freed.
 * @param Pointer to input data
 * @param length of input in samples
 * @param Cutoff frequency
 * @param Sample frequency
 * @return libc errno
 */
int Iir::high_pass(float **out, const float *in, const size_t len,
		const float cutoff_freq, const float sample_freq) {
	calc_simple_coefficients(cutoff_freq, sample_freq);

	float *d1 = 0;
	calc(&d1, in, len);

	float *d2 = static_cast<float *>(mm_malloc(len * sizeof(float)));
	if (d2 == NULL) {
		PERROR("malloc");
		return errno;
	}

	for (size_t i = 0; i < len; ++i) {
		d2[i] = in[i] - d1[i];
	}

	mm_free(d1);

	if (*out) {
		mm_free(*out);
	}

	*out = d2;

	return 0;
}

/**
 * Nth order Chebyshev high-pass filter
 *
 * @param Pointer to pointer to output memory. The memory is allocated internally. Preallocated memory is freed.
 * @param Pointer to input data
 * @param length of input in samples
 * @param Cutoff frequency
 * @param Allowed ripple percentage
 * @parma Number of poles
 * @param Sample frequency
 * @return libc errno
 */
int Iir::high_pass(float **out, const float *in, const size_t len,
		const float cutoff_freq, float const ripple_percent,
		const size_t number_of_poles, const float sample_freq) {
	calc_chebyshev_coefficients(cutoff_freq, true, ripple_percent,
			number_of_poles, sample_freq);
	return calc(out, in, len);
}

/**
 * Narrow pass filter
 *
 * @param Pointer to pointer to output memory. The memory is allocated internally. Preallocated memory is freed.
 * @param Pointer to input data
 * @param length of input in samples
 * @param Center frequency
 * @param Passband around the center frequency
 * @param Sample frequency
 * @return libc errno
 */
int Iir::narrow_pass(float **out, const float *in, const size_t len,
		const float center_freq, const float bandwidth,
		const float sample_freq) {

	calc_narrow_pass_coefficients(bandwidth, center_freq);
	return calc(out, in, len);
}

/**
 * First order bandpass filter
 *
 * @param Pointer to pointer to output memory. The memory is allocated internally. Preallocated memory is freed.
 * @param Pointer to input data
 * @param length of input in samples
 * @param Cutoff low frequency
 * @param Cutoff high frequency
 * @param Sample frequency
 * @return libc errno
 */
int Iir::bandpass(float **out, const float *in, const size_t len,
		const float low_freq, const float high_freq, const float sample_freq) {

	if (low_freq == 0) {
		if (high_freq > 0.5 * sample_freq) {
			float *d = static_cast<float *>(mm_malloc(len * sizeof(float)));
			if (d == NULL) {
				PERROR("malloc");
				return errno;
			}

			memcpy(d, in, len * sizeof(float));

			if (*out) {
				mm_free(*out);
			}

			*out = d;
		}
		low_pass(out, in, len, high_freq, sample_freq);
	}

	if (high_freq > 0.5 * sample_freq) {
		high_pass(out, in, len, low_freq, sample_freq);
		return 0;
	}

	float *d = 0;
	low_pass(&d, in, len, high_freq, sample_freq);
	high_pass(out, d, len, low_freq, sample_freq);

	mm_free(d);

	return 0;
}

/**
 * Nth order Chebyshev bandpass filter
 *
 * @param Pointer to pointer to output memory. The memory is allocated internally. Preallocated memory is freed.
 * @param Pointer to input data
 * @param length of input in samples
 * @param Cutoff low frequency
 * @param Cutoff high frequency
 * @param Allowed ripple percentage
 * @parma Number of poles
 * @param Sample frequency
 * @return libc errno
 */
int Iir::bandpass(float **out, const float *in, const size_t len,
		const float low_freq, const float high_freq, float const ripple_percent,
		const size_t number_of_poles, const float sample_freq) {

	if (low_freq == 0) {
		if (high_freq > 0.5 * sample_freq) {
			float *d = static_cast<float *>(mm_malloc(len * sizeof(float)));
			if (d == NULL) {
				PERROR("malloc");
				return errno;
			}

			memcpy(d, in, len * sizeof(float));

			if (*out) {
				mm_free(*out);
			}

			*out = d;
		}

		low_pass(out, in, len, high_freq, ripple_percent, number_of_poles,
				sample_freq);
		return 0;
	}

	if (high_freq > 0.5 * sample_freq) {
		high_pass(out, in, len, low_freq, ripple_percent, number_of_poles,
				sample_freq);
		return 0;
	}

	float *d = 0;
	low_pass(&d, in, len, high_freq, ripple_percent, number_of_poles,
			sample_freq);
	high_pass(out, d, len, low_freq, ripple_percent, number_of_poles,
			sample_freq);
	mm_free(d);

	return 0;
}

/**
 * Nth order Chebyshev bandpass filter
 *
 * @param Pointer to pointer to output memory. The memory is allocated internally. Preallocated memory is freed.
 * @param Pointer to input data
 * @param length of input in samples
 * @param Cutoff low frequency
 * @param Cutoff high frequency
 * @param Allowed ripple percentage
 * @parma Number of poles
 * @param Sample frequency
 * @return libc errno
 */
int Iir::bandpass(double **out, const double *in, const size_t len,
		const float low_freq, const float high_freq, float const ripple_percent,
		const size_t number_of_poles, const float sample_freq) {

	float *tmp = static_cast<float *>(mm_malloc(len * sizeof(float)));
	for (size_t i = 0; i < len; ++i) {
		tmp[i] = in[i];
	}
	float *_out = 0;

	bandpass(&_out, tmp, len, low_freq, high_freq, ripple_percent,
			number_of_poles, sample_freq);

	mm_free(tmp);

	if (*out) {
		mm_free(*out);
	}

	*out = static_cast<double *>(mm_malloc(len * sizeof(double)));
	for (size_t i = 0; i < len; ++i) {
		(*out)[i] = _out[i];
	}

	mm_free(_out);

	return 0;
}

/**
 * Narrow-pass coefficient calculator
 *
 * @param Bandwidth
 * @param Center frequency
 */
void Iir::calc_narrow_pass_coefficients(const float bandwidth,
		const float center_freq) {
#ifdef VERBOSE
	fprintf(stdout, "  Calculating coefficients..\n");
#endif

	float cosw = cos(2 * M_PI * center_freq);

	float R = 1 - 3 * bandwidth;
	float K = (1 - 2 * R * cosw + R * R) / (2 - 2 * cosw);

	coeff.a.resize(3);
	coeff.b.resize(3);

	coeff.a[0] = 1 - K;
	coeff.a[1] = 2 * (K - R) * cosw;
	coeff.a[2] = R * R - K;
	coeff.b[1] = 2 * R * cosw;
	coeff.b[2] = -R * R;
}

/**
 * First order coefficient calculator
 *
 * @param Cutoff frequency
 * @param Sample frequency
 */
void Iir::calc_simple_coefficients(float const cutoff_freq,
		float const sample_freq) {
	float f = cutoff_freq / sample_freq;
	float x = exp(-14.445 * f);

	coeff.a.resize(5);
	coeff.b.resize(5);

	coeff.a[0] = pow(1 - x, 4);
	coeff.a[1] = 0;
	coeff.a[2] = 0;
	coeff.a[3] = 0;
	coeff.a[4] = 0;

	coeff.b[1] = 4.0 * x;
	coeff.b[2] = -6.0 * pow(x, 2);
	coeff.b[3] = 4.0 * pow(x, 3);
	coeff.b[4] = -pow(x, 4);
}

/**
 *
 * @param Cut-off frequency. Must be in range 0 to 0.5 times the sampling frequency.
 * @param Are the desired coefficients for high-pass (or low pass) filter
 * @param Desired ripple percentage in range 0 to 29
 * @param Number of poles, an event integer in berween 2 and 20
 * @oaram Desired sampling frequency domain
 */
void Iir::calc_chebyshev_coefficients(float const cutoff_freq,
		bool const is_high_pass, float const ripple_percent,
		const int number_of_poles, const float sample_freq) {
	float const _cutoff_freq = cutoff_freq / sample_freq;

	char cachefilename[FILENAME_MAX];
	sprintf(cachefilename, "/tmp/Iir.chebyshev_coefficients.%f_%u_%f_%i.dat",
			_cutoff_freq, is_high_pass, ripple_percent, number_of_poles);
	{
		struct stat sb;
		if (stat(cachefilename, &sb) == 0) {
			FILE *f = fopen(cachefilename, "r");
			if (f == 0) {
				PERROR("fopen");
				throw errno;
			}

			int n = 0;
			if (fread(&n, sizeof(int), 1, f) != 1) {
				PERROR("fread");
				throw errno;
			}

			coeff.a.resize(n);
			coeff.b.resize(n);

			for (std::vector<float>::iterator it = coeff.a.begin();
					it != coeff.a.end(); ++it) {
				if (fread(&(*it), sizeof(float), 1, f) != 1) {
					PERROR("fread");
					throw errno;
				}
			}

			for (std::vector<float>::iterator it = coeff.b.begin();
					it != coeff.b.end(); ++it) {
				if (fread(&(*it), sizeof(float), 1, f) != 1) {
					PERROR("fread");
					throw errno;
				}
			}

			fclose(f);

			return;
		}
	}

	int n = 1 + number_of_poles + 2;

	coeff.a.clear();
	coeff.b.clear();

	coeff.a.resize(n);
	coeff.b.resize(n);

	coeff.a[2] = 1, coeff.b[2] = 1;

	{
		struct coefficients tmp;
		tmp.a.resize(n);
		tmp.b.resize(n);

		for (size_t i = 0, _len = number_of_poles / 2; i < _len; ++i) {
			struct coefficients tmp2;
			tmp2.a.clear();
			tmp2.b.clear();

			chebyshev_coefficient_iterator(_cutoff_freq, is_high_pass,
					ripple_percent, number_of_poles, i, tmp2);

			for (int i = 0; i < n; ++i) {
				tmp.a[i] = coeff.a[i], tmp.b[i] = coeff.b[i];
			}

			for (int i = 2; i < n; ++i) {
				coeff.a[i] = tmp2.a[0] * tmp.a[i] + tmp2.a[1] * tmp.a[i - 1]
						+ tmp2.a[2] * tmp.a[i - 2];
				coeff.b[i] = tmp.b[i] - tmp2.b[1] * tmp.b[i - 1]
						- tmp2.b[2] * tmp.b[i - 2];
			}
		}
	}

	coeff.b[2] = 0;

	n = 1 + number_of_poles;

	for (int i = 0; i < n; ++i) {
		coeff.a[i] = coeff.a[i + 2];
		coeff.b[i] = -coeff.b[i + 2];
	}

	// normalize gain
	float sa = 0, sb = 0;
	if (is_high_pass) {
		for (int i = 0; i < n; ++i) {
			sa += coeff.a[i] * pow(-1, i);
			sb += coeff.b[i] * pow(-1, i);
		}
	} else {
		for (int i = 0; i < n; ++i) {
			sa += coeff.a[i];
			sb += coeff.b[i];
		}
	}

	float gain = sa / (1 - sb);
	for (int i = 0; i < n; ++i) {
		coeff.a[i] /= gain;
	}

	coeff.a.resize(n);
	coeff.b.resize(n);

	FILE *f = fopen(cachefilename, "w");
	if (f == 0) {
		PERROR("fopen");
		throw errno;
	}

	if (fwrite(&n, sizeof(int), 1, f) != 1) {
		PERROR("fwrite");
		throw errno;
	}

	for (std::vector<float>::const_iterator it = coeff.a.begin();
			it != coeff.a.end(); ++it) {
		if (fwrite(&(*it), sizeof(float), 1, f) != 1) {
			PERROR("fwrite");
			throw errno;
		}
	}

	for (std::vector<float>::const_iterator it = coeff.b.begin();
			it != coeff.b.end(); ++it) {
		if (fwrite(&(*it), sizeof(float), 1, f) != 1) {
			PERROR("fwrite");
			throw errno;
		}
	}

	fclose(f);

	return;
}

/**
 *
 * @param Cut-off frequency in fractions of sampling frequency. Must be in range 0 to 0.5
 * @param Are the desired coefficients for high-pass (or low pass) filter
 * @param Desired ripple percentage in range 0 to 29
 * @param Number of poles, an event integer in berween 2 and 20
 * @param Iteration round index
 */
void Iir::chebyshev_coefficient_iterator(float const cutoff_freq,
		bool const is_high_pass, float const ripple_percent,
		int const number_of_poles, int ii, struct coefficients &out) {

	// calculate the pole location on the unit circle
	float re = -cos(
			M_PI / (2.0 * number_of_poles)
					+ static_cast<float>(ii) * M_PI / number_of_poles);
	float im = sin(
			M_PI / (2.0 * number_of_poles)
					+ static_cast<float>(ii) * M_PI / number_of_poles);

	// wrap from a circle to an ellipse
	if (ripple_percent) {
		float es = sqrt(pow(100.0 / (100.0 - ripple_percent), 2) - 1.0);
		float vx = (1.0 / number_of_poles)
				* log((1.0 / es) + sqrt(1.0 / pow(es, 2) + 1));
		float kx = (1.0 / number_of_poles)
				* log((1.0 / es) + sqrt(1.0 / pow(es, 2) - 1));
		kx = (exp(kx) + exp(-kx)) / 2;

		re *= ((exp(vx) - exp(-vx)) / 2) / kx;
		im *= ((exp(vx) + exp(-vx)) / 2) / kx;
	}

	// s-domain to z-domain conversion
	float const t = 2 * tan(0.5);
	float const t_pow2 = pow(t, 2);
	float const w = 2 * M_PI * cutoff_freq;
	float const m = pow(re, 2) + pow(im, 2);
	float d = 4 - 4 * re * t + m * t_pow2;

	float const x0 = t_pow2 / d;
	float const x1 = 2 * t_pow2 / d;
	float const x2 = t_pow2 / d;
	float const y1 = (8.0 - 2.0 * m * t_pow2) / d;
	float const y2 = (-4.0 - 4.0 * re * t - m * t_pow2) / d;

	// lp to lp, or lp to hp transform
	float k =
			is_high_pass ?
					-cos(w / 2 + 0.5) / cos(w / 2 - 0.5) :
					sin(0.5 - w / 2) / sin(0.5 + w / 2);

	d = 1 + y1 * k - y2 * pow(k, 2);

	out.a.resize(3);
	out.b.resize(3);

	out.a[0] = (x0 - x1 * k + x2 * pow(k, 2)) / d;
	out.a[1] = (is_high_pass ? -1 : 1)
			* ((-2 * x0 * k + x1 + x1 * pow(k, 2) - 2 * x2 * k) / d);
	out.a[2] = (x0 * pow(k, 2) - x1 * k + x2) / d;
	out.b[1] = (is_high_pass ? -1 : 1)
			* ((2 * k + y1 + y1 * pow(k, 2) - 2 * y2 * k) / d);
	out.b[2] = (-pow(k, 2) - y1 * k + y2) / d;
}

} /* namespace DSP */
