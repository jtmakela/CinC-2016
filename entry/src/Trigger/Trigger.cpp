/* 
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
 * Trigger.cpp
 *
 *  Created on: Aug 24, 2015
 *      Author: jtmakela
 */

#include "Trigger.h"

namespace Accbpm {

/*
 * A rough energy slope trigger. Derived from the Accbpm project.
 *
 * @param cl_float const * Pointer to energy signal
 * @param size_t const & Length of energy signal in samples
 * @param doube const & Sample rate in time domain (eq. samples per second)
 */
Trigger::Trigger(cl_float const *energy, size_t const &len,
		double const &sample_freq) :
		sample_freq(sample_freq), energy(energy), len(len) {
	skip_offset = 0.5 * sample_freq; // There is often some noise in the beginning of the data we want to skip..
	define_threshold_value(energy, len, skip_offset, sample_freq, &limit);
	ev.clear();
	trig_do();
}

Trigger::~Trigger() {
}

/*
 * Fast implementation for finding kth smallest element in array.
 *
 * @param cl_float [] An array of values to sort. Notice that this implementation modifies the array.
 * @param ssize_t const Length of array
 * @param ssize_t const kth smallest value to extract from the array
 * @return double Value fo the kth smallest item in array
 */
double Trigger::kth_smallest(cl_float a[], ssize_t const n, ssize_t const k) {
	cl_float const &x = a[k];
	for (ssize_t l = 0, m = n - 1; l < m;) {
		ssize_t i = l;
		ssize_t j = m;
		do {
			while (a[i] < x) {
				++i;
			}
			while (x < a[j]) {
				--j;
			}
			if (i <= j) {
				double tmp = a[i]; // swap
				a[i++] = a[j];
				a[j--] = tmp;
			}
		} while (i <= j);

		if (j < k) {
			l = i;
		}
		if (k < i) {
			m = j;
		}
	}
	return x;
}

/*
 * Define trigger threshold value by energy minmaxminmax..
 *
 * A bit too complex threshold function for simple algorithm -- a copy from ECG trig algorithm
 *
 * TODO more suitable
 *
 * @param cl_float const * Pointer to the source signal
 * @param int Lenght of source signal in samples
 * @param int Number of samples to skip from the beginning of source signal
 * @param double Sample rate of the source signal in time domain (eq. samples per second)
 * @param double *threshold Pointer to output
 * @return int Always returns 0
 */
int Trigger::define_threshold_value(cl_float const *energy, int len,
		int skip_offset, double freq, double *threshold) {
	int base_len = freq * 0.100;
	int max_rr = freq * 3.0;
	int n_step = (len - skip_offset - base_len) / max_rr;
	cl_float *max, *base, base_max;
	ssize_t n, i, j, set;
	cl_float peak_estimate, base_estimate;

	if (n_step < 1) {
		return (-1); // too short data
	}
	max = static_cast<cl_float *>(mm_malloc(n_step * sizeof(cl_float)));
	base = static_cast<cl_float *>(mm_malloc(n_step * sizeof(cl_float)));

	for (n = 0; n < n_step; n++) {
		set = skip_offset + n * max_rr;
		max[n] = energy[set];
		for (i = set; i < set + max_rr; i++) {
			if (energy[i] > max[n]) {
				max[n] = energy[i];
			}
			base_max = energy[i];
			for (j = i + 1; j < i + base_len && j < len; j++) {
				if (energy[j] > base_max) {
					base_max = energy[j];
				}
			}
			if (i == set || base_max < base[n]) {
				base[n] = base_max;
			}
		}
	}
	peak_estimate = kth_smallest(max, n_step, 0.1 * n_step);
	base_estimate = kth_smallest(base, n_step, 0.9 * n_step);
	*threshold = base_estimate + 0.125 * (peak_estimate - base_estimate);

	mm_free(max);
	mm_free(base);
	return (0);
}

/*
 * The actual trigger
 *
 * Trig on energy signal above limit. Expect the signal to remain above limit for 400ms
 * tolerating 60ms continuous sinking below limit. Skips 200ms after each trig.
 *
 * @return int Always returns zero
 */
int Trigger::trig_do() {
	int const max_tolerance = 0.060 * sample_freq, max_above_len = 0.400
			* sample_freq;
	int const dead_time = 0.200 * sample_freq; // 0.250
	cl_float max;
	size_t max_i;
	size_t below_limit;

	cl_float low_limit_factor = 0.1;
	cl_float low_limit = 0.0;

	cl_float min = FLT_MAX;
	size_t k = 0;

	ev.reserve(512);

	for (size_t i = skip_offset, _len = len - skip_offset; i < _len; i++) {
		if (energy[i] < min) {
			min = energy[i];
		}

		// has energy been below the low limit
		if (energy[i] < low_limit) {
			continue;
		}

		if (energy[i] < energy[i - 1]) {
			continue;
		}

		if (energy[i] - min < 0.2 * limit) {
			continue;
		}

		if (i < k) {
			continue;
		}

		if (energy[i] > limit) {
			max_i = i;
			max = energy[i];
			below_limit = 0;
			for (size_t j = i + 1; j < len && j < i + max_above_len; j++) {
				if (energy[j] < limit) {
					below_limit++;
					if (below_limit > max_tolerance) {
						break;
					}
				} else {
					below_limit = 0;
					if (energy[j] > max) {
						max = energy[j];
						max_i = j;
					}
				}
			}

			min = max;

			struct ref_event e = { max_i };
			ev.push_back(e);

			k = max_i + dead_time;
			low_limit = low_limit_factor * max;
		}
	}
	return (0);
}

/*
 * Getter for trigged events
 *
 * @return std::vector<struct ref_event> const & Reference to vector of trigged events
 */
ref_ev const &Trigger::get_events() const {
	return ev;
}

} /* namespace Accbpm */
