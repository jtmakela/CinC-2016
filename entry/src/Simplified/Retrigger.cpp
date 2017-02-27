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
 * Retrigger.cpp
 *
 *  Created on: Mar 27, 2016
 *      Author: jtmakela
 */

#include "Retrigger.h"

namespace Simplified {

/**
 * Setup routines for retigger
 *
 * Uses cached Blackman window (/tmp/Retrigger.blackman,window.dat). Must be manually deleted
 * if window langht or sample rate is altered.
 *
 * @param double const & Energy function window length in fraction of sample time (seconds)
 */
Retrigger::Retrigger(double const &window_length_in_fractions_of_sample_time) :
		raw(0), len(0), len_in_bytes(0), in(0), energy(0), blackman_window(0), retrig_convolution_kernel(
				0), evx_offset(0) {
	energy_window.len = window_length_in_fractions_of_sample_time * sample_freq;
	energy_window.offset = 0;
	energy_window.len_in_bytes = energy_window.len * sizeof(cl_float);

	convolution_window.len = 0;
	convolution_window.offset = 0;
	convolution_window.len_in_bytes = 0;

	lookaround_window.len = ceil(0.1 * sample_freq);

	lookaround_window.offset = lookaround_window.len / 2;
	lookaround_window.len_in_bytes = lookaround_window.len * sizeof(cl_float);

	correlation_window.len = 0.4 * sample_freq;
	correlation_window.offset = correlation_window.len / 2;
	correlation_window.len_in_bytes = correlation_window.len * sizeof(cl_float);

	blackman_window = static_cast<cl_float *>(mm_malloc(
			energy_window.len_in_bytes));

	{
		struct stat sb;
		char cachefilename[FILENAME_MAX];
		sprintf(cachefilename, "/tmp/Retrigger.blackman window.dat");
		if (stat(cachefilename, &sb) == 0) {
			FILE *f = fopen(cachefilename, "r");
			if (f == 0) {
				PERROR("fopen");
				throw errno;
			}

			if (fread(blackman_window, sizeof(cl_float), energy_window.len, f)
					!= energy_window.len) {
				PERROR("fread");
				throw errno;
			}

			fclose(f);
		} else {
			cl_float const twopiperm = 2 * M_PI / (energy_window.len - 1);
			for (size_t i = 0; i < energy_window.len; ++i) {
				// strict blackman
				blackman_window[i] = 0.42 - 0.5 * cos(i * twopiperm)
						+ 0.08 * cos(2 * i * twopiperm);
			}

			FILE *f = fopen(cachefilename, "w");
			if (f == 0) {
				PERROR("fopen");
				throw errno;
			}

			if (fwrite(blackman_window, sizeof(cl_float), energy_window.len, f)
					!= energy_window.len) {
				PERROR("fwrite");
				throw errno;
			}

			fclose(f);
		}
	}
}

Retrigger::~Retrigger() {
	if (blackman_window) {
		mm_free(blackman_window);
	}

	if (retrig_convolution_kernel) {
		mm_free(retrig_convolution_kernel);
	}

	mm_free(in);
	mm_free(energy);
}

/**
 * Set data and generates the energy norm
 *
 * Energy norm is defined by:
 * 1. IIR bandpass filter;
 * 2. FIR filter with externally set convolution kernel; and
 * 3. Smooth with Blackman window
 *
 *  @param data_raw_t const * Pointer to data
 *  @param size_t const & Length of data in samples
 */
void Retrigger::set_data(data_raw_t const *raw, size_t const &len) {
	this->raw = raw;
	this->len = len;

	if (in) {
		mm_free(in);
		in = 0;
	}

	DSP::Iir iir;
	float *tmp = static_cast<float *>(mm_malloc(len * sizeof(float)));
	for (size_t i = 0; i < len; ++i) {
		*(tmp + i) = *(raw + i);
	}

	iir.bandpass(&in, tmp, len, 10.0, 500.0, 0.5, 4, sample_freq);
	mm_free(tmp);

	calculate_energy();
}

/**
 * Set convolution kernel for retirgger routines
 *
 * Simplified to PhysioNet environment from OpenCL enabled environment.
 *
 * @param cl_float const * Pointer to data (Convolution kernel)
 * @param size_t const & Length of data in samples
 */
void Retrigger::set_convolution_kernel(cl_float const *kernel,
		size_t const &len) {
	convolution_window.len = len;
	convolution_window.offset = 0;
	convolution_window.len_in_bytes = sizeof(cl_float) * len;

	if (retrig_convolution_kernel) {
		mm_free(retrig_convolution_kernel);
	}

	retrig_convolution_kernel = reinterpret_cast<cl_float *>(mm_malloc(
			convolution_window.len_in_bytes));

	memcpy(retrig_convolution_kernel, kernel, convolution_window.len_in_bytes);
}

#define POW2(_a) ((_a) * (_a))

/**
 * A convolution routine.
 *
 * See retrig.cl
 *
 * @param cl_float * Pointer to output memory space. Must be properly allocated.
 * @param cl_float * Pointer to longer convolution component (stream)
 * @param ulong const Length of longer convolution component
 * @param cl_float * Pointer to shorter convolution component (window)
 * @param ulong const Length of shorter convolution component
 */
static void retrig_convolution_program(cl_float *out, cl_float *a,
		ulong const a_len, cl_float *b, ulong const b_len) {
	ulong const len = b_len / 2;

	memset(out, -0.0, len * sizeof(cl_float));
	memset(out + a_len - len - 1, -0.0, len * sizeof(cl_float));

	for (ulong id = len, _len = a_len - len; id < _len; ++id) {
		cl_float const *p = a + id - len;
		register cl_float d = 0.0;
		for (ulong i = 0; i < b_len; ++i) {
			d += *(p + i) * *(b + i);
		}
		out[id] = d / b_len;
	}
}

/**
 * An energy norm routine
 *
 * See retrig.cl
 *
 * @param cl_float * Pointer to output memory space. Must be properly allocated.
 * @param cl_float * Pointer to data (stream)
 * @param ulong const Length of data in samples
 * @param cl_float * Pointer to window function memory space
 * @param ulong const Length of window function in bytes
 */
static void retrig_energy_program(cl_float *out, cl_float *a, ulong const a_len,
		cl_float *b, ulong const b_len) {
	ulong const len = b_len / 2;

	memset(out, -0.0, len * sizeof(cl_float));
	memset(out + a_len - len - 1, -0.0, len * sizeof(cl_float));

	for (ulong id = len, _len = a_len - len; id < _len; ++id) {
		cl_float const *p = a + id - len;
		register cl_float d = 0.0;
		for (ulong i = 0; i < b_len; ++i) {
			d += POW2(*(p + i) * *(b + i));
		}
		out[id] = d;
	}
}

/**
 * Calculate the internal energy signal from internally set data signal
 */
void Retrigger::calculate_energy() {
	if (energy) {
		mm_free(energy);
	}

	cl_float *tmp = static_cast<cl_float *>(mm_calloc(len, sizeof(cl_float)));
	if (tmp == 0) {
		PERROR("calloc");
		throw errno;
	}

	// first pass
	// convolve with the convolution kernel
	retrig_convolution_program(tmp, in, len, retrig_convolution_kernel,
			convolution_window.len);

	energy = static_cast<cl_float *>(mm_calloc(len, sizeof(cl_float)));
	if (energy == 0) {
		PERROR("calloc");
		throw errno;
	}

	// second pass
	// calculate energy
	retrig_energy_program(energy, tmp, len, blackman_window, energy_window.len);

	mm_free(tmp);
}

/**
 * Find local energy maxima near reference events and set locally stored temporary events accordingly
 *
 * @param struct ref_event const & Reference event
 * @param off_t const & Bias in sample space
 */
void Retrigger::calc(struct ref_event const &e, off_t const &_offset) {
	size_t const offset = e.offset + _offset;

	if (offset < lookaround_window.offset
			|| offset - lookaround_window.offset + lookaround_window.len
					> len) {
		fprintf(stderr, "[%s:%u] partial window.\n", __FILE__,
		__LINE__);
		return;
	}

	size_t max_at = offset;
	cl_float max = energy[max_at];
	for (size_t i = offset - lookaround_window.offset, _len = i
			+ lookaround_window.len; i < _len; ++i) {
		cl_float const &d = energy[i];
		if (d > max) {
			max_at = i;
			max = d;
		}
	}

	struct ref_event_ext _e = { &e, max_at, static_cast<off_t>(max_at - offset),
			max, 0 };
	evx.push_back(_e);

#ifdef DEBUG
	printf("max at %lu, change: %li\n", max_at, (int64_t) max_at - offset);
#endif
}

/**
 * Calculate average
 *
 * @param cl_float * Pointer to data
 * @param ulong const & Length of data in samples
 * @return cl_float Average of len data points
 */
static cl_float avg(cl_float const *d, ulong const &len) {
	cl_float _avg = 0.0;
	for (ulong i = 0; i < len; ++i) {
		_avg += *d++;
	}
	return _avg / len;
}

/**
 * A custom correlation signal routine
 *
 * Divider is set as the greatest standard deviation. See the article for ruther reasoning.
 *
 * See retrig.cl
 *
 * @param cl_float * Pointer to output memory space. Must be properly allocated.
 * @param cl_float * Pointer to longer correlation component (stream)
 * @param ulong const Length of longer correlation component
 * @param cl_float * Pointer to shorter correlation component (window)
 * @param ulong const Length of shorter correlation component
 * @param cl_float const & Precalculated average value of the shorter correlation component (window)
 * @param cl_float const & Precalculated standard deviation value of the shorter correlation component (window)
 */
static void retrig_correlation_program(cl_float *out, cl_float *a,
		ulong const &a_len, cl_float *b, ulong const &b_len,
		cl_float const &avg_b, cl_float const &stdev_b) {
	ulong const len = b_len / 2;
	for (ulong id = 0; id < a_len; ++id) {
		cl_float const *p = a + id - len;

		cl_float const avg_a = avg(p, b_len);
		cl_float stdev_a = 0.0;

		cl_float conv = 0.0;
		for (ulong i = 0; i < b_len; ++i) {
			register cl_float _a = *(p + i) - avg_a;
			register cl_float _b = *(b + i) - avg_b;
			conv += _a * _b;
			stdev_a += POW2(_a);
		}

		out[id] = conv / (stdev_a > stdev_b ? stdev_a : stdev_b);
	}
}

struct range {
	off_t set;
	off_t end;
};

/**
 * Calculate cross correlations near each local temporary event
 */
void Retrigger::calculate_correlations() {
	std::vector<struct range> ranges;
	struct range r = { };
	for (ref_ev_ext_it it = evx.begin(); it != evx.end(); ++it) {
		off_t set = it->offset - correlation_window.offset;
		if (set < 0) {
			continue;
		}
		off_t end = it->offset + correlation_window.len;
		if (len < end) {
			continue;
		}

		if (end > r.end) {
			if (r.set) {
				r.end = end;
				ranges.push_back(r);
			}
			r.set = set;
			r.end = end;
		} else {
			r.end = end;
		}
	}

	for (ref_ev_ext_it it = evx.begin(); it != evx.end(); ++it) {
		it->correlation.signal = (cl_float *) mm_calloc(len, sizeof(cl_float));
		if (it->correlation.signal == 0) {
			PERROR("calloc");
			throw errno;
		}

		cl_float avg_b = avg(in + it->offset - correlation_window.offset,
				correlation_window.len);
		cl_float stdev_b = 0.0;
		for (ulong i = 0; i < correlation_window.len; ++i) {
			register cl_float _b = *(in + it->offset - correlation_window.offset
					+ i) - avg_b;
			stdev_b += POW2(_b);
		}

		// calculate correlation signal locally for each event
		for (std::vector<struct range>::const_iterator it2 = ranges.begin();
				it2 != ranges.end(); ++it2) {
			retrig_correlation_program(it->correlation.signal + it2->set,
					in + it2->set, it2->end - it2->set,
					in + it->offset - correlation_window.offset,
					correlation_window.len, avg_b, stdev_b);
		}
	}
}

/**
 * Form cluster from event
 *
 * @param std::vector<ref_ev_ext_it>::iterator & Reference to trunk event
 * @return std::vector<struct retrig_event> Vector of offsets and correlations of autocorrelating events
 */
retrig_ev Retrigger::form_cluster(
		std::vector<ref_ev_ext_it>::iterator & trunk_it) {
	std::vector<std::pair<size_t, double> > _ev;
	retrig_ev ev;

	for (double _d = 0.95; ev.size() < 3; _d -= 0.025) {
		if (_d < 0.8) {
			errno = ENOSYS;
			fprintf(stderr,
					"[%s:%u] failed to get good enough template events. Refusing to continue.\n",
					__FILE__, __LINE__);
			throw errno;
		}

		ev.clear();
		_ev.clear();

		size_t ref_offset = (*trunk_it)->offset;
		for (std::vector<ref_ev_ext_it>::iterator it =
				(*trunk_it)->cluster.stack.begin();
				it != (*trunk_it)->cluster.stack.end(); ++it) {

			if ((*it)->cluster.assigned) {
				continue;
			}

			off_t offset = 0;
			cl_float dd = 0.0;
			cl_float const *d = (*it)->correlation.signal;
			for (size_t i = ref_offset - lookaround_window.offset, _len = i
					+ lookaround_window.len; i < len && i < _len; ++i) {
				if (d[i] > dd) {
					offset = i;
					dd = d[offset];
				}
			}

			if (dd < 0.6) {
				continue;
			}

			offset -= ref_offset;
#ifdef DEBUG
			printf("template event at %lu, p-value: %f, template offset: %li\n",
					(*it)->offset, dd, offset);
#endif
			(*it)->correlation.offset = offset;

			for (size_t i = 0, _len = len - correlation_window.len; i < _len;
					++i) {
				cl_float const &d = (*it)->correlation.signal[i];
				if (d > _d) {
					cl_float best_d = d;
					size_t offset = i;
					for (size_t j = 1; j < correlation_window.len; ++j) {
						cl_float const &d = (*it)->correlation.signal[i + j];
						if (d > best_d) {
							best_d = d;
							offset = i + j;
						}
					}

					_ev.push_back(
							std::pair<size_t, double>(
									offset - (*it)->correlation.offset,
									(*it)->correlation.p));
				}
			}
		}

		std::sort(_ev.begin(), _ev.end());

		// join nearby events
		static size_t const jitter = 100;
		for (std::vector<std::pair<size_t, double> >::const_iterator it =
				_ev.begin(); it != _ev.end(); ++it) {
			if (it->first == 0) {
				continue;
			}

			size_t _offset = it->first;
			double _p = it->second;
			size_t n = 1;
			for (std::vector<std::pair<size_t, double> >::iterator it2 =
					_ev.begin(); it2 != _ev.end(); ++it2) {
				if (it == it2 || it2->first == 0) {
					continue;
				}

				if (it2->first > it->first - jitter
						&& it2->first < it->first + jitter) {
					_offset += it2->first;
					_p += it2->second;
					++n;
					it2->first = 0;
				}
			}

			struct retrig_event e = { _offset / n, _p / n, *(energy
					+ _offset / n) };
			ev.push_back(e);

#ifdef DEBUG
			printf("retriggered event at %lu, p = %.03f\n", e.offset,
					e.correlation.p);
#endif
		}
	}

	for (std::vector<ref_ev_ext_it>::iterator it =
			(*trunk_it)->cluster.stack.begin();
			it != (*trunk_it)->cluster.stack.end(); ++it) {
		(*it)->cluster.assigned = true;
	}

	return ev;
}

/**
 * A sort function
 *
 * Sorts by size weighted correlation
 *
 * @param std::vector<struct ref_event_ext>::iterator
 * @param std::vector<struct ref_event_ext>::iterator
 * @return bool
 */
bool corr_p_sort(ref_ev_ext_it i, ref_ev_ext_it j) {
	return (i->correlation.p * i->cluster.stack.size()
			> j->correlation.p * j->cluster.stack.size());
}

/**
 * Run autocorrelation for ref_events and form clusters from events
 *
 * Chooses the cluster with the best total correlation as the base event
 * If two strong clusters can be found, they are labeled as S1 ans S2 based on
 * event order.
 * If a distinct separation to S1 and S2 can't be made, the best cluster is labeled as EV.
 * If no cluster has more than three events, all clusters are rejected.
 *
 * @param double const & Cut-off correlation limit
 */
void Retrigger::calc_correlations(double const &correlation_limit) {
	{
		clock_t ref = clock();
		calculate_correlations();
		clock_t t = clock();
		printf("calculated correlations in %.3f s\n",
				static_cast<double>(t - ref) / CLOCKS_PER_SEC);
	}

	std::vector<ref_ev_ext_it> corr_map;

	// sort all correlations by average the highest correlations near the references
	// store correlation offset

	for (ref_ev_ext_it it = evx.begin(); it != evx.end(); ++it) {
		cl_float const *b = it->correlation.signal;
		cl_float dd = 1.0; // self correlation
		size_t _n = 1;

		for (ref_ev_ext_it it2 = evx.begin(); it2 != evx.end(); ++it2) {
			if (it == it2) {
				continue;
			}

			cl_float _dd = 0.0;
			for (size_t i = it2->offset - correlation_window.offset, _len = i
					+ correlation_window.len; i < len && i < _len; ++i) {
				if (b[i] > _dd) {
					_dd = b[i];
				}
			}

			if (_dd < correlation_limit) {
				continue;
			}

			dd += _dd;
			++_n;

			it->cluster.stack.push_back(it2);
		}

		if (_n < 3) {
			continue;
		}
		it->correlation.p = dd / _n;
		it->cluster.assigned = false;

		corr_map.push_back(it);
	}

	if (corr_map.size() == 0) {
		return;
	}

	// select best total correlations
	std::sort(corr_map.begin(), corr_map.end(), corr_p_sort);

//	corr_map.resize(ceil(corr_map.size() / 2));

	std::vector<ref_ev_ext_it>::iterator best_it = corr_map.begin();
	(*best_it)->cluster.assigned = true;

	evs.clear();
	for (std::vector<ref_ev_ext_it>::iterator it = corr_map.begin();
			it != corr_map.end(); ++it) {
		if (it != best_it && (*it)->cluster.assigned) {
			continue;
		}

		printf("template event at %lu, p-value: %f, cluster size: %lu\n",
				(*it)->offset, (*it)->correlation.p,
				(*it)->cluster.stack.size());
		try {
			retrig_ev _ev = form_cluster(it);

			if (evs.size() == 0) {
				evs.push_back(_ev);
				ev = _ev;
				continue;
			}

			// if distance to a better pool is less than 100ms, skip cluster (don't try to join)
			off_t distance = std::numeric_limits<off_t>::max();
			for (std::vector<retrig_ev>::const_iterator it2 = evs.begin();
					it2 != evs.end(); ++it2) {
				for (retrig_ev_it it3 = (*it2).begin(); it3 != (*it2).end();
						++it3) {
					for (retrig_ev_it it4 = _ev.begin(); it4 != _ev.end();
							++it4) {
						off_t d = it4->offset - it3->offset;
						if (labs(d) < labs(distance)) {
							distance = d;
						}
					}
				}
			}

			printf("smallest distance: %.3f s\n",
					static_cast<double>(distance) / sample_freq);

			// calculate rr estimate
			size_t rr = std::numeric_limits<off_t>::max();
			for (retrig_ev_it prev = ev.begin(), it = prev + 1; it != ev.end();
					++it) {
				size_t _rr = it->offset - prev->offset;
				if (_rr < rr) {
					rr = _rr;
				}
				prev = it;
			}

			printf("best cluster rr: %lu (%.03f s)\n", rr, rr / sample_freq);

			if (labs(distance) < static_cast<off_t>(sample_freq * 0.2)) {
				puts("too close to existing, skip");
			} else if (labs(distance) > static_cast<off_t>(sample_freq * 0.5)
					|| labs(distance) > 0.8 * rr) {
				puts("too far to existing, skip");
			} else {
				evs.push_back(_ev);

				if (distance < std::numeric_limits<off_t>::max()) {
					if (distance > 0) {
						s1 = *(evs.begin());
						s2 = _ev;
					} else {
						s1 = _ev;
						s2 = *(evs.begin());
					}
					puts("s1 and s2 models found.");

					// verify that s1 doesn't contain s2 triggers
					if (labs(distance) < sample_freq) {
						retrig_ev_it prev = s1.begin();
						for (retrig_ev::iterator it = s1.begin();
								it != s1.end(); ++it) {
							if (it == prev) {
								continue;
							}

							size_t rr = it->offset - prev->offset;
							prev = it;

							if (rr < 1.25 * distance) {
								printf("bad s1 at %lu, rr = %lu\n", it->offset,
										rr);
								it = s1.erase(it);
								if (it == s1.end()) {
									break;
								}
								prev = it;
							}
						}
					}
					break;
				}
			}
		} catch (int e) {

		}
	}

	for (ref_ev_ext_it it = evx.begin(); it != evx.end(); ++it) {
		mm_free(it->correlation.signal);
	}
	evx.clear();
}

/**
 * Getter for events
 *
 * @return const std::vector<struct retrig_event> & Events in EV cluster
 */
const retrig_ev &Retrigger::get_events() const {
	return ev;
}

/**
 * Getter for events
 *
 * @return const std::vector<struct retrig_event> & Events in S1 cluster
 */
const retrig_ev &Retrigger::get_s1_events() const {
	if (evs.size() < 2) {
		throw EINVAL;
	}
	return s1;
}

/**
 * Getter for events
 *
 * @return const std::vector<struct retrig_event> & Events in S2 cluster
 */
const retrig_ev &Retrigger::get_s2_events() const {
	if (evs.size() < 2) {
		throw EINVAL;
	}
	return s2;
}

/**
 * Getter for energy signal
 *
 * @return const cl_float * Pointer to the energy signal
 */
const cl_float *Retrigger::get_energy() const {
	return energy;
}

/**
 * A sort function
 *
 * Sorts by event offset
 *
 * @param struct ref_event_ext
 * @param struct ref_event_ext
 * @return bool
 */
bool offset_sort(struct ref_event_ext i, struct ref_event_ext j) {
	return (i.offset < j.offset);
}

/**
 * Copies an external reference event vector to an internal event structure.
 *
 * If the number of reference events exceeds /ref_ev_limit/, only the
 * middlemost /ref_ev_events/ are copied.
 *
 * @param std::vector<struct ref_event> const & Vector of reference events
 */
void Retrigger::set_ref_ev(ref_ev const &ev) {
	if (ev.size() < ref_ev_limit) {
		for (ref_ev_it it = ev.begin(); it != ev.end(); ++it) {
			size_t const &max_at = it->offset;
			cl_float const &max = energy[max_at];

			struct ref_event_ext _e = { &(*it), max_at, 0, max, { 0, 0, 0 } };
			evx.push_back(_e);
		}
		return;
	}

	printf("limiting count to %lu\n", ref_ev_limit);
	evx.reserve(ref_ev_limit);

	off_t center = ev.size() / 2;

	{
		ref_event const &e = ev.at(center);
		size_t const &max_at = e.offset;
		cl_float const &max = energy[max_at];

		struct ref_event_ext _e = { &e, max_at, 0, max, { 0, 0, 0 } };
		evx.push_back(_e);
	}

	for (size_t i = 1, _len = ref_ev_limit / 2; i < _len; ++i) {
		{
			ref_event const &e = ev.at(center - i);
			size_t const &max_at = e.offset;
			cl_float const &max = energy[max_at];

			struct ref_event_ext _e = { &e, max_at, 0, max, { 0, 0, 0 } };
			evx.push_back(_e);
		}

		{
			ref_event const &e = ev.at(center + i);
			size_t const &max_at = e.offset;
			cl_float const &max = energy[max_at];

			struct ref_event_ext _e = { &e, max_at, 0, max, { 0, 0, 0 } };
			evx.push_back(_e);
		}
	}
	std::sort(evx.begin(), evx.end(), offset_sort);
}

/**
 * Setter for convolution lookaround window
 *
 * @param double const & Window length in samples
 * @param double const & Window center offset in samples
 */
void Retrigger::set_lookaround_window(double const &len, double const &offset) {
	lookaround_window.len = len * sample_freq;
	lookaround_window.offset = offset * sample_freq;
	lookaround_window.len_in_bytes = lookaround_window.len * sizeof(cl_float);
}

/**
 * Setter for correlation window
 *
 * @param double const & Window length in samples
 * @param double const & Window center offset in samples
 */
void Retrigger::set_correlation_window(double const &len,
		double const &offset) {
	correlation_window.len = len * sample_freq;
	correlation_window.offset = offset * sample_freq;
	correlation_window.len_in_bytes = correlation_window.len * sizeof(cl_float);
}

} /* namespace OpenCL */
