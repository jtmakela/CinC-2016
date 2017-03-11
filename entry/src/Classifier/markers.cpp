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
 * markers.cpp
 *
 *  Created on: Apr 6, 2016
 *      Author: hvaanane
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <unistd.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#include "../utils/memory_manager.h"
#include "../myDSP/Iir.h"
#include "markers.h"
#include "classifier.h"

#if 1
#ifndef POW2
#define POW2(s) ((s)*(s))
#endif

namespace Classifier {
namespace Markers {

static double constexpr conf_s_start = -0.100;
static double constexpr conf_s_end = 0.100;
static double constexpr conf_margin = 0.050;
static double constexpr conf_moving_std_len = 0.100;
static double constexpr conf_default_rr = 0.800;
static double constexpr conf_default_s1s2_dur = 0.400;
static double constexpr conf_win_len = 3.0;
static double constexpr conf_ignore_from_start = 1.0;

struct double_array {
	double *data;
	size_t len;
};

struct markers_s1s2 {
	double s1;
	double s2;
	double s1s2;
	double s2s1;
	double s2s1_min;
	double s2s1_max;
	double s2s1_minmax;
};

struct markers_s {
	double s;
	double ss;
	double ss_min;
	double ss_max;
	double ss_minmax;
};

static inline void set_double_array_len(struct double_array *data, int len,
		int uninitialized) {
	if (uninitialized) {
		data->data =
				static_cast<double *>(mm_malloc(data->len * sizeof(double)));
	} else if (data->len < len) {
		data->data = static_cast<double *>(mm_realloc(data->data,
				len * sizeof(double)));
	}
	data->len = len;
}

static double get_kth_biggest(double a[], size_t n, size_t k) {
	long int i, j, l, m;
	double x, tmp;
	l = 0;
	m = n - 1;
	while (l < m) {
		x = a[k];
		i = l;
		j = m;
		do {
			while (a[i] < x)
				i++;
			while (x < a[j])
				j--;
			if (i <= j) {
				tmp = a[i]; // swap
				a[i++] = a[j];
				a[j--] = tmp;
			}
		} while (i <= j);
		if (j < k)
			l = i;
		if (k < i)
			m = j;
	}
	return a[k];
}

static double get_std(struct double_array *data, size_t start, int std_len) {
	double power = 0.0;
	double ave = 0.0;
	size_t i;
	if (start + std_len > data->len) {
		std_len = data->len - start;
		if (std_len < 3) {
			printf("Oops, too little data for std\n");
			return (-1.0);
		}
	}
	for (i = start; i < start + std_len; i++) {
		power += POW2(data->data[i]);
		ave += data->data[i];
	}
	ave /= (double) (std_len);
	power = power / (double) (std_len) - POW2(ave);
	if (power > 0.0) { // this check is needed due to rounding errors when using floating point arithmetic..
		return (sqrt(power));
	}
	return (0.0);
}

static int define_moving_std(struct double_array *data,
		struct double_array *std, int std_len) {
	size_t i;
	double sum, sum2;
	double ave;
	double inv_win_len = 1.0 / (double) (std_len);
	int half_win = std_len / 2;

	set_double_array_len(std, data->len, 0);
	memset(std->data, 0, std->len * sizeof(double));

	sum = sum2 = 0.0;
	for (i = 0; i < std_len; i++) {
		sum += data->data[i];
		sum2 += POW2(data->data[i]);
		ave = sum / (double) (i + 1);
		if (i >= half_win) {
			std->data[i - half_win] = sum2 / (double) (i + 1) - POW2(ave);
			if (std->data[i - half_win] > 0.0) {
				std->data[i - half_win] = sqrt(std->data[i - half_win]);
			}
		}
	}
	for (; i < data->len; i++) {
		sum += data->data[i] - data->data[i - std_len];
		sum2 += POW2(data->data[i]) - POW2(data->data[i - std_len]);
		ave = sum * inv_win_len;
		std->data[i - half_win] = sum2 * inv_win_len - POW2(ave);
		if (std->data[i - half_win] > 0.0) {
			std->data[i - half_win] = sqrt(std->data[i - half_win]);
		}
	}
	for (; i < data->len + half_win; i++) {
		sum -= data->data[i - std_len];
		sum2 -= POW2(data->data[i - std_len]);
		ave = sum / (double) (std_len - (i - data->len));
		std->data[i - half_win] = sum2
				/ (double) (std_len - (i - data->len))-POW2(ave);
		if (std->data[i - half_win] > 0.0) {
			std->data[i - half_win] = sqrt(std->data[i - half_win]);
		}
	}
	return (0);
}

static int get_events_stds(struct double_array *data,
		Simplified::retrig_ev const *events, int win_start, int win_end,
		struct double_array *std) {
	ssize_t j, n, start, end;

	set_double_array_len(std, events->size(), 0);

	for (n = 0, j = 0; j < events->size(); j++) {
		start = (*events)[j].offset + win_start;
		if (start < 0) {
			continue;
		}
		end = (*events)[j].offset + win_end;
		if (end > data->len) {
			continue;
		}
		std->data[n++] = get_std(data, start, end - start);
	}
	std->len = n;
	if (n == 0) {
		printf(
				"ERROR_std: there is data with no events far enough from data bounders\n Please remove it from the list!");
	}
	return (0);
}

static int get_events_absmax(struct double_array *data,
		Simplified::retrig_ev const *events, int win_start, int win_end,
		struct double_array *absmax) {
	ssize_t i, j, n, start, end;

	set_double_array_len(absmax, events->size(), 0);

	for (n = 0, j = 0; j < events->size(); j++) {
		start = (*events)[j].offset + win_start;
		if (start < 0) {
			continue;
		}
		end = (*events)[j].offset + win_end;
		if (end > data->len) {
			continue;
		}
		absmax->data[n] = fabs(data->data[start]);
		for (i = start + 1; i < end; i++) {
			if (absmax->data[n] < fabs(data->data[i])) {
				absmax->data[n] = fabs(data->data[i]);
			}
		}
		n++;
	}
	absmax->len = n;
	if (n == 0) {
		printf(
				"ERROR_extreme: there is data with no events far enough from data bounders\n Please remove it from the list!");
	}
	return (0);
}

static int get_events_width(struct double_array *data,
		Simplified::retrig_ev const *events, int win_start, int win_end,
		double limit, struct double_array *width, double freq) {
	ssize_t j, n, start, end;
	size_t event_start, event_end;

	set_double_array_len(width, events->size(), 0);

	for (n = 0, j = 0; j < events->size(); j++) {
		start = (*events)[j].offset + win_start;
		if (start < 0) {
			continue;
		}
		end = (*events)[j].offset + win_end;
		if (end > data->len) {
			continue;
		}
		event_start = event_end = 0;
		for (event_start = start;
				event_start < end && data->data[event_start] < limit;
				event_start++)
			;
		for (event_end = end - 1;
				event_end > start && data->data[event_end] < limit; event_end--)
			;
		width->data[n] = ((double) (event_end) - (double) (event_start)) / freq;
		n++;
	}
	width->len = n;
	if (n == 0) {
		printf(
				"ERROR_extreme: there is data with no events far enough from data bounders\n Please remove it from the list!");
	}
	return (0);
}

static int get_events_extreme_values(struct double_array *data,
		Simplified::retrig_ev const *events, int win_start, int win_end,
		struct double_array *min, struct double_array *max,
		struct double_array *min_max) {
	ssize_t i, j, n, start, end;

	set_double_array_len(min, events->size(), 0);
	set_double_array_len(max, events->size(), 0);
	set_double_array_len(min_max, events->size(), 0);

	for (n = 0, j = 0; j < events->size(); j++) {
		start = (*events)[j].offset + win_start;
		if (start < 0) {
			continue;
		}
		end = (*events)[j].offset + win_end;
		if (end > data->len) {
			continue;
		}
		min->data[n] = max->data[n] = data->data[start];
		for (i = start + 1; i < end; i++) {
			if (min->data[n] > data->data[i]) {
				min->data[n] = data->data[i];
			} else if (max->data[n] < data->data[i]) {
				max->data[n] = data->data[i];
			}
		}
		min_max->data[n] = max->data[n] - min->data[n];
		n++;
	}
	min->len = n;
	max->len = n;
	min_max->len = n;
	if (n == 0) {
		printf(
				"ERROR_extreme: there is data with no events far enough from data bounders\n Please remove it from the list!");
	}
	return (0);
}

inline static int is_too_saturated(size_t n_saturated, size_t dur_saturated,
		size_t tot_n_samp) {
	//	return(sqrt(n_saturated)*dur_saturated>0.01*tot_n_samp);
	return (sqrt(10000.0 * n_saturated / tot_n_samp) * dur_saturated
			> 0.01 * tot_n_samp);
}

int is_data_saturated(struct data *data, size_t *tot_n_saturated,
		size_t *tot_len_saturated) {
	int ch;
	size_t i;
	double range_min, range_max;
	size_t n_min_saturated, n_max_saturated;
	int min_n_saturated = 2;
	int skip_from_start = 500.0;
	(*tot_n_saturated) = 0;
	(*tot_len_saturated) = 0;

	// specify limits from data
	range_min = range_max = data->ch[0].raw[0];
	for (ch = 0; ch < data->number_of_channels; ch++) {
		for (i = 0; i < data->samples_per_channel; i++) {
			if (range_max < data->ch[ch].raw[i]) {
				range_max = data->ch[ch].raw[i];
			} else if (range_min > data->ch[ch].raw[i]) {
				range_min = data->ch[ch].raw[i];
			}
		}
	}
	// Calculate the number and total duration of the saturated segments
	for (ch = 0; ch < data->number_of_channels; ch++) {
		for (i = skip_from_start; i < data->samples_per_channel; i++) {
			if (range_max
					== data->ch[ch].raw[i]/* && n_max_saturated<min_n_saturated*/) {
				for (n_max_saturated = 1;
						range_max <= data->ch[ch].raw[i + n_max_saturated];
						n_max_saturated++, i++)
					;
				if (n_max_saturated > min_n_saturated) {
					(*tot_n_saturated)++;
					(*tot_len_saturated) += n_max_saturated;
				}
			} else if (range_min
					== data->ch[ch].raw[i]/* && n_min_saturated<min_n_saturated*/) {
				for (n_min_saturated = 1;
						range_min >= data->ch[ch].raw[i + n_min_saturated];
						n_min_saturated++, i++)
					;
				if (n_min_saturated > min_n_saturated) {
					(*tot_n_saturated)++;
					(*tot_len_saturated) += n_min_saturated;
				}
			}
		}
	}
	return (is_too_saturated((*tot_n_saturated), (*tot_len_saturated),
			data->samples_per_channel));
}

static int get_repeating_extreme_values(struct double_array *data,
		int ignore_from_start, int win_len, struct double_array *min,
		struct double_array *max, struct double_array *min_max) {
	size_t n_win = (data->len - ignore_from_start) / win_len;
	size_t i, n, start;

	set_double_array_len(min, n_win, 0);
	set_double_array_len(max, n_win, 0);
	set_double_array_len(min_max, n_win, 0);

	for (n = 0; n < n_win; n++) {
		start = ignore_from_start + n * win_len;
		min->data[n] = max->data[n] = data->data[start];
		for (i = start + 1; i < start + win_len; i++) {
			if (min->data[n] > data->data[i]) {
				min->data[n] = data->data[i];
			} else if (max->data[n] < data->data[i]) {
				max->data[n] = data->data[i];
			}
		}
		min_max->data[n] = max->data[n] - min->data[n];
	}
	return (0);
}

void free_data(struct data *data) {
	for (int ch = 0; ch < data->number_of_channels; ch++) {
		mm_free(data->ch[ch].raw);
	}
	mm_free(data->ch);
	data->ch = NULL;
	data->number_of_channels = 0;

}

static inline double define_s_rr(Simplified::retrig_ev const *events,
		double freq) {
	size_t min_rr = 0.6 * freq;
	size_t max_rr = 2.2 * freq;
	size_t i, n = 0;
	size_t this_rr;
	double *rr = static_cast<double *>(mm_malloc(
			events->size() * sizeof(double)));
	double median_rr = -1.0;
	for (i = 1; i < events->size(); i++) {
		this_rr = (*events)[i].offset - (*events)[i - 1].offset;
		if (this_rr > min_rr && this_rr < max_rr) {
			rr[n++] = this_rr;
		}
	}
	if (n > 0) {
		if (n == 1) {
			median_rr = rr[0];
		} else {
			median_rr = get_kth_biggest(rr, n, n / 2);
		}
	}
	mm_free(rr);
	return (median_rr / freq);
}

static inline double get_ss(Simplified::retrig_ev const *events, double freq) {
	double dur = define_s_rr(events, freq);
	if (dur < 0.0) {
		dur = conf_default_rr;
	}
	return (dur);
}

static inline double define_s1s2_dur(Simplified::retrig_ev const *s1_events,
		Simplified::retrig_ev const *s2_events, double freq) {
	size_t min_dur = 0.200 * freq;
	size_t max_dur = 0.600 * freq;
	size_t i, j, n = 0;
	size_t this_dur;
	double *dur = static_cast<double *>(mm_malloc(
			s1_events->size() * sizeof(double)));
	double median_dur = -1.0;
	for (i = 0, j = 0; i < s1_events->size(); i++) {
		for (;
				j < s2_events->size()
						&& (*s2_events)[j].offset < (*s1_events)[i].offset; j++)
			;
		this_dur = (*s2_events)[j].offset - (*s1_events)[i].offset;
		if (this_dur > min_dur && this_dur < max_dur) {
			dur[n++] = this_dur;
		}
	}
	if (n > 0) {
		if (n == 1) {
			median_dur = dur[0];
		} else {
			median_dur = get_kth_biggest(dur, n, n / 2);
		}
	}
	mm_free(dur);
	return (median_dur / freq);
}

static inline double get_s1s2_dur(Simplified::retrig_ev const *s1_events,
		Simplified::retrig_ev const *s2_events, double freq) {
	double dur = define_s1s2_dur(s1_events, s2_events, freq);
	if (dur < 0.0) {
		dur = conf_default_s1s2_dur;
	}
	return (dur);
}

#define MARKERS_N_WIDTH_LEVELS 5

static int get_named_value(char const *where, char const *how, struct double_array *ddata,
		double sfreq, Simplified::retrig_ev const *s1_events,
		Simplified::retrig_ev const *s2_events, double *marker_value) {
	static struct double_array tmp_array = { }; // set as static, so it doesn't have to be realloced again every time (TODO free).
	static struct double_array std_array = { }; // set as static, so it doesn't have to be realloced again every time (TODO free).
	static struct double_array min = { }; // set as static, so it doesn't have to be realloced again every time (TODO free).
	static struct double_array max = { }; // set as static, so it doesn't have to be realloced again every time (TODO free).
	static struct double_array minmax = { }; // set as static, so it doesn't have to be realloced again every time (TODO free).
	double s1s2_dur, ss_dur;

	if (!strcmp(how, "all")) {
		if (!strcmp(where, "s1")) {
			get_events_stds(ddata, s1_events, conf_s_start * sfreq,
					conf_s_end * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
		} else if (!strcmp(where, "s2")) {
			get_events_stds(ddata, s2_events, conf_s_start * sfreq,
					conf_s_end * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
		} else if (!strcmp(where, "s")) {
			get_events_stds(ddata, s1_events, conf_s_start * sfreq,
					conf_s_end * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
			if (s2_events != NULL) { // define ave of s1 and s2 values
				get_events_stds(ddata, s2_events, conf_s_start * sfreq,
						conf_s_end * sfreq, &tmp_array);
				*marker_value = (get_kth_biggest(tmp_array.data, tmp_array.len,
						tmp_array.len / 2) + *marker_value) / 2.0;
			}
		} else if (!strcmp(where, "s1s2")) {
			s1s2_dur = get_s1s2_dur(s1_events, s2_events, sfreq);
			get_events_stds(ddata, s1_events,
					(conf_s_end + conf_margin) * sfreq,
					(s1s2_dur - conf_margin) * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
		} else if (!strcmp(where, "s2s1")) {
			s1s2_dur = get_s1s2_dur(s1_events, s2_events, sfreq);
			ss_dur = get_ss(s1_events, sfreq);
			get_events_stds(ddata, s2_events,
					(-conf_margin - ss_dur + s1s2_dur) * sfreq,
					(-conf_margin) * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
		} else if (!strcmp(where, "ss")) {
			ss_dur = get_ss(s1_events, sfreq);
			get_events_stds(ddata, s1_events,
					(conf_s_end + conf_margin) * sfreq,
					(ss_dur - conf_margin) * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
		} else if (!strcmp(where, "base")) {
			get_events_stds(ddata, s1_events, (-0.125) * sfreq,
					(-0.075) * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
		} else if (!strcmp(where, "q1")) {
			get_events_stds(ddata, s1_events, -0.125 * sfreq, -0.075 * sfreq,
					&tmp_array);
			double base_val = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
			s1s2_dur = get_s1s2_dur(s1_events, s2_events, sfreq);
			get_events_stds(ddata, s1_events, (s1s2_dur / 4.0 - 0.025) * sfreq,
					(s1s2_dur / 4.0 + 0.025) * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2) - base_val;
		} else if (!strcmp(where, "q2")) {
			get_events_stds(ddata, s1_events, -0.125 * sfreq, -0.075 * sfreq,
					&tmp_array);
			double base_val = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
			s1s2_dur = get_s1s2_dur(s1_events, s2_events, sfreq);
			get_events_stds(ddata, s1_events,
					(2.0 * s1s2_dur / 4.0 - 0.025) * sfreq,
					(2.0 * s1s2_dur / 4.0 + 0.025) * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2) - base_val;
		} else if (!strcmp(where, "q3")) {
			get_events_stds(ddata, s1_events, -0.125 * sfreq, -0.075 * sfreq,
					&tmp_array);
			double base_val = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
			s1s2_dur = get_s1s2_dur(s1_events, s2_events, sfreq);
			get_events_stds(ddata, s1_events,
					(3.0 * s1s2_dur / 4.0 - 0.025) * sfreq,
					(3.0 * s1s2_dur / 4.0 + 0.025) * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2) - base_val;
		} else if (!strcmp(where, "q5")) {
			get_events_stds(ddata, s1_events, -0.125 * sfreq, -0.075 * sfreq,
					&tmp_array);
			double base_val = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
			s1s2_dur = get_s1s2_dur(s1_events, s2_events, sfreq);
			get_events_stds(ddata, s2_events, (s1s2_dur / 4.0 - 0.025) * sfreq,
					(s1s2_dur / 4.0 + 0.025) * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2) - base_val;
		} else if (!strcmp(where, "q6")) {
			get_events_stds(ddata, s1_events, -0.125 * sfreq, -0.075 * sfreq,
					&tmp_array);
			double base_val = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
			s1s2_dur = get_s1s2_dur(s1_events, s2_events, sfreq);
			get_events_stds(ddata, s2_events,
					(2.0 * s1s2_dur / 4.0 - 0.025) * sfreq,
					(2.0 * s1s2_dur / 4.0 + 0.025) * sfreq, &tmp_array);
			*marker_value = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2) - base_val;
		}

		else if (!strcmp(where, "-")) {

		} else {
			printf("Oops! Unsupported where (%s) for %s\n", where, how);
			return (-1);
		}
		return (0);
	} else { // how == min, max or minmax
		define_moving_std(ddata, &std_array, conf_moving_std_len * sfreq);

		if (!strcmp(where, "s1s2")) {
			s1s2_dur = get_s1s2_dur(s1_events, s2_events, sfreq);
			get_events_extreme_values(&std_array, s1_events,
					(conf_s_end + conf_margin) * sfreq,
					(s1s2_dur - conf_margin) * sfreq, &min, &max, &minmax);
			if (!strcmp(how, "min")) {
				*marker_value = get_kth_biggest(min.data, min.len, min.len / 2); // take the median value
			} else if (!strcmp(how, "max")) {
				*marker_value = get_kth_biggest(max.data, max.len, max.len / 2); // take the median value
			} else if (!strcmp(how, "minmax")) {
				*marker_value = get_kth_biggest(minmax.data, minmax.len,
						minmax.len / 2); // take the median value
			} else {
				printf("Oops! Unsupported how (%s) for %s\n", how, where);
			}
		} else if (!strcmp(where, "s2s1")) {
			s1s2_dur = get_s1s2_dur(s1_events, s2_events, sfreq);
			ss_dur = get_ss(s1_events, sfreq);
			get_events_extreme_values(&std_array, s1_events,
					(-conf_margin - ss_dur + s1s2_dur) * sfreq,
					(-conf_margin) * sfreq, &min, &max, &minmax);
			if (!strcmp(how, "min")) {
				*marker_value = get_kth_biggest(min.data, min.len, min.len / 2); // take the median value
			} else if (!strcmp(how, "max")) {
				*marker_value = get_kth_biggest(max.data, max.len, max.len / 2); // take the median value
			} else if (!strcmp(how, "minmax")) {
				*marker_value = get_kth_biggest(minmax.data, minmax.len,
						minmax.len / 2); // take the median value
			} else {
				printf("***Oops! Unsupported how (%s) for %s\n", how, where);
			}
		} else if (!strcmp(where, "ss")) {
			ss_dur = get_ss(s1_events, sfreq);
			get_events_extreme_values(&std_array, s1_events,
					(conf_s_end + conf_margin) * sfreq,
					(ss_dur - conf_margin) * sfreq, &min, &max, &minmax);
			if (!strcmp(how, "min")) {
				*marker_value = get_kth_biggest(min.data, min.len, min.len / 2); // take the median value
			} else if (!strcmp(how, "max")) {
				*marker_value = get_kth_biggest(max.data, max.len, max.len / 2); // take the median value
			} else if (!strcmp(how, "minmax")) {
				*marker_value = get_kth_biggest(minmax.data, minmax.len,
						minmax.len / 2); // take the median value
			} else {
				printf("Oops! Unsupported how (%s) for %s\n", how, where);
			}
		} else if (!strcmp(where, "untrigged")) {
			get_repeating_extreme_values(&std_array,
					conf_ignore_from_start * sfreq, conf_win_len * sfreq, &min,
					&max, &minmax);
			if (!strcmp(how, "min")) {
				*marker_value = get_kth_biggest(min.data, min.len, min.len / 2); // take the median value
			} else if (!strcmp(how, "max")) {
				*marker_value = get_kth_biggest(max.data, max.len, max.len / 2); // take the median value
			} else if (!strcmp(how, "minmax")) {
				*marker_value = get_kth_biggest(minmax.data, minmax.len,
						minmax.len / 2); // take the median value
			} else {
				printf("Oops! Unsupported how (%s) for %s\n", how, where);
			}
		}

		else {
			printf("Oops! Unsupported where (%s) for %s\n", where, how);
			return (-1);
		}
		return (0);
	}
	return (0);
}

int create_named(char const *name, struct data *data,
		Simplified::retrig_ev const *s1_events,
		Simplified::retrig_ev const *s2_events, double *marker_value) {
	char what[10], where[10], to[10], how[10];
	double freq1, freq2;
	double help_val;
	struct double_array filtered = { 0 };
	static struct double_array tmp_array = { 0 }; // set as static, so it doesn't have to be realloced again every time (TODO free).

	if (sscanf(name, "%[^_]_%[^_]_%[^_]_%[^_]_%lf_%lf", what, where, to, how,
			&freq1, &freq2) < 6) {
		printf("Unknown marker %s. Not enough name parts.\n", name);
		return (-1);
	}

	filtered.len = 0;
	filtered.data = NULL;
	if (freq2 <= 0.0) { // no filtering
		set_double_array_len(&filtered, data->samples_per_channel, 0);
		memcpy(filtered.data, data->ch[0].raw,
				data->samples_per_channel * sizeof(double));
	} else {
		DSP::Iir iir;
		iir.bandpass(&(filtered.data), data->ch[0].raw,
				data->samples_per_channel, freq1, freq2, 0.5, 4,
				data->sample_freq);
		filtered.len = data->samples_per_channel;
	}

	if (!strcmp(what, "abs")) {
		get_named_value(where, how, &filtered, data->sample_freq, s1_events,
				s2_events, marker_value);
	} else if (!strcmp(what, "rel")) {
		get_named_value(where, how, &filtered, data->sample_freq, s1_events,
				s2_events, marker_value);
		get_named_value(to, "all", &filtered, data->sample_freq, s1_events,
				s2_events, &help_val);
		if (help_val != 0.0) {
			(*marker_value) /= help_val;
		} else {
			printf("Oops! Trying to divide by zero when defining %s\n", name);
			(*marker_value) *= 10000000000.0;
		}
	} else if (!strcmp(what, "corr")) {
		get_named_value(where, how, &filtered, data->sample_freq, s1_events,
				s2_events, marker_value);
		get_named_value("base", "all", &filtered, data->sample_freq, s1_events,
				s2_events, &help_val);
		(*marker_value) -= help_val;
	} else if (!strcmp(what, "relcorr")) {
		get_named_value(where, how, &filtered, data->sample_freq, s1_events,
				s2_events, marker_value);
		get_named_value("base", "all", &filtered, data->sample_freq, s1_events,
				s2_events, &help_val);
		(*marker_value) -= help_val;
		get_named_value(to, "all", &filtered, data->sample_freq, s1_events,
				s2_events, &help_val);
		if (help_val != 0.0) {
			(*marker_value) /= help_val;
		} else {
			printf("Oops! Trying to divide by zero when defining %s\n", name);
			(*marker_value) *= 10000000000.0;
		}
	} else if (!strcmp(what, "norm")) {
		get_named_value(where, how, &filtered, data->sample_freq, s1_events,
				s2_events, marker_value);
		memcpy(filtered.data, data->ch[0].raw,
				data->samples_per_channel * sizeof(double)); // create not filtered for creating all band normalization reference to bandpassed marker
		get_named_value(where, how, &filtered, data->sample_freq, s1_events,
				s2_events, &help_val);
		if (help_val != 0.0) {
			(*marker_value) /= help_val;
		} else {
			printf("Oops! Trying to divide by zero when defining %s\n", name);
			(*marker_value) *= 10000000000.0;
		}
	} else if (!strcmp(what, "dur")) {
		if (!strcmp(where, "s1s2")) {
			(*marker_value) = get_s1s2_dur(s1_events, s2_events,
					data->sample_freq);
		} else if (!strcmp(where, "ss")) {
			(*marker_value) = get_ss(s1_events, data->sample_freq);
		} else {
			printf("Unknown marker %s. Dur us defined only for s1s2 or ss.\n",
					name);
		}
	} else if (!strcmp(what, "width")) {
		double perc_lvl;
		if (!strcmp(where, "s1") || !strcmp(where, "s")) {
			get_events_absmax(&filtered, s1_events,
					conf_s_start * data->sample_freq,
					conf_s_end * data->sample_freq, &tmp_array);
			help_val = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2); // define absmax value
			if (sscanf(to, "%lf", &perc_lvl) < 1) {
				printf("Unknown marker %s. Can't read to (width lvl value).\n",
						name);
				return (-1);
			}
			get_events_width(&filtered, s1_events,
					(conf_s_start - conf_margin) * data->sample_freq,
					(conf_s_end + conf_margin) * data->sample_freq,
					perc_lvl / 100.0 * help_val, &tmp_array, data->sample_freq);
			(*marker_value) = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
		} else if (!strcmp(where, "s2")) {
			get_events_absmax(&filtered, s2_events,
					conf_s_start * data->sample_freq,
					conf_s_end * data->sample_freq, &tmp_array);
			help_val = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2); // define absmax value
			if (sscanf(to, "%lf", &perc_lvl) < 1) {
				printf("Unknown marker %s. Can't read to (width lvl value).\n",
						name);
				return (-1);
			}
			get_events_width(&filtered, s2_events,
					(conf_s_start - conf_margin) * data->sample_freq,
					(conf_s_end + conf_margin) * data->sample_freq,
					perc_lvl / 100.0 * help_val, &tmp_array, data->sample_freq);
			(*marker_value) = get_kth_biggest(tmp_array.data, tmp_array.len,
					tmp_array.len / 2);
		} else {
			printf(
					"Unknown marker %s. Width is defined only for s1,s2 and s.\n",
					name);
			return (-1);
		}
	} else if (!strcmp(what, "ext")) {
		printf("Oops! External markers not defined quite yet..\n"); // TODO: ...
	} else {
		printf("Unknown marker %s. Can't understand what.\n", name);
		return (-1);
	}
	mm_free(filtered.data);
	return (0);
}

} // namespace Markers
} // namespace Classifier

#endif

