/**
 * Retrigger.h
 *
 *  Created on: Mar 27, 2016
 *      Author: jtmakela
 */

#ifndef OPENCL_RETRIGGER_H_
#define OPENCL_RETRIGGER_H_

#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>
#include <map>
#include <algorithm>    // std::sort
#include <limits> // std::numeric
#include <time.h> // clock

#include "../utils/memory_manager.h"
#include "../myDSP/Iir.h"
#include "types.data.h"
#include "types.event.h"

namespace Simplified {

struct ref_event_ext {
	struct ref_event const *ref;
	size_t offset;
	off_t change;
	cl_float energy;
	struct {
		cl_float *signal;
		cl_float p;
		size_t offset;
	} correlation;

	struct {
		std::vector<std::vector<struct ref_event_ext>::iterator> stack;
		bool assigned;
	} cluster;
};

struct window {
	size_t len;
	off_t offset;
	size_t len_in_bytes;
};

typedef std::vector<struct ref_event_ext> ref_ev_ext;
typedef ref_ev_ext::iterator ref_ev_ext_it;

typedef std::vector<struct retrig_event> retrig_ev;
typedef retrig_ev::const_iterator retrig_ev_it;

class Retrigger {
public:
	static double const sample_freq = 2000.0;
	static size_t const ref_ev_limit = 100;

	Retrigger(double const &window_length_in_fractions_of_sample_time);
	virtual ~Retrigger();

	virtual void set_data(data_raw_t const *raw, size_t const &len);
	virtual void set_convolution_kernel(cl_float const *kernel,
			size_t const &len);
	virtual void set_ref_ev(ref_ev const &ev);

	virtual void set_lookaround_window(double const &len, double const &offset);
	virtual void set_correlation_window(double const &len,
			double const &offset);

	virtual void calc(struct ref_event const &e, off_t const &offset);
	void calc_correlations();
	void calc_correlations(double const &correlation_limit);

	ref_ev_ext const &get_extended_ref_events() const;
	retrig_ev const &get_events() const;
	retrig_ev const &get_s1_events() const;
	retrig_ev const &get_s2_events() const;

	cl_float const *get_energy() const;

protected:
	data_raw_t const *raw;
	size_t len;
	size_t len_in_bytes;

	cl_float *in;
	cl_float *energy;

	void calculate_energy();
	void calculate_correlations();

private:
	struct window energy_window;
	struct window convolution_window;
	struct window lookaround_window;
	struct window correlation_window;

	cl_float *blackman_window;
	cl_float *retrig_convolution_kernel;

	ref_ev_ext evx; // temporary resource, lifespan from calc() to calc_correlations()
	off_t evx_offset; // average offset to ref_event offset

	retrig_ev ev;
	std::vector<retrig_ev> evs;

	retrig_ev s1;
	retrig_ev s2;

	retrig_ev form_cluster(std::vector<ref_ev_ext_it>::iterator &trunk_it);
};

} /* namespace OpenCL */

#endif /* OPENCL_RETRIGGER_H_ */
