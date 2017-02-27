/**
 * Trigger.h
 *
 *  Created on: Aug 24, 2015
 *      Author: jtmakela
 */

#ifndef SRC_TRIGGER_H_
#define SRC_TRIGGER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "../utils/memory_manager.h"
#include "types.data.h"
#include "types.event.h"

namespace Accbpm {

class Trigger {
public:
	Trigger(cl_float const *energy, size_t const &len, double const &sample_freq);
	virtual ~Trigger();

	ref_ev const &get_events() const;
private:
	size_t skip_offset;
	double limit;
	ref_ev ev;

	double const sample_freq;
	cl_float const *energy;
	size_t len;

	int trig_do();

	double kth_smallest(cl_float a[], ssize_t const n, ssize_t const k);
	int define_threshold_value(cl_float const *energy, int len, int skip_offset,
			double freq, double *threshold);
};

} /* namespace Accbpm */

#endif /* SRC_TRIGGER_H_ */
