/*
 * types.event.h
 *
 *  Created on: Mar 24, 2016
 *      Author: jtmakela
 */

#ifndef INCLUDE_TYPES_EVENT_H_
#define INCLUDE_TYPES_EVENT_H_

#include <vector>

typedef float cl_float;
typedef std::vector<struct ref_event> ref_ev;
typedef std::vector<struct ref_event>::const_iterator ref_ev_it;

struct ref_event {
	size_t offset;
};

struct retrig_event {
	size_t offset;

	struct {
		double p;
	} correlation;

	cl_float nominal_energy;
};

#endif /* INCLUDE_TYPES_EVENT_H_ */
