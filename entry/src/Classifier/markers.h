/**
 * markers.h
 *
 *  Created on: Apr 15, 2016
 *      Author: hvaanane
 */

#ifndef MARKERS_H_
#define MARKERS_H_

#include "../Simplified/Retrigger.h"

#define MARKERS_S1S2 1
#define MARKERS_EVENTS 2
#define MARKERS_UNTRIGGED 3

namespace Classifier {
namespace Markers {
int create_named(char const *name, struct data *data,
		Simplified::retrig_ev const *ev1, Simplified::retrig_ev const *ev2,
		double *marker_value);
} // namespace Markers
} // namespace Classifier

#endif /* MARKERS_H_ */
