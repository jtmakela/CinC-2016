/**
 * classifier.h
 *
 *  Created on: Apr 28, 2016
 *      Author: hvaanane
 */

#ifndef CLASSIFIER_H_
#define CLASSIFIER_H_

#include "../Simplified/Retrigger.h"

namespace Classifier {

#define MAX_MARKERNAME_LEN 60
struct string_node {
	double split_value;
	char marker_name[MAX_MARKERNAME_LEN];
	int up;
	int left;
	int right;
};

struct string_tree {
	struct string_node *nodes;
	int n_classes;
	int n_nodes;
};

enum result_e {
	normal = -3, abnormal = -2, unknown = -1,
};

enum result_e classify_this(struct data *data, Simplified::retrig_ev const *ev1,
		Simplified::retrig_ev const *ev2, char *tree_name);

} // namespace Classifier

#endif /* CLASSIFIER_H_ */
