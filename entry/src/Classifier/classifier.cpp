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
 * classifier.cpp
 *
 *  Created on: Apr 4, 2016
 *      Author: hvaanane
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <float.h>

#include "markers.h"
#include "../utils/memory_manager.h"
#include "classifier.h"
//#define VERBOSE 1

namespace Classifier {

__attribute_used__ int free_string_tree(struct string_tree *tree) {
	mm_free(tree->nodes);
	tree->nodes = NULL;
	tree->n_nodes = 0;
	return (0);
}

__attribute_used__ int save_string_tree(char const *filename,
		struct string_tree const *tree) {
	FILE *fp = fopen(filename, "w");
	fwrite(tree, sizeof(struct string_tree), 1, fp);
	fwrite(tree->nodes, sizeof(struct string_node), tree->n_nodes, fp);
	fclose(fp);
	return (0);
}

__attribute_used__ int load_string_tree(char const *filename,
		struct string_tree *tree) {
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("Oops, cannot open the treefile %s\n", filename);
		return (-1);
	}
	fread(tree, sizeof(struct string_tree), 1, fp);
	tree->nodes = static_cast<struct string_node *>(mm_malloc(
			tree->n_nodes * sizeof(struct string_node)));
	if (fread(tree->nodes, sizeof(struct string_node), tree->n_nodes, fp)
			!= tree->n_nodes) {
		printf("Oops, treefile %s is corrupted\n", filename);
		return (-1);
	}
	fclose(fp);
	return (0);
}

__attribute_used__ int load_txt_string_tree(char const *filename,
		struct string_tree *tree) {
	struct stat sb;
	if (stat(filename, &sb)) {
		printf("Oops, cannot open the treefile %s\n", filename);
		return (-1);
	}

	char cachefile[FILENAME_MAX];
	sprintf(cachefile, "/tmp/Classifier.tree.%lu.dat", sb.st_ino);
	if (stat(cachefile, &sb) == 0) {
		return load_string_tree(cachefile, tree);
	}

	FILE *fp = fopen(filename, "r");
	size_t i;
	if (fp == NULL) {
		printf("Oops, cannot open the treefile %s\n", filename);
		return (-1);
	}
	fscanf(fp, "%d\t%d\n", &(tree->n_nodes), &(tree->n_classes));
	tree->nodes = static_cast<struct string_node *>(mm_malloc(
			tree->n_nodes * sizeof(struct string_node)));
	for (i = 0; i < tree->n_nodes; i++) {
		if (fscanf(fp, "%s\t%lf\t%d\t%d\t%d\n", tree->nodes[i].marker_name,
				&(tree->nodes[i].split_value), &(tree->nodes[i].up),
				&(tree->nodes[i].left), &(tree->nodes[i].right)) < 5) {
			printf("Oops, treefile %s is corrupted in line %lu\n", filename,
					i + 1);
			return (-1);
		}
	}
	fclose(fp);

	save_string_tree(cachefile, tree);

	return (0);
}

__attribute_used__ int do_with_strings(struct data *data,
		Simplified::retrig_ev const *ev1, Simplified::retrig_ev const *ev2,
		struct string_tree *tree, int this_node) {
	double marker_value;
	int next_node;
	if (Markers::create_named(tree->nodes[this_node].marker_name, data, ev1,
			ev2, &marker_value) < 0) { // error
		return (-tree->n_classes); // return CLASSIFIER_UNKNOWN
	}
	if (tree->nodes[this_node].split_value >= marker_value) { // go left?
		next_node = tree->nodes[this_node].left;
	} else { // go right
		next_node = tree->nodes[this_node].right;
	}
	if (next_node > 0) { // continue to next node?
		return (do_with_strings(data, ev1, ev2, tree, next_node));
	} else { // we found the leaf
		return (next_node); // the actual return value is coded (multipled with -1) to the negative node values
	}
}

/*
 * A dynamically created tree classifier. See the article for further details.
 *
 * @param struct data * Pointer to ECG legacy data structure
 * @param std::vector<struct retrig_event> const * Pointer to mandatory first event cluster
 * @param std::vector<struct retrig_event> const * Pointer to optional second event cluster
 * @param char * Marker tree filename to be used
 * @return enum result_e Classifier result; normal, abnormal or unknown
 */
enum result_e classify_this(struct data *data, Simplified::retrig_ev const *ev1,
		Simplified::retrig_ev const *ev2, char *tree_name) {
	struct string_tree tree = { };

	load_txt_string_tree(tree_name, &tree);
	enum result_e result = static_cast<enum result_e>(do_with_strings(data, ev1,
			ev2, &tree, 0));
	free_string_tree(&tree);
	return (result);
}

} // namespace Classifier
