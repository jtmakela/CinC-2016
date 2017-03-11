/*
 * types.data.h
 *
 *  Created on: Apr 8, 2015
 *      Author: jtmakela
 */

#ifndef INCLUDE_TYPES_DATA_H_
#define INCLUDE_TYPES_DATA_H_

#include <cstdlib>

#ifndef ssize_t
typedef int64_t ssize_t;
#endif

#ifndef cl_float
typedef float cl_float;
#endif

typedef double data_raw_t;
//typedef int16_t data_raw_t;

struct data_channel {
	data_raw_t *raw;
	double samples_per_mv;
};

#ifndef STRUCT_DATA_IS_DEFINED
struct data {
#define STRUCT_DATA_IS_DEFINED
#else
	struct decoder_data {
#endif
	size_t samples_per_channel;
	size_t number_of_channels;
	double sample_freq;
	time_t time_offset;

	struct data_channel *ch;
};

#endif /* INCLUDE_TYPES_DATA_H_ */
