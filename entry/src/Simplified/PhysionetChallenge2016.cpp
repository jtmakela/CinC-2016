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
 * PhysionetChallenge2016.cpp
 *
 *  Created on: Apr 28, 2016
 *      Author: jtmakela
 *  Based on PhysioNet WAV reader by hvaanane
 */

#include <cstdio>
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "macro.h"
#include "../utils/memory_manager.h"

#include "PhysionetChallenge2016.h"

namespace Signal {

/**
 * Read PhysioNet / CinC 2016 Wav file to buffer
 *
 * Corrects corrupted data baseline and amplitude scale on demand
 *
 * @param Input filename
 */
PhysionetChallenge2016::PhysionetChallenge2016(char const *basename) {
	char filename[FILENAME_MAX];
	sprintf(filename, "%s.wav", basename);

	struct stat sb;
	if (stat(filename, &sb)) {
		PERROR("stat");
		throw errno;
	}

	FILE *f = fopen(filename, "rb");
	if (f == NULL) {
		PERROR("fopen");
		throw errno;
	}

 	// FIXME: supports only one channel little-endian 16bit PCM data
	dat.sample_freq = 2000.0;
	dat.samples_per_channel = (sb.st_size - 44) / sizeof(int16_t);
	dat.number_of_channels = 1;
	dat.ch = static_cast<struct data_channel *>(mm_malloc(
			sizeof(struct data_channel)));

	// read and convert
	int16_t *data_buffer = static_cast<int16_t *>(mm_malloc(
			dat.samples_per_channel * sizeof(int16_t)));
	fseek(f, 44, SEEK_SET);
	size_t nread = fread(data_buffer, sizeof(int16_t), dat.samples_per_channel,
			f);
	fclose(f);

	if (nread < dat.samples_per_channel) {
		errno = EIO;
		throw errno;
	}

	dat.ch->raw = static_cast<data_raw_t *>(mm_malloc(
			dat.samples_per_channel * sizeof(data_raw_t)));
	dat.ch->samples_per_mv = 1000;
	{
		data_raw_t *a = dat.ch->raw;
		int16_t *b = data_buffer;
		for (size_t i = 0; i < dat.samples_per_channel; ++i) {
			*a++ = *b++;
		}
	}
	mm_free(data_buffer);

	// correct average
	if (1) {
		double ave = 0.0;
		size_t start;
		if (dat.samples_per_channel > 4000) {
			start = 1000;
		}

		data_raw_t *a = dat.ch->raw;
		for (size_t i = start; i < dat.samples_per_channel; i++) {
			ave += *a++;
		}
		ave /= (double) (dat.samples_per_channel - start);

		a = dat.ch->raw;
		for (size_t i = 0; i < dat.samples_per_channel; i++) {
			*(a + i) -= ave;
		}
	}

	// correct amplitude scale
	if (1) {
		size_t start, end;
		double min, max, scale;
		double max_range = 2000.0;
		if (dat.samples_per_channel < 4000) {
			start = 0;
			end = dat.samples_per_channel;
		} else if (dat.samples_per_channel < 10000) {
			start = 1000;
			end = dat.samples_per_channel - 1000;
		} else {
			start = 2000;
			end = 9000;
		}
		min = max = dat.ch->raw[start];
		for (size_t i = start; i < end; i++) {
			data_raw_t &d = *(dat.ch->raw + i);

			if (min > d) {
				min = d;
			}
			if (max < d) {
				max = d;
			}
		}
		scale = max_range / (max - min);
		data_raw_t *a = dat.ch->raw;
		for (size_t i = 0; i < dat.samples_per_channel; i++) {
			*(a + i) *= scale;
		}
	}
}

PhysionetChallenge2016::~PhysionetChallenge2016() {
	mm_free(dat.ch->raw);
}

/**
 * @return Pointer to struct data_raw_t
 */
data_raw_t *PhysionetChallenge2016::get_signal() const {
	return dat.ch->raw;
}
/**
 * @return Number of samples in data
 */
size_t const &PhysionetChallenge2016::size() const {
	return dat.samples_per_channel;
}

} /* namespace Signal */
