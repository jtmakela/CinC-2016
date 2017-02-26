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
 * Csv2kernel.cpp
 *
 *  Created on: Apr 11, 2016
 *      Author: jtmakela
 */

#include "Csv2kernel.h"

namespace Trigger {

Csv2kernel::Csv2kernel(char const *filename) :
		kernel(0), kernel_len(0) {
	struct stat sb;

	if (stat(filename, &sb)) {
		PERROR("stat");
		throw errno;
	}

	char cachefilename[FILENAME_MAX];
	sprintf(cachefilename, "/tmp/Csv2kernel.%lu.dat", sb.st_ino);
	if (stat(cachefilename, &sb) == 0) {
		kernel_len = sb.st_size / sizeof(cl_float);
		kernel = static_cast<cl_float *>(mm_malloc(
				kernel_len * sizeof(cl_float)));

		FILE *f = fopen(cachefilename, "r");
		if (f == 0) {
			PERROR("fopen");
			throw errno;
		}

		if (fread(kernel, sizeof(cl_float), kernel_len, f) != kernel_len) {
			PERROR("fread");
			throw errno;
		}

		fclose(f);
		return;
	}

	char *c = (char *) mm_malloc(sb.st_size);
	if (c == 0) {
		PERROR("malloc");
		throw errno;
	}

	FILE *f = fopen(filename, "r");
	if (f == 0) {
		PERROR("fopen");
		throw errno;
	}

	try {
		// count lines;
		for (;; ++kernel_len) {
			if (fgets(c, sb.st_size, f) == 0) {
				if (feof(f)) {
					break;
				}
				PERROR("fgets");
				throw errno;
			}
		}

		kernel = static_cast<cl_float *>(mm_calloc(kernel_len,
				sizeof(cl_float)));

		fseek(f, 0, SEEK_SET);

		for (cl_float *p = kernel;; ++p) {
			if (fgets(c, sb.st_size, f) == 0) {
				if (feof(f)) {
					break;
				}
				PERROR("fgets");
				throw errno;
			}

			*p = atof(c);
		}
	} catch (typeof(errno) e) {
		fclose(f);
		mm_free(c);
		throw e;
	}

	fclose(f);
	mm_free(c);

	f = fopen(cachefilename, "w");
	if (f == 0) {
		PERROR("fopen");
		throw errno;
	}

	if (fwrite(kernel, sizeof(cl_float), kernel_len, f) != kernel_len) {
		PERROR("fwrite");
		throw errno;
	}

	fclose(f);
}

Csv2kernel::~Csv2kernel() {
	if (kernel) {
		mm_free(kernel);
	}
}

cl_float const *Csv2kernel::get_data() const {
	return kernel;
}
size_t const &Csv2kernel::size() const {
	return kernel_len;
}

} /* namespace Trigger */
