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
 * main.cpp
 *
 *  Created on: Mar 23, 2016
 *      Author: jtmakela
 */

#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>
#include <sys/stat.h>
#include <time.h>

#include "macro.h"

#include "Trigger/Csv2kernel.h"
#include "Trigger/Trigger.h"
#include "Simplified/Retrigger.h"
#include "Simplified/PhysionetChallenge2016.h"
#include "myDSP/Iir.h"
#include "Classifier/classifier.h"

int main(int argc, char *argv[]) {
	if (argc < 3) {
		fprintf(stderr,
				"[%s:%u] usage: %s <trigger convolution kernel.csv> <file base id, eg. a0123>\n",
				__FILE__, __LINE__, argv[0]);
		exit (EXIT_FAILURE);
	}

	char const *convolution_kernel = argv[1];
	char const *data_filename = argv[2];

	clock_t ref = clock();

	Signal::PhysionetChallenge2016 dat(data_filename);
	printf("data: %s\n", data_filename);

	Classifier::result_e result = Classifier::unknown;
	Simplified::Retrigger retrig(0.25);

	// A bit inconsistent naming convention.. these are required to get the energy signal
	Trigger::Csv2kernel kernel(convolution_kernel);
	retrig.set_convolution_kernel(kernel.get_data(), kernel.size());

	retrig.set_data(dat.get_signal(), dat.size());
	cl_float const *energy = retrig.get_energy();

	// A simple minmaxminmax trigger, will fire both on S1 and S2
	Accbpm::Trigger trig(energy, dat.size(), 2000.0);
	try {
		ref_ev const &ev = trig.get_events();
		printf("trigged %lu events\n", ev.size());

		retrig.set_ref_ev(ev);
		retrig.set_correlation_window(0.25, 0.125);
		retrig.set_lookaround_window(0.05, 0.025);

		// auto correlate, jut for fun!
		retrig.calc_correlations(0.8);
	} catch (int e) {
		exit (EXIT_FAILURE);
	}

	{
		struct data data = { };
		data.sample_freq = 2000.0;
		data.number_of_channels = 1;
		data.samples_per_channel = dat.size();

		struct data_channel data_ch[data.number_of_channels];
		memset(&data_ch, 0, sizeof(data_ch));
		data_ch[0].raw = dat.get_signal();

		data.ch = data_ch;
		try {
			Simplified::retrig_ev const &ev1 = retrig.get_s1_events();
			Simplified::retrig_ev const &ev2 = retrig.get_s2_events();
			printf(
					"got %lu S1 events and %lu S2 events after auto correlation\n",
					ev1.size(), ev2.size());

			result = Classifier::classify_this(&data, &ev1, &ev2,
					"params/s1s2.txt");
		} catch (int e) {
			Simplified::retrig_ev const &ev = retrig.get_events();
			printf("got %lu events after auto correlation\n", ev.size());

			if (ev.size()) {
				result = Classifier::classify_this(&data, &ev, 0,
						"params/ev.txt");
			} else {
				result = Classifier::classify_this(&data, 0, 0,
						"params/rest.txt");
			}

		}
	}

	{
		double t = static_cast<double>(clock() - ref) / CLOCKS_PER_SEC;
		printf("total execution time \e[38;5;%im%.3f\e[0m s\n",
				t < 9 ? 82 + static_cast<int>(t) : 196, t);
	}

	{
#define RESULT_CASE(_a) case _a: printf(" (" #_a ")\n"); break;
		printf("verdict: %i", result);
		switch (result) {
		RESULT_CASE(Classifier::unknown)
		RESULT_CASE(Classifier::normal)
		RESULT_CASE(Classifier::abnormal)
		}
	}

	// write result
	{
		FILE *f = fopen("answers.txt", "a");
		if (f == 0) {
			PERROR("fopen");
			return errno;
		}

		int r = 0;
		switch (result) {
		case Classifier::abnormal:
			r = 1;
			break;
		case Classifier::normal:
			r = -1;
			break;
		}

		fprintf(f, "%s,%i\n", data_filename, r);

		fclose(f);
	}
	return 0;
}
