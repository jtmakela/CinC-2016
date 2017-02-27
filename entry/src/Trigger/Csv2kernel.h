/**
 * Csv2kernel.h
 *
 *  Created on: Apr 11, 2016
 *      Author: jtmakela
 */

#ifndef TRIGGER_CSV2KERNEL_H_
#define TRIGGER_CSV2KERNEL_H_

#include <sys/types.h>
#include <sys/stat.h>
#include <libgen.h>
#include <unistd.h>
#include "macro.h"
#include "types.event.h"
#include "../utils/memory_manager.h"

namespace Trigger {

class Csv2kernel {
public:
	Csv2kernel(char const *filename);
	virtual ~Csv2kernel();

	cl_float const *get_data() const;
	size_t const &size() const;

private:
	cl_float *kernel;
	size_t kernel_len;
};

} /* namespace Trigger */

#endif /* TRIGGER_CSV2KERNEL_H_ */
