/**
 * Csv2kernel.h
 *
 *  Created on: Apr 11, 2016
 *      Author: jtmakela
 */

#ifndef TRIGGER_CSV2KERNEL_H_
#define TRIGGER_CSV2KERNEL_H_

#include "types.data.h"

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
