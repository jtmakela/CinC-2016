/*
 * PhysionetChannelge2016.h
 *
 *  Created on: Apr 28, 2016
 *      Author: jtmakela
 */

#ifndef ENTRY_SRC_SIMPLIFIED_PHYSIONETCHALLENGE2016_H_
#define ENTRY_SRC_SIMPLIFIED_PHYSIONETCHALLENGE2016_H_

#include <stdio.h>
#include <errno.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "macro.h"
#include "types.data.h"
#include "../utils/memory_manager.h"

namespace Signal {

class PhysionetChannelge2016 {
public:
	explicit PhysionetChannelge2016(char const *basename);
	virtual ~PhysionetChannelge2016();

	data_raw_t *get_signal() const;
	size_t const &size() const;
private:
	struct data dat;
};

} /* namespace Signal */

#endif /* ENTRY_SRC_SIMPLIFIED_PHYSIONETCHALLENGE2016_H_ */
