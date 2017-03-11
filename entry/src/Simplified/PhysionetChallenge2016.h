/**
 * PhysionetChallenge2016.h
 *
 *  Created on: Apr 28, 2016
 *      Author: jtmakela
 */

#ifndef ENTRY_SRC_SIMPLIFIED_PHYSIONETCHALLENGE2016_H_
#define ENTRY_SRC_SIMPLIFIED_PHYSIONETCHALLENGE2016_H_

#include "types.data.h"

namespace Signal {

class PhysionetChallenge2016 {
public:
	explicit PhysionetChallenge2016(char const *basename);
	virtual ~PhysionetChallenge2016();

	data_raw_t *get_signal() const;
	size_t const &size() const;
private:
	struct data dat;
};

} /* namespace Signal */

#endif /* ENTRY_SRC_SIMPLIFIED_PHYSIONETCHALLENGE2016_H_ */
