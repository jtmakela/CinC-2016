/*
 * macro.h
 *
 *  Created on: Feb 26, 2016
 *      Author: jtmakela
 */

#ifndef INCLUDE_MACRO_H_
#define INCLUDE_MACRO_H_

#include <cstdio>
#include <cstring>
#include <errno.h>

#define PERROR(_e) { fprintf(stderr, "[%s:%u] %s: %s\n", __FILE__, __LINE__, _e, strerror(errno)); throw errno; }

#endif /* INCLUDE_MACRO_H_ */
