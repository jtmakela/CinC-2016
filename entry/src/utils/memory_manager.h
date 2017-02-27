/**
 * memory_manager.h
 *
 *  Created on: Apr 8, 2015
 *      Author: jtmakela
 */

#ifndef SRC_UTILS_MEMORY_MANAGER_H_
#define SRC_UTILS_MEMORY_MANAGER_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <errno.h>
#include <pthread.h>
#include <unistd.h>

#include "macro.h"

#ifdef USE_MEMORY_MANAGER
void *_mm_malloc(size_t size, char const * const file, size_t const line);
void _mm_free(void *ptr, char const * const file, size_t const line);
void *_mm_calloc(size_t nmemb, size_t size, char const * const file, size_t const line);
void *_mm_realloc(void *ptr, size_t size, char const * const file, size_t const line);
char *_mm_strdup(char const *ptr, char const * const file, size_t const line);
void memory_manager_atexit(void);

#define mm_malloc(_a)		(_mm_malloc((_a), __FILE__, __LINE__))
#define mm_calloc(_a, _b)	(_mm_calloc((_a), (_b), __FILE__, __LINE__))
#define mm_realloc(_a, _b)	(_mm_realloc((_a), (_b), __FILE__, __LINE__))
#define mm_free(_a)			(_mm_free((_a), __FILE__, __LINE__))

#define mm_strdup(_a)		(_mm_strdup((_a), __FILE__, __LINE__))

#else
#ifdef SHARED_OBJ
void *mm_malloc(size_t size);
void mm_free(void *ptr);
void *mm_calloc(size_t nmemb, size_t size);
void *mm_realloc(void *ptr, size_t size);
#else
#define mm_malloc(size)			(malloc(size))
#define mm_free(ptr)			(free(ptr))
#define mm_calloc(nmemb, size)	(calloc(nmemb, size))
#define mm_realloc(ptr, size)	(realloc(ptr, size))
#define mm_strdup(_a)			(strdup((_a)))
#endif /* SHARED_OBJ */
#endif /* USE_MEMORY_MANAGER */

#endif /* SRC_UTILS_MEMORY_MANAGER_H_ */
