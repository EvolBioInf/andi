/*
 * Copyright (c) 2015-2020, Fabian Klötzl <fabian-pfasta@kloetzl.info>
 *
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 *
 */

#include <assert.h>
#include <ctype.h>
#include <err.h>
#include <errno.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "pfasta.h"

#define VERSION "v15"

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __STDC_NO_THREADS__
#define thread_local
#else
#include <threads.h>
#endif

/** The following is the maximum length of an error string. It has to be
 * carefully chosen, so that all calls to PF_FAIL_STR succeed. For instance,
 * the line number can account for up to 20 characters.
 */
#define PF_ERROR_STRING_LENGTH 100

thread_local char errstr_buffer[PF_ERROR_STRING_LENGTH];

void *pfasta_reallocarray(void *ptr, size_t nmemb, size_t size);

#define BUFFER_SIZE 16384

#define LIKELY(X) __builtin_expect((intptr_t)(X), 1)
#define UNLIKELY(X) __builtin_expect((intptr_t)(X), 0)

enum { NO_ERROR, E_EOF, E_ERROR, E_ERRNO, E_BUBBLE, E_STR, E_STR_CONST };

#define PF_FAIL_ERRNO(PP)                                                      \
	do {                                                                       \
		(void)strerror_r(errno, errstr_buffer, PF_ERROR_STRING_LENGTH);        \
		(PP)->errstr = errstr_buffer;                                          \
		return_code = E_ERRNO;                                                 \
		goto cleanup;                                                          \
	} while (0)

#define PF_FAIL_BUBBLE_CHECK(PP, CHECK)                                        \
	do {                                                                       \
		if (UNLIKELY(CHECK)) {                                                 \
			return_code = CHECK;                                               \
			goto cleanup;                                                      \
		}                                                                      \
	} while (0)

#define PF_FAIL_BUBBLE(PP)                                                     \
	do {                                                                       \
		if (UNLIKELY((PP)->errstr)) {                                          \
			return_code = E_BUBBLE;                                            \
			goto cleanup;                                                      \
		}                                                                      \
	} while (0)

#define PF_FAIL_STR_CONST(PP, STR)                                             \
	do {                                                                       \
		(PP)->errstr = (STR);                                                  \
		return_code = E_STR_CONST;                                             \
		goto cleanup;                                                          \
	} while (0)

#define PF_FAIL_STR(PP, ...)                                                   \
	do {                                                                       \
		(void)snprintf(errstr_buffer, PF_ERROR_STRING_LENGTH, __VA_ARGS__);    \
		(PP)->errstr = errstr_buffer;                                          \
		return_code = E_STR;                                                   \
		goto cleanup;                                                          \
	} while (0)

int pfasta_read_name(struct pfasta_parser *pp, struct pfasta_record *pr);
int pfasta_read_comment(struct pfasta_parser *pp, struct pfasta_record *pr);
int pfasta_read_sequence(struct pfasta_parser *pp, struct pfasta_record *pr);

static inline char *buffer_begin(struct pfasta_parser *pp);
static inline char *buffer_end(struct pfasta_parser *pp);
static inline int buffer_advance(struct pfasta_parser *pp, size_t steps);
static inline int buffer_is_empty(const struct pfasta_parser *pp);
static inline int buffer_is_eof(const struct pfasta_parser *pp);
static inline int buffer_peek(struct pfasta_parser *pp);
static inline int buffer_read(struct pfasta_parser *pp);

typedef struct dynstr {
	char *str;
	size_t capacity, count;
} dynstr;

static inline char *dynstr_move(dynstr *ds);
static inline int dynstr_init(dynstr *ds, struct pfasta_parser *pp);
static inline size_t dynstr_len(const dynstr *ds);
static inline void dynstr_free(dynstr *ds);
static inline int dynstr_append(dynstr *ds, const char *str, size_t length,
                                struct pfasta_parser *pp);

static inline int my_isspace(int c) {
	// ascii whitespace
	return (c >= '\t' && c <= '\r') || (c == ' ');
}

const char *pfasta_version(void) { return VERSION; }

int buffer_init(struct pfasta_parser *pp) {
	int return_code = 0;

	pp->buffer = malloc(BUFFER_SIZE);
	if (!pp->buffer) PF_FAIL_ERRNO(pp);

	int check = buffer_read(pp);
	PF_FAIL_BUBBLE_CHECK(pp, check);

cleanup:
	return return_code;
}

int buffer_read(struct pfasta_parser *pp) {
	int return_code = NO_ERROR;
	ssize_t count = read(pp->file_descriptor, pp->buffer, BUFFER_SIZE);

	if (UNLIKELY(count < 0)) PF_FAIL_ERRNO(pp);
	if (UNLIKELY(count == 0)) { // EOF
		pp->fill_ptr = pp->buffer;
		pp->read_ptr = pp->buffer + 1;
		pp->errstr = "EOF (maybe error)"; // enable bubbling
		return E_EOF;
	}

	pp->read_ptr = pp->buffer;
	pp->fill_ptr = pp->buffer + count;

cleanup:
	return return_code;
}

int buffer_peek(struct pfasta_parser *pp) {
	return LIKELY(pp->read_ptr < pp->fill_ptr) ? *(unsigned char *)pp->read_ptr
	                                           : EOF;
}

char *buffer_begin(struct pfasta_parser *pp) { return pp->read_ptr; }

char *buffer_end(struct pfasta_parser *pp) { return pp->fill_ptr; }

inline int buffer_advance(struct pfasta_parser *pp, size_t steps) {
	int return_code = 0;

	pp->read_ptr += steps;
	if (UNLIKELY(pp->read_ptr >= pp->fill_ptr)) {
		assert(pp->read_ptr == pp->fill_ptr);
		int check = buffer_read(pp); // resets pointers
		PF_FAIL_BUBBLE_CHECK(pp, check);
	}

cleanup:
	return return_code;
}

int buffer_is_empty(const struct pfasta_parser *pp) {
	return pp->read_ptr == pp->fill_ptr;
}

int buffer_is_eof(const struct pfasta_parser *pp) {
	return pp->read_ptr > pp->fill_ptr;
}

char *find_first_space(const char *begin, const char *end) {
	size_t offset = 0;
	size_t length = end - begin;

#ifdef __SSE2__

	typedef __m128i vec_type;
	static const size_t vec_size = sizeof(vec_type);

	const vec_type all_tab = _mm_set1_epi8('\t' - 1);
	const vec_type all_carriage = _mm_set1_epi8('\r' + 1);
	const vec_type all_space = _mm_set1_epi8(' ');

	size_t vec_offset = 0;
	size_t vec_length = (end - begin) / vec_size;

	for (; vec_offset < vec_length; vec_offset++) {
		vec_type chunk;
		memcpy(&chunk, begin + vec_offset * vec_size, vec_size);

		// isspace: \t <= char <= \r || char == space
		vec_type v1 = _mm_cmplt_epi8(all_tab, chunk);
		vec_type v2 = _mm_cmplt_epi8(chunk, all_carriage);
		vec_type v3 = _mm_cmpeq_epi8(chunk, all_space);

		unsigned int vmask = (_mm_movemask_epi8(v1) & _mm_movemask_epi8(v2)) |
		                     _mm_movemask_epi8(v3);

		if (UNLIKELY(vmask)) {
			offset += __builtin_ctz(vmask);
			offset += vec_offset * vec_size;
			return (char *)begin + offset;
		}
	}

	offset += vec_offset * vec_size;
#endif

	for (; offset < length; offset++) {
		if (my_isspace(begin[offset])) break;
	}
	return (char *)begin + offset;
}

char *find_first_not_space(const char *begin, const char *end) {
	size_t offset = 0;
	size_t length = end - begin;

	for (; offset < length; offset++) {
		if (!my_isspace(begin[offset])) break;
	}
	return (char *)begin + offset;
}

size_t count_newlines(const char *begin, const char *end) {
	size_t offset = 0;
	size_t length = end - begin;
	size_t newlines = 0;

	for (; offset < length; offset++) {
		if (begin[offset] == '\n') newlines++;
	}

	return newlines;
}

static int copy_word(struct pfasta_parser *pp, dynstr *target) {
	int return_code = 0;

	int c;
	while (c = buffer_peek(pp), c != EOF && LIKELY(!my_isspace(c))) {
		char *end_of_word = find_first_space(buffer_begin(pp), buffer_end(pp));
		size_t word_length = end_of_word - buffer_begin(pp);

		assert(word_length > 0);

		int check = dynstr_append(target, buffer_begin(pp), word_length, pp);
		PF_FAIL_BUBBLE_CHECK(pp, check);

		check = buffer_advance(pp, word_length);
		PF_FAIL_BUBBLE_CHECK(pp, check);
	}

cleanup:
	return return_code;
}

static int skip_whitespace(struct pfasta_parser *pp) {
	int return_code = 0;

	while (my_isspace(buffer_peek(pp))) {
		char *split = find_first_not_space(buffer_begin(pp), buffer_end(pp));

		// advance may clear the buffer. So count first …
		size_t newlines = count_newlines(buffer_begin(pp), split);
		int check = buffer_advance(pp, split - buffer_begin(pp));
		PF_FAIL_BUBBLE_CHECK(pp, check);

		// … and then increase the counter.
		pp->line_number += newlines;
	}

cleanup:
	return return_code;
}

struct pfasta_parser pfasta_init(int file_descriptor) {
	int return_code = 0;
	struct pfasta_parser pp = {0};
	pp.line_number = 1;

	pp.file_descriptor = file_descriptor;
	int check = buffer_init(&pp);
	if (check && check != E_EOF) PF_FAIL_BUBBLE_CHECK(&pp, check);

	if (buffer_is_empty(&pp) || buffer_is_eof(&pp)) {
		PF_FAIL_STR(&pp, "File is empty.");
	}

	if (buffer_peek(&pp) != '>') {
		PF_FAIL_STR(&pp, "File must start with '>'.");
	}

cleanup:
	// free buffer if necessary
	if (return_code) {
		pfasta_free(&pp);
	}
	pp.done = return_code || buffer_is_eof(&pp);
	return pp;
}

struct pfasta_record pfasta_read(struct pfasta_parser *pp) {
	int return_code = 0;
	struct pfasta_record pr = {0};

	int check = pfasta_read_name(pp, &pr);
	PF_FAIL_BUBBLE_CHECK(pp, check);

	check = pfasta_read_comment(pp, &pr);
	PF_FAIL_BUBBLE_CHECK(pp, check);

	check = pfasta_read_sequence(pp, &pr);
	PF_FAIL_BUBBLE_CHECK(pp, check);

cleanup:
	if (return_code) {
		pfasta_record_free(&pr);
		pfasta_free(pp);
	}
	pp->done = return_code || buffer_is_eof(pp);
	return pr;
}

int pfasta_read_name(struct pfasta_parser *pp, struct pfasta_record *pr) {
	int return_code = 0;

	dynstr name;
	dynstr_init(&name, pp);
	PF_FAIL_BUBBLE(pp);

	assert(!buffer_is_empty(pp));
	if (buffer_peek(pp) != '>') {
		PF_FAIL_STR(pp, "Expected '>' but found '%c' on line %zu.",
		            buffer_peek(pp), pp->line_number);
	}

	int check = buffer_advance(pp, 1); // skip >
	if (check == E_EOF)
		PF_FAIL_STR(pp, "Unexpected EOF in name on line %zu.", pp->line_number);
	PF_FAIL_BUBBLE(pp);

	check = copy_word(pp, &name);
	if (check == E_EOF)
		PF_FAIL_STR(pp, "Unexpected EOF in name on line %zu.", pp->line_number);
	PF_FAIL_BUBBLE(pp);

	if (dynstr_len(&name) == 0)
		PF_FAIL_STR(pp, "Empty name on line %zu.", pp->line_number);

	pr->name_length = dynstr_len(&name);
	pr->name = dynstr_move(&name);

cleanup:
	if (return_code) {
		dynstr_free(&name);
	}
	return return_code;
}

int pfasta_read_comment(struct pfasta_parser *pp, struct pfasta_record *pr) {
	int return_code = 0;

	if (buffer_peek(pp) == '\n') {
		pr->comment_length = 0;
		pr->comment = NULL;
		return 0;
	}

	dynstr comment;
	dynstr_init(&comment, pp);
	PF_FAIL_BUBBLE(pp);

	assert(!buffer_is_empty(pp));

	int check = buffer_advance(pp, 1); // skip first whitespace
	if (check == E_EOF) goto label_eof;
	PF_FAIL_BUBBLE(pp);

	assert(!buffer_is_empty(pp));

	// get comment
	while (buffer_peek(pp) != '\n') {
		check = dynstr_append(&comment, buffer_begin(pp), 1, pp);
		PF_FAIL_BUBBLE_CHECK(pp, check);

		check = buffer_advance(pp, 1);
		if (check == E_EOF) goto label_eof;
		PF_FAIL_BUBBLE_CHECK(pp, check);
	}

label_eof:
	if (buffer_is_eof(pp))
		PF_FAIL_STR(pp, "Unexpected EOF in comment on line %zu.",
		            pp->line_number);

	pr->comment_length = dynstr_len(&comment);
	pr->comment = dynstr_move(&comment);

cleanup:
	if (return_code) {
		dynstr_free(&comment);
	}
	return return_code;
}

int pfasta_read_sequence(struct pfasta_parser *pp, struct pfasta_record *pr) {
	int return_code = 0;

	dynstr sequence;
	dynstr_init(&sequence, pp);
	PF_FAIL_BUBBLE(pp);

	assert(!buffer_is_empty(pp));
	assert(!buffer_is_eof(pp));
	assert(buffer_peek(pp) == '\n');

	int check = skip_whitespace(pp);
	if (check == E_EOF)
		PF_FAIL_STR(pp, "Empty sequence on line %zu.", pp->line_number);
	PF_FAIL_BUBBLE_CHECK(pp, check);

	// Assume a line begins only with alpha, -, *, or more spaces
	char c;
	while (c = buffer_peek(pp), LIKELY(isalpha(c) || c == '-' || c == '*')) {
		int check = copy_word(pp, &sequence);
		if (UNLIKELY(check == E_EOF)) break;
		PF_FAIL_BUBBLE_CHECK(pp, check);

		// optimize for more common case
		ptrdiff_t length = buffer_end(pp) - buffer_begin(pp);
		if (LIKELY(length >= 2 && buffer_begin(pp)[0] == '\n' &&
		           buffer_begin(pp)[1] > ' ')) {
			pp->read_ptr++; // nasty hack
			pp->line_number += 1;
		} else {
			check = skip_whitespace(pp);
			if (UNLIKELY(check == E_EOF)) break;
			PF_FAIL_BUBBLE_CHECK(pp, check);
		}
	}

	if (dynstr_len(&sequence) == 0)
		PF_FAIL_STR(pp, "Empty sequence on line %zu.", pp->line_number);

	pr->sequence_length = dynstr_len(&sequence);
	pr->sequence = dynstr_move(&sequence);
	pp->errstr = NULL; // reset error

cleanup:
	if (return_code) {
		dynstr_free(&sequence);
	}
	return return_code;
}

void pfasta_record_free(struct pfasta_record *pr) {
	if (!pr) return;
	free(pr->name);
	free(pr->comment);
	free(pr->sequence);
	pr->name = pr->comment = pr->sequence = NULL;
}

void pfasta_free(struct pfasta_parser *pp) {
	if (!pp) return;
	free(pp->buffer);
	pp->buffer = NULL;
}

/** @brief Creates a new string that can grow dynamically.
 *
 * @param ds - A reference to the dynstr container.
 *
 * @returns 0 iff successful.
 */
static inline int dynstr_init(dynstr *ds, struct pfasta_parser *pp) {
	int return_code = 0;

	*ds = (dynstr){NULL, 0, 0};
	ds->str = malloc(61);
	if (!ds->str) PF_FAIL_ERRNO(pp);

	ds->str[0] = '\0';
	ds->capacity = 61;
	ds->count = 0;

cleanup:
	return return_code;
}

/** @brief A append more than one character to a string.
 *
 * @param ds - A reference to the dynstr container.
 * @param str - The new characters.
 * @param length - number of new characters to append
 *
 * @returns 0 iff successful.
 */
static inline int dynstr_append(dynstr *ds, const char *str, size_t length,
                                struct pfasta_parser *pp) {
	int return_code = 0;
	size_t required = ds->count + length;

	if (UNLIKELY(required >= ds->capacity)) {
		char *neu = pfasta_reallocarray(ds->str, required / 2, 3);
		if (UNLIKELY(!neu)) {
			dynstr_free(ds);
			PF_FAIL_ERRNO(pp);
		}
		ds->str = neu;
		ds->capacity = (required / 2) * 3;
	}

	memcpy(ds->str + ds->count, str, length);
	ds->count = required;

cleanup:
	return return_code;
}

/** @brief Frees a dynamic string. */
static inline void dynstr_free(dynstr *ds) {
	if (!ds) return;
	free(ds->str);
	*ds = (dynstr){NULL, 0, 0};
}

/** @brief Returns the string as a standard `char*`. The internal reference is
 * then deleted. Hence the name *move* as in *move semantics*.
 *
 * @param ds - The dynamic string to move from.
 *
 * @returns a `char*` to a standard null-terminated string.
 */
static inline char *dynstr_move(dynstr *ds) {
	char *out = pfasta_reallocarray(ds->str, ds->count + 1, 1);
	if (!out) {
		out = ds->str;
	}
	out[ds->count] = '\0';
	*ds = (dynstr){NULL, 0, 0};
	return out;
}

/** @brief Returns the current length of the dynamic string. */
static inline size_t dynstr_len(const dynstr *ds) { return ds->count; }

__attribute__((weak)) void *reallocarray(void *ptr, size_t nmemb, size_t size);

/**
 * @brief Unsafe fallback in case reallocarray isn't provided by the stdlib.
 */
void *pfasta_reallocarray(void *ptr, size_t nmemb, size_t size) {
	if (reallocarray == NULL) {
		return realloc(ptr, nmemb * size);
	} else {
		return reallocarray(ptr, nmemb, size);
	}
}
