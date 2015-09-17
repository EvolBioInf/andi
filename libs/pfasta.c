/*
 * Copyright (c) 2015, Fabian Klötzl <fabian-pfasta@kloetzl.info>
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
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#include "pfasta.h"

#define BUFFERSIZE 4096

/** @file
 *
 * Welcome to the code of `pfasta`, the pedantic FASTA parser. For future
 * reference I here explain some general, noteworthy things about the code.
 *
 *  - Most functions returning an `int` follow the zero-errors convention. On
 *    success a `0` is returned. A negative number indicates an error. Positive
 *    numbers can be used to signal a different exceptional state (i.e. EOF).
 *  - All functions use pointers to objects for their parameters.
 *  - The low-level Unix `read(2)` function is used to grab bytes from a file
 *    descriptor. These are stored in a buffer. The contents of this buffer is
 *    checked one char at a time and then appended to string buffer. Finally,
 *    the resulting strings are returned.
 *  - To work around the fact that C has no exceptions, I declare some nifty
 *    macros below.
 *  - As the length of the sequences are not known in advance I implemented a
 *    simple structure for growable strings called `dynstr`. One character at
 *    a time can be appended using `dynstr_put`. Internally an array is
 *    realloced with growth factor 1.5.
 *  - The functions which do the actual parsing are prefixed `pfasta_read_*`.
 */

#define PF_EXIT_ERRNO()                                                        \
	do {                                                                       \
		pf->errno__ = errno;                                                   \
		pf->errstr = NULL;                                                     \
		return errno;                                                          \
	} while (0)

#define PF_EXIT_FORWARD() return -1

#define PF_FAIL_FORWARD()                                                      \
	do {                                                                       \
		return_code = -1;                                                      \
		goto cleanup;                                                          \
	} while (0)

#define PF_FAIL_ERRNO()                                                        \
	do {                                                                       \
		pf->errno__ = errno;                                                   \
		pf->errstr = NULL;                                                     \
		return_code = errno;                                                   \
		goto cleanup;                                                          \
	} while (0)

#define PF_FAIL_STR(str)                                                       \
	do {                                                                       \
		pf->errno__ = 0;                                                       \
		pf->errstr = str;                                                      \
		return_code = -1;                                                      \
		goto cleanup;                                                          \
	} while (0)

static int buffer_init(pfasta_file *pf);
static inline int buffer_peek(const pfasta_file *pf);
static inline int buffer_adv(pfasta_file *pf);
static int buffer_read(pfasta_file *pf);

typedef struct dynstr {
	char *str;
	size_t capacity, count;
} dynstr;

static inline void *reallocarray(void *optr, size_t nmemb, size_t size);
static inline int dynstr_init(dynstr *ds);
static inline int dynstr_put(dynstr *ds, char c);
static inline void dynstr_free(dynstr *ds);
static inline char *dynstr_move(dynstr *ds);
static inline size_t dynstr_len(const dynstr *ds);

int pfasta_read_name(pfasta_file *pf, pfasta_seq *ps);
int pfasta_read_comment(pfasta_file *pf, pfasta_seq *ps);
int pfasta_read_seq(pfasta_file *pf, pfasta_seq *ps);

/*
 * When reading from a buffer, basically three things can happen.
 *
 *  1. Bytes are read (success)
 *  2. Low-level error (fail)
 *  3. No more bytes (EOF)
 *
 * A low-level error is indicated by `buffer_adv` returning non-zero. As
 * end-of-file is technically a successful read, EOF is instead signalled by
 * `buffer_peek`.
 */

/** @brief Initialises the read buffer. First, memory is allocated and then
 * filled. Both operations may fail. But this ensure we detect problems with
 * an unreadable file at the the initialisation of the parser!
 *
 * @param pf - The parser we want to initialise.
 * @returns 0 iff successful.
 */
static int buffer_init(pfasta_file *pf) {
	char *buffer = malloc(BUFFERSIZE);
	if (!buffer) PF_EXIT_ERRNO();

	pf->buffer = pf->readptr = pf->fillptr = buffer;
	if (buffer_read(pf) < 0) PF_EXIT_FORWARD();
	return 0;
}

/** @brief Returns the current character or EOF.
 *
 * @param pf - The parser to read from.
 * @returns The current character or EOF.
 */
static inline int buffer_peek(const pfasta_file *pf) {
	if (pf->readptr < pf->fillptr) {
		return (int)*(pf->readptr);
	}
	return EOF;
}

/** @brief Advances the read pointer in the buffer to the next character. If
 * needed, the buffer is filled with fresh bytes. A non-zero value is returned
 * if reading fails.
 *
 * @param pf - The parser which should be advanced.
 * @returns 0 iff successful.
 */
static inline int buffer_adv(pfasta_file *pf) {
	if (pf->readptr < pf->fillptr - 1) {
		pf->readptr++;
		return 0;
	}

	if (buffer_read(pf) < 0) PF_EXIT_FORWARD();

	return 0;
}

/** @brief Fills the buffer with new data.
 *
 * @param pf - The parser which should be updated.
 * @returns 0 iff successful.
 */
static int buffer_read(pfasta_file *pf) {
	ssize_t count = read(pf->fd, pf->buffer, BUFFERSIZE);
	if (count < 0) PF_EXIT_ERRNO();
	if (count == 0) { // EOF
		pf->fillptr = pf->buffer;
		pf->readptr = pf->buffer + 1;
		return 1;
	}

	pf->readptr = pf->buffer;
	pf->fillptr = pf->buffer + count;

	return 0;
}

/** @brief Frees all data associated with a parser. Also nulls pointers to avoid
 * a potential double-free.
 *
 * @param pf - The parser that shall be freed.
 */
void pfasta_free(pfasta_file *pf) {
	if (!pf) return;
	free(pf->buffer);
	pf->buffer = pf->readptr = pf->fillptr = pf->errstr = NULL;
	pf->errno__ = 0;
	pf->fd = -1;
}

/** @brief Creates a new parser for the given file. This includes allocating the
 * buffer and reading the first few bytes. These are then used to break on empty
 * or non-FASTA files.
 *
 * @param pf - A pointer to the parser structure we intend to initialise. No
 * assumption is made about the referenced memory except its existence.
 * @returns 0 iff successful.
 */
int pfasta_parse(pfasta_file *pf, int file_descriptor) {
	assert(pf && file_descriptor >= 0);
	int return_code = 0;

	pf->errno__ = 0;
	pf->buffer = pf->readptr = pf->fillptr = pf->errstr = NULL;
	pf->fd = file_descriptor;

	if (buffer_init(pf) != 0) PF_FAIL_FORWARD();

	int c = buffer_peek(pf);
	if (c == EOF) PF_FAIL_STR("Empty file");
	if (c != '>') PF_FAIL_STR("File does not start with '>'");

cleanup:
	return return_code;
}

/** @brief Frees the memory of a FASTA sequence. **/
void pfasta_seq_free(pfasta_seq *ps) {
	if (!ps) return;
	free(ps->name);
	free(ps->comment);
	free(ps->seq);
	ps->name = ps->comment = ps->seq = NULL;
}

/** @brief Reads the next sequence from the parser into the memory pointed to by
 * the parameter `ps`. This may fail for various reasons. No matter, what
 * happens, always free `ps` after usage!
 *
 * @param pf - The parser to read from.
 * @param ps - A reference to memory for the sequence data.
 *
 * @returns 0 if successful, 1 if the end of the file was reached and a negative
 * number on error.
 */
int pfasta_read(pfasta_file *pf, pfasta_seq *ps) {
	assert(pf && ps && pf->buffer);
	*ps = (pfasta_seq){NULL, NULL, NULL};
	int return_code = 0;

	if (buffer_peek(pf) == EOF) return 1;
	if (buffer_peek(pf) != '>') PF_FAIL_STR("Expected '>'");

	if (pfasta_read_name(pf, ps) < 0) PF_FAIL_FORWARD();
	if (isblank(buffer_peek(pf))) {
		if (pfasta_read_comment(pf, ps) < 0) PF_FAIL_FORWARD();
	}
	if (pfasta_read_seq(pf, ps) < 0) PF_FAIL_FORWARD();

	// Skip blank lines
	while (buffer_peek(pf) == '\n') {
		if (buffer_adv(pf) != 0) PF_FAIL_FORWARD();
	}

cleanup:
	return return_code;
}

/** @brief Reads the sequence name and saves it into the structure.
 *
 * @param pf - The parser used for reading.
 * @param ps - The structure used to hold the name, later.
 *
 * @returns 0 iff successful
 */
int pfasta_read_name(pfasta_file *pf, pfasta_seq *ps) {
	int return_code = 0;
	dynstr name;
	if (dynstr_init(&name) != 0) PF_FAIL_ERRNO();

	while (1) {
		if (buffer_adv(pf) != 0) PF_FAIL_FORWARD();

		int c = buffer_peek(pf);
		if (c == EOF) PF_FAIL_STR("Unexpected EOF in sequence name");
		if (!isgraph(c)) break;

		if (dynstr_put(&name, c) != 0) PF_FAIL_ERRNO();
	}

	if (dynstr_len(&name) == 0) PF_FAIL_STR("Empty name");

	ps->name = dynstr_move(&name);

cleanup:
	dynstr_free(&name);
	return return_code;
}

/** @brief Reads the sequence comment and saves it into the structure.
 *
 * @param pf - The parser used for reading.
 * @param ps - The structure used to hold the comment, later.
 *
 * @returns 0 iff successful
 */
int pfasta_read_comment(pfasta_file *pf, pfasta_seq *ps) {
	int return_code = 0;
	dynstr comment;
	if (dynstr_init(&comment) != 0) PF_FAIL_ERRNO();

	while (1) {
		if (buffer_adv(pf) != 0) PF_FAIL_FORWARD();

		int c = buffer_peek(pf);
		if (c == EOF) PF_FAIL_STR("Unexpected EOF in sequence comment");
		if (c == '\n') break;

		if (dynstr_put(&comment, c) != 0) PF_FAIL_ERRNO();
	}

	ps->comment = dynstr_move(&comment);

cleanup:
	dynstr_free(&comment);
	return return_code;
}

/** @brief Reads the sequence data and saves it into the structure.
 *
 * @param pf - The parser used for reading.
 * @param ps - The structure used to hold the data, later.
 *
 * @returns 0 iff successful
 */
int pfasta_read_seq(pfasta_file *pf, pfasta_seq *ps) {
	int return_code = 0;
	dynstr seq;
	if (dynstr_init(&seq) != 0) PF_FAIL_ERRNO();

	while (1) {
		assert(buffer_peek(pf) == '\n');

		// deal with the first character explicitly
		if (buffer_adv(pf) != 0) PF_FAIL_FORWARD();

		int c = buffer_peek(pf);
		if (c == EOF || c == '>' || c == '\n') break;

		goto regular;

		// read line
		while (1) {
			if (buffer_adv(pf) != 0) PF_FAIL_FORWARD();

			c = buffer_peek(pf);
			if (c == '\n') break;

		regular:
			if (!isgraph(c)) PF_FAIL_STR("Unexpected character in sequence");
			if (dynstr_put(&seq, c) != 0) PF_FAIL_ERRNO();
		}
	}

	if (dynstr_len(&seq) == 0) PF_FAIL_STR("Empty sequence");
	ps->seq = dynstr_move(&seq);

cleanup:
	dynstr_free(&seq);
	return return_code;
}

/** @brief Returns an explanatory string for encountered errors. */
const char *pfasta_strerror(const pfasta_file *pf) {
	if (!pf) return NULL;
	if (pf->errno__ == 0) {
		return pf->errstr;
	} else {
		return strerror(pf->errno__);
	}
}

/** @brief Creates a new string that can grow dynamically.
 *
 * @param ds - A reference to the dynstr container.
 *
 * @returns 0 iff successful.
 */
static inline int dynstr_init(dynstr *ds) {
	*ds = (dynstr){NULL, 0, 0};
	ds->str = malloc(61);
	if (!ds->str) return -1;
	ds->str[0] = '\0';
	ds->capacity = 61;
	return 0;
}

/** @brief A append a character to a string.
 *
 * @param ds - A reference to the dynstr container.
 * @param c - The new character.
 *
 * @returns 0 iff successful.
 */
static inline int dynstr_put(dynstr *ds, char c) {
	if (ds->count >= ds->capacity - 1) {
		char *neu = reallocarray(ds->str, ds->capacity / 2, 3);
		if (!neu) {
			dynstr_free(ds);
			return -1;
		}
		ds->str = neu;
		ds->capacity = (ds->capacity / 2) * 3;
	}

	ds->str[ds->count++] = c;
	return 0;
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
	char *out = ds->str;
	out[ds->count] = '\0';
	*ds = (dynstr){NULL, 0, 0};
	return out;
}

/** @brief Returns the current length of the dynamic string. */
static inline size_t dynstr_len(const dynstr *ds) { return ds->count; }

/***************
 * The following code was shamelessly stolen from the OpenBSD-libc.
 * It is released under MIT-style license.
 ***************/

/*
 * Copyright (c) 2008 Otto Moerbeek <otto@drijf.net>
 *
 * Permission to use, copy, modify, and distribute this software for any
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
 */

/*
 * This is sqrt(SIZE_MAX+1), as s1*s2 <= SIZE_MAX
 * if both s1 < MUL_NO_OVERFLOW and s2 < MUL_NO_OVERFLOW
 */
#define MUL_NO_OVERFLOW ((size_t)1 << (sizeof(size_t) * 4))

static inline void *reallocarray(void *optr, size_t nmemb, size_t size) {
	if ((nmemb >= MUL_NO_OVERFLOW || size >= MUL_NO_OVERFLOW) && nmemb > 0 &&
	    SIZE_MAX / nmemb < size) {
		errno = ENOMEM;
		return NULL;
	}
	return realloc(optr, size * nmemb);
}
