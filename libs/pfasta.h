/*
 * Copyright (c) 2015-2018, Fabian Klötzl <fabian-pfasta@kloetzl.info>
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

#ifndef PFASTA_H
#define PFASTA_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

/**
 * There is no magic to this structure. Its just a container of three strings.
 * Feel free to duplicate or move them. But don't forget to free the data after
 * usage!
 */
struct pfasta_record {
	char *name, *comment, *sequence;
	size_t name_length, comment_length, sequence_length;
};

/**
 * This structure holds a number of members to represent the state of the FASTA
 * parser. Please make sure that it is properly initialized before usage.
 * Always free this structure when the parser is done.
 */
struct pfasta_parser {
	const char *errstr;
	int done;

	/*< private -- do not touch! >*/
	int file_descriptor;
	char *buffer;
	char *read_ptr, *fill_ptr;
	size_t line_number;
};

/**
 * This function initializes a `pfasta_parser` struct with a parser bound to a
 * specific file descriptor. Iff an error occurred `errstr` is set to contain a
 * suitable message. Otherwise you can read data from it as long as `done` isn't
 * set. The parser should be freed after usage.
 *
 * Please note that the user is responsible for opening the file descriptor as
 * readable and closing after usage.
 */
struct pfasta_parser pfasta_init(int file_descriptor);

/**
 * Using a properly initialized parser, this function can read FASTA sequences.
 * These are stored in the simple structure and returned. On error, the `errstr`
 * property of the parser is set.
 */
struct pfasta_record pfasta_read(struct pfasta_parser *pp);

/**
 * This function frees the resources held by a pfasta record.
 */
void pfasta_record_free(struct pfasta_record *pr);

/**
 * This function frees the resources held by a pfasta parser.
 */
void pfasta_free(struct pfasta_parser *pp);

/**
 * Get a string defining the version of the pfasta library.
 */
const char *pfasta_version(void);

#ifdef __STDC_NO_THREADS__
/** If the preprocessor macro `PFASTA_NO_THREADS` is defined, the parser is not
 * fully thread safe. */
#define PFASTA_NO_THREADS
#endif

#ifdef __cplusplus
}
#endif

#endif /* PFASTA_H */
