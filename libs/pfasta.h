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

#ifndef PFASTA_H
#define PFASTA_H

/** The following is the maximum length of an error string. It has to be
 * carefully chosen, so that all calls to PF_FAIL_STR succeed. For instance,
 * the line number can account for up to 20 characters.
 */
#define PF_ERROR_STRING_LENGTH 100

typedef struct pfasta_file {
	char *buffer, *readptr, *fillptr;
	char *errstr;
	int errno__;
	int fd;
	size_t line;
	char errstr_buf[PF_ERROR_STRING_LENGTH];
	char unexpected_char;
} pfasta_file;

typedef struct pfasta_seq {
	char *name, *comment, *seq;
} pfasta_seq;

int pfasta_parse(pfasta_file *, int file_descriptor);
void pfasta_free(pfasta_file *);
void pfasta_seq_free(pfasta_seq *);
int pfasta_read(pfasta_file *, pfasta_seq *);

const char *pfasta_strerror(const pfasta_file *);

#endif /* PFASTA_H */
