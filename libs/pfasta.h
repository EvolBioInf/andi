#ifndef PFASTA_H
#define PFASTA_H

typedef struct pfasta_file {
	char *buffer, *readptr, *fillptr;
	char *errstr;
	int errno__;
	int fd;
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
