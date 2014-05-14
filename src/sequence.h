#ifndef _SEQUENCE_H_
#define _SEQUENCE_H_

typedef struct {
	char *S;
	char *RS;
	size_t len;
	size_t RSlen;
	char *name;
} seq_t;

void freeSeq( seq_t *S);
char *revcomp( const char *str, size_t len);
char *catcomp( char *s , size_t len);

#endif

