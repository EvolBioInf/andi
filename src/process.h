#ifndef _PROCESS_H_
#define _PROCESS_H_

#include "sequence.h"

void printDistMatrix( seq_t* sequences, int n);

char *revcomp( const char *str, size_t len);
char *catcomp( char *s , size_t len);

#endif


