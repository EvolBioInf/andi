#ifndef _GLOBAL_H_
#define _GLOBAL_H_


extern int FLAGS;
extern int CORES;
extern int STRATEGY;

enum {
	F_NONE = 0,
	F_RAW = 1,
	F_SINGLE = 2,
	F_VERBOSE = 4,
	F_EXTRA_VERBOSE = 8
};

enum {
	S_SIMPLE,
	S_INC,
	S_WINDOW,
	S_CHECK,
	S_SOPHISTICATED
};

#endif

