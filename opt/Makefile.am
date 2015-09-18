noinst_LIBRARIES= libcompat.a
libcompat_a_SOURCES= compat-string.h compat-stdlib.h

if !HAVE_STRCHRNUL
libcompat_a_SOURCES+= strchrnul.c
endif

if !HAVE_REALLOCARRAY
libcompat_a_SOURCES+= reallocarray.c
endif
