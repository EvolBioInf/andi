/* @brief Here follows a simple implementation of the GNU function `strchrnul`.
 * Please check the gnulib manual for details.
 */
#include <string.h>

char *strchrnul(const char *s, int c){
	char *p = strchr(s,c);

	return p != NULL ? p : strchr(s, '\0');
}
