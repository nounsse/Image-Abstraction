/*
 * error.h
 */

#ifndef _ERROR_H_
#define _ERROR_H_


#ifdef __cplusplus
extern "C" {
#endif



/* src/error.c */
void mwdebug(char *fmt, ...);
void mwerror(int code, int exit_code, char *fmt, ...);

#ifdef __cplusplus
}
#endif


#endif /* !_ERROR_H_ */
