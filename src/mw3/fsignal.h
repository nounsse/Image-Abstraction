/*
 * fsignal.h
 */

#ifndef _FSIGNAL_H_
#define _FSIGNAL_H_

#ifdef __cplusplus
extern "C" {
#endif


/* src/fsignal.c */
Fsignal mw_new_fsignal(void);
Fsignal mw_alloc_fsignal(Fsignal signal, int N);
void mw_delete_fsignal(Fsignal signal);
Fsignal mw_change_fsignal(Fsignal signal, int N);
void mw_clear_fsignal(Fsignal signal, float v);
void mw_copy_fsignal_values(Fsignal in, Fsignal out);
void mw_copy_fsignal_header(Fsignal in, Fsignal out);
void mw_copy_fsignal(Fsignal in, Fsignal out);

#ifdef __cplusplus
}
#endif


#endif /* !_FSIGNAL_H_ */
