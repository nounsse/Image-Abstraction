/*
 * ccmovie.h
 */

#ifndef _CCMOVIE_H_
#define _CCMOVIE_H_

#ifdef __cplusplus
extern "C" {
#endif

/* src/ccmovie.c */
Ccmovie mw_new_ccmovie();
void mw_delete_ccmovie(Ccmovie movie);
Ccmovie mw_change_ccmovie(Ccmovie movie);

#ifdef __cplusplus
}
#endif


#endif /* !_CCMOVIE_H_ */
