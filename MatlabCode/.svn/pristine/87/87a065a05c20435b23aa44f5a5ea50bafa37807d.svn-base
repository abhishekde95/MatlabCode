/* WhiteNoise.h - WhiteNoise specific ecodes */
#ifndef __WHITENOISE_INCLUDED
#define __WHITENOISE_INCLUDED

/* Codes that get dropped in the header (on the first trial only) */
#define FIXXCD          8000
#define FIXYCD          8001
#define FPSIZECD        8002
#define FPRCD           8003
#define FPGCD           8004
#define FPBCD           8005
#define EYEWINXCD       8006
#define EYEWINYCD       8007
#define RFXCD           8008
#define RFYCD           8009
#define NPIXPERSTIXCD   8010
#define NSTIXPERSIDECD  8011
#define GAMMATABLECD    8012
#define MONSPDCD        8013
#define FRAMERATECD     8014
#define GAUSSLOCUTCD    8015
#define GAUSSHICUTCD    8016
#define FUNDAMENTALSCD  8017
#define PARADIGMIDENTCD 8999

/* Codes that get dropped on every trial */
#define SEEDCD       7000
#define MU1CD        7001
#define MU2CD        7002
#define MU3CD        7003
#define SIGMA1CD     7004
#define SIGMA2CD     7005
#define SIGMA3CD     7006
#define NFRAMESCD    7007
#define BKGNDRCD     7008
#define BKGNDGCD     7009
#define BKGNDBCD     7010
#define SYNTHIMAGECD 7011
#define NOISETYPECD  7012
#define MASKCD       7013

typedef struct {
    int flagcodeval;
    int casttype;
    int nelements;
    union {
        long lval;
        double dval;
        long *lpval;
    } val;
} ecode_struct;

enum {INTERLEAVE, GUN, CONE, RGCONE}; // noise types

#define INT 0
#define LONG 1
#define DOUBLE 2
#define LONGPOINTER 3
#define PARADIGMID 100

#endif /* __WHITENOISE_INCLUDED */
