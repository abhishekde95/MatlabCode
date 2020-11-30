#ifndef FIXSTIM_INCLUDED
#define FIXSTIM_INCLUDED

#define INT 0
#define LONG 1
#define DOUBLE 2
#define FLOAT 3
#define PARADIGMID 105

/* Codes that get dropped in the header (on the first trial only) */
#define FIXXCD          8000
#define FIXYCD          8001
#define FPSIZECD        8002
#define FPRCD           8003
#define FPGCD           8004
#define FPBCD           8005
#define EYEWINXCD       8006
#define EYEWINYCD       8007
#define TARGSIZECD      8008
#define TARGRCD         8009
#define TARGGCD         8010
#define TARGBCD         8011
#define TARGWINXCD      8012
#define TARGWINYCD      8013
#define RFXCD           8014
#define RFYCD           8015
#define NTRLSPERCONDCD  8016
#define BKGNDRGBCD      8017
#define GAMMATABLECD    8018
#define FUNDAMENTALSCD  8019
#define MONSPDCD        8020
#define MODULATORCD     8021
#define PARADIGMIDENTCD 8999

/* Codes that get dropped on every trial */
#define TARGSHOWNCD     7000
#define STIMTYPECD      7001
#define TTLSTIMFREQCD   7002
#define STIMFREQCD      7002
#define ELECSTIMFREQCD  7002
#define LASERSTIMFREQCD 7002
#define LASERPOWERCD    7003
#define CONTRASTCD      7004
#define STIMHEIGHTCD    7005
#define STIMWIDTHCD     7006
#define LCCCD           7007
#define MCCCD           7008
#define SCCCD           7009

typedef struct {
    int flagcodeval;
    int casttype;
    int nelements;
    union {
        long lval;
        float fval;
        double dval;
        long *lpval;
    } val;
} ecode_struct;

#define PARADIGM_PRELUDE( name )                                      \
do { static int _____do_once = 1;                                     \
     if ( _____do_once ) {                                            \
         dprintf(name " compiled on %s at %s\n", __DATE__, __TIME__); \
         _____do_once = 0; }}                                         \
while ( 0 )

#endif /* FIXSTIM_INCLUDED */
