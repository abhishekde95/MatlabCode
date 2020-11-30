#ifndef CHROMCAT_INCLUDED
#define CHROMCAT_INCLUDED

#define INT 0
#define LONG 1
#define DOUBLE 2
#define FLOAT 3
#define PARADIGMID 158

// header codes
#define FIXXCD     8000
#define FIXYCD     8001
#define FPSIZECD   8002
#define EYEWINWCD  8003
#define EYEWINHCD  8004
#define TARGXCD    8005
#define TARGYCD    8006
#define TARGSIZECD 8007
#define TARGWINWCD 8008
#define TARGWINHCD 8009
#define STIMDURCD  8010
#define NBLOCKSCD  8013
#define MMATRIXCD  8014
#define BKGNDRGBCD 8015

#define PARADIGMIDENTCD 8999

// trial codes
#define CORRECTTRIALCD  7001
#define THRESHMULTCD    7002
#define LCCCD           7003
#define MCCCD           7004
#define CORRTARGXCD     7005
#define CORRTARGYCD     7006
#define RIGHTWARDSACCCD 7007
#define TFCD            7008

typedef struct {
    int flagcodeval;
    int casttype;
    int nelements;
    union {
        long lval;
        float fval;
        double dval;
        long* lpval;
    } val;
} ecode_struct;

void updateEyewin(int, int, int, int, int);

#define PARADIGM_PRELUDE( name ) \
do { \
    static int _____do_once = 1; \
    if ( _____do_once ) { \
        dprintf(name " compiled on %s at %s\n", __DATE__, __TIME__); \
        _____do_once = 0; \
    } \
} \
while ( 0 )

#endif
