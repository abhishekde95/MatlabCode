#ifndef __LMTF_INCLUDED
#define __LMTF_INCLUDED

#define INT 0
#define LONG 1
#define DOUBLE 2
#define FLOAT 3
#define PARADIGMID 157 // change this to your paradigm id

// Codes that get dropped in the header (on the first trial only; 8000s range)
#define FIXXCD           8000
#define FIXYCD           8001
#define FPSIZECD         8002
#define EYEWINWCD        8003
#define EYEWINHCD        8004
#define TARGSIZECD       8005
#define TARGWINWCD       8006
#define TARGWINHCD       8007
#define STIMXCD          8008
#define STIMYCD          8009
#define NSTIMCD          8010
#define TRIALSPERBLOCKCD 8011
#define NBLOCKSCD        8012
#define FRAMERATECD      8013
#define PIXPERDEGCD      8014
#define BKGNDRGBCD       8015
#define GAMMATABLECD     8016
#define FUNDSCD          8017
#define MONSPDCD         8018
#define NHYPCD           8019
#define HYPCD            8020
#define PLATFRAMESCD     8021
#define RAMPFRAMESCD     8022
#define THETACD          8023
#define SFCD             8024
#define PHICD            8025
#define GAMMACD          8026
#define SIGMACD          8027
#define NSIGMASCD        8028
#define HEADERREADYCD    8998
#define PARADIGMIDENTCD  8999

// Codes that get dropped on every trial (in the 7000s range)
#define ACTSTIMXCD     7000
#define ACTSTIMYCD     7001
#define LCCCD          7002
#define MCCCD          7003
#define TFCD           7004
#define OOGCD          7005
#define STIMIDXCD      7006
#define CORRECTTRIALCD 7007

typedef struct {
    int flagcodeval;
    int casttype;
    int nelements;
    union {
        long lval;
        float fval;
        double dval;
    } val;
} ecode_struct;

void updateEyewin(int, int, int, int, int);

#define PARADIGM_PRELUDE(name) \
do { static int do__once_ = 1; \
    if ( do__once_ ) { \
        dprintf(name " compiled on %s at %s\n", __DATE__, __TIME__); \
        do__once_ = 0; } \
} while (0)

#define DECL_UDPSEND_FUNC(funcname, reqstr, dest) int funcname() { udp_send_to(reqstr, dest); return 0; }

// the variable names in the following macros must match the declarations in checkMsg()
#define SEND_DVAL_TO_STREAM(CODE) \
do { \
    dsingle = hex2double(message.contents, NULL); \
    sendToPlexDataStream("double", (void*) &dsingle, message.size, CODE, 0); \
} while (0)

#define SEND_DARY_TO_STREAM(CODE) \
do { \
    hex2double(message.contents, (void*) &darray); \
    sendToPlexDataStream("double", (void*) &darray, message.size, CODE, 1); \
} while (0)

#endif /* __LMTF_INCLUDED */
