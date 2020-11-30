#ifndef ISOSAMP_INCLUDED
#define ISOSAMP_INCLUDED

#define INT 0
#define LONG 1
#define DOUBLE 2
#define FLOAT 3
#define PARADIGMID 107 // change this to your paradigm id

// Codes that get dropped in the header (on the first trial only; 8000s range)
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
#define MONKEYCD        8010 // not used after June 2014 but keep it here for reference
#define NUMSTIMSCD      8011
#define BKGNDRGBCD      8012
#define GAMMATABLECD    8013
#define MONSPDCD        8014
#define FRAMERATECD     8015
#define FUNDAMENTALSCD  8016
#define PIXPERDEGCD     8017
#define MODELPARAMSCD   8018
#define PLOTLIMITSCD    8019 // not used after June 2014 but keep it here for reference
#define MODELREVCD      8020
#define MODELSTRCD      8021
#define RSCALESCD       8022
#define THBOUNDSCD      8023
#define ZBOUNDSCD       8024
#define HEADERREADYCD   8998
#define PARADIGMIDENTCD 8999

// Codes that get dropped on every trial (in the 7000s range)
#define GABFLASHTIMECD 7000
#define GABPHICD       7001
#define GABTHETACD     7002
#define GABSFCD        7003
#define GABGAMMACD     7004
#define GABTFCD        7005
#define GABNSIGMASCD   7006
#define GABSIGMACD     7007
#define GABLCD         7008
#define GABMCD         7009
#define GABSCD         7010
#define CORRECTTRIALCD 7011
#define STIMXCD        7012
#define STIMYCD        7013

#define REQSTIMUPDATECD 1040

// Parameters that the user can change
typedef struct {
    int fp_x, fp_y;
    int rf_x, rf_y;
    int behavior, targ_size;
    int eyewin_w, eyewin_h;
    int targwin_w, targwin_h;
    int fp_size, fp_rgb[3];
    int target_nstims;
    int nblocks_max;
    float r[2];
    float th[2];
    float z[2];
    int special_case, blanks, prev_stims;
} params_public;

// Housekeeping (private) globals
typedef struct {
    int stim_x, stim_y;
    int nblocks_done, first_trial;
    int trials_remaining_block, nstims_block;
    int stimon, fpon, done;
    int mac_ready, plx_ready;
    int correct_choice;
} params_private;

// Housekeeping globals for the gabor
typedef struct {
    int flash_time;
    float phi, gamma, tf;
    float num_sigmas, sigma;
    float theta, sf;
    double lms[3];
} params_gabor;

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

#define PARADIGM_PRELUDE( name ) \
do { \
    static int _____do_once = 1; \
    if ( _____do_once ) { \
        dprintf(name " compiled on %s at %s\n", __DATE__, __TIME__); \
        _____do_once = 0; \
    } \
} \
while ( 0 )

#endif /* ISOSAMP_INCLUDED */
