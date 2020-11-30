#ifndef __GRIDLMSUBUNIT_INCLUDED
#define __GRIDLMSUBUNIT_INCLUDED

#define INT 0
#define LONG 1
#define DOUBLE 2
#define FLOAT 3
#define PARADIGMID 112

// Codes that get dropped in the header (on the first trial only; 8000s range)
#define FIXXCD         			8000
#define FIXYCD          		8001
#define FPSIZECD				8002
#define FPRCD           		8003
#define FPGCD           		8004
#define FPBCD           		8005
#define EYEWINXCD       		8006
#define EYEWINYCD       		8007
#define FRAMERATECD   			8008
#define BKGNDRGBCD    			8009
#define PIXPERDEGCD    			8010
#define GAMMATABLECD 			8011
#define FUNDAMENTALSCD 			8012
#define MONSPDCD 				8013
#define RFXCD					8014
#define RFYCD					8015
#define HEADERCMPLTCD			8998
#define PARADIGMIDENTCD			8999

// Codes that get dropped on every trial (in the 7000s range)
#define LCCCD						7000
#define MCCCD						7001
#define SCCCD						7002
#define GRIDXCD						7003
#define GRIDYCD						7004
#define	NSTIXGRIDCD					7005
#define DVAPERSTIXCD				7006
#define STIMDURCD					7007
#define RFCORRCD					7008
#define CHOSECORR					7009

// Some macros
#define DECL_UDPSEND_FUNC( funcname, reqstr, dest ) int funcname##() { udp_send_to(reqstr, dest); return 0; }

// Parameters that the user can change
typedef struct {
	int fp_x, fp_y;
	float rf_x, rf_y;
	int eyewin_x, eyewin_y;
	int fp_size, fp_rgb[3];
	double stim_lmscc[3];
	int grid_x[600], grid_y[600];
	int nstixgrid;
	float dvaperstix;
	float stimdur;
} params_public;

typedef struct {
	int done;
	int firsttrial;
	int macready;
	int plexready;
	int rexready;
	int paradigmdone;
	int getnewstim;
	int fpon;
	int fpoff;
	int grid_len;
	float bkgndrgb[3];
} params_private;

// Housekeeping globals for the current trial
typedef struct {
	long timer_begin, timer_end, fpoff_t; // timers used in state set
} params_trial;

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

#endif /* __GRIDLMSUBUNIT_INCLUDED */