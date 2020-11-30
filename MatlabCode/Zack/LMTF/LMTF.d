/*
 * LMTF: adaptively sample an observer's detection surface in L-/M-cone contrast and temporal frequency space.
 */

//#define TESTING

#include <string.h>
#include <math.h>
#include "ldev.h"
#include "labcodes.h"
#include "UDPtransfer.h"
#include "rigconsts.h"
#include "LMTF.h"

// "public" globals
int gl_fp_x = 0;
int gl_fp_y = 0;
int gl_stim_x = -50;
int gl_stim_y = 0;
int gl_eyewin_w = 5;
int gl_eyewin_h = 5;
int gl_fp_size = 2;
int gl_targ_size = 2;
int gl_targwin_w = 13;
int gl_targwin_h = 13;
int gl_min_tf = 1;
int gl_max_tf = 25;
int gl_human_subject = 0;
int gl_trials_per_block = 8;
int gl_nblocks = 5;

// "private" globals
double gl_stimulus[5]; // Lcc, Mcc, and TF, stim idx, and hemifield
int gl_done = 0;
int gl_trial_count = 0;
int gl_first_trial = 1;
int gl_stimon = 0;
int gl_mac_ready = 0;
int gl_plx_ready = 0;
int gl_stim_location = 0; // 1 -> actual stim location the same as gl_stim_x,y
int gl_correct_choice = 0;
int gl_nstim = 0;
int gl_is_human = 0;
int gl_trials_per_stim = 0;
int gl_trials_left = 0;
int gl_condidx_updated = 0;
//int gl_platframes = 24;
//int gl_rampframes = 13;
int gl_flash_time = 666;
double gl_theta = 0.0f;
double gl_sf = 1.0;
double gl_phi = 0.0;
double gl_gamma = 1.0;
double gl_sigma = 0.15;
double gl_nsigmas = 2.0;
double gl_framerate = 0.0;

#define WINDFIXN 0
#define WINDCORR 1
#define WINDWRNG 2

/*
 * Stuff for TrialParams array
 */

// Indices for gl_trialparams
#define TP_ACTSTIMX  0
#define TP_ACTSTIMY  1
#define TP_LCC       2
#define TP_MCC       3
#define TP_TF        4
#define TP_OOG       5
#define TP_STIMIDX   6
#define TP_CORRECTCD 7
#define NTRIALPARAMS 8

int gl_trialparamcodes[NTRIALPARAMS] = {ACTSTIMXCD, ACTSTIMYCD, LCCCD, MCCCD, TFCD, OOGCD, STIMIDXCD, CORRECTTRIALCD};
int gl_trialparamtypes[NTRIALPARAMS] = {INT, INT, DOUBLE, DOUBLE, DOUBLE, INT, INT, INT};
ecode_struct gl_trialparams[NTRIALPARAMS];

/*
 * Functions for opening/closing UDP channels
 */

int udpOpen()
{
    rigconstants(whichrig());
    udp_open_2con(REXIP, MACIP, PLEXIP, UDPPORT);
    return 0;
}

int udpClose()
{
    udp_close();
    return 0;
}

/*
 * Functions for sending requests via UDP
 */

const char FRAMERATE_REQ[] = "sendToRex(udpCom, gl.framerate, 'double', message);";
const char PIXPERDEG_REQ[] = "sendToRex(udpCom, gl.pixperdeg, 'double', message);";
const char BKGNDRGB_REQ[] = "sendToRex(udpCom, gl.bkgndrgb, 'double', message);";
const char GAMMATABLE_REQ[] = "sendToRex(udpCom, gl.cal.gammaTable, 'double', message);";
const char FUNDAMENTALS_REQ[] = "sendToRex(udpCom, gl.cal.fundamentals, 'double', message);";
const char MONSPD_REQ[] = "sendToRex(udpCom, gl.cal.monSpd, 'double', message);";
const char INITLMTF_REQ[] = "sendToRex(udpCom, init_LMTF(%d,%d,%d,%d), 'double', message);";
const char GETSTIM_REQ[] = "sendToRex(udpCom, next_stimulus(), 'double', message);";

DECL_UDPSEND_FUNC(startSlavePrgm, "LMTFSlave();", "MAC");
DECL_UDPSEND_FUNC(stopSlavePrgm, "return;", "MAC");
DECL_UDPSEND_FUNC(startOnlinePrgm, "LMTFOnline();", "PLEXON");
DECL_UDPSEND_FUNC(stopOnlinePrgm, "return;", "PLEXON");
DECL_UDPSEND_FUNC(getFrameRate, FRAMERATE_REQ, "MAC");
DECL_UDPSEND_FUNC(getPixPerDeg, PIXPERDEG_REQ, "MAC");
DECL_UDPSEND_FUNC(getBkgndrgb, BKGNDRGB_REQ, "MAC");
DECL_UDPSEND_FUNC(getGammaTable, GAMMATABLE_REQ, "MAC");
DECL_UDPSEND_FUNC(getFundamentals, FUNDAMENTALS_REQ, "MAC");
DECL_UDPSEND_FUNC(getMonspd, MONSPD_REQ, "MAC");

/*
 * Mac commands via UDP
 */

int showFP()
{
    char buf[256];
    updateEyewin(WINDFIXN, gl_fp_x, gl_fp_y, gl_eyewin_w, gl_eyewin_h);
    sprintf(buf, "ShowFP(%d, %d, %d)", gl_fp_x, gl_fp_y, gl_fp_size);
    udp_send_to(buf, "MAC");
    return 0;
}

int hideFP()
{
    udp_send_to("HideFP()", "MAC");
    return 0;
}

const char FUNDS_FILENAME[] = "'T_cones_smj10.mat'";
int setupMac()
{
    char buf[512];
    gl_mac_ready = 0;
    sprintf(buf, "InitDisplay(%d, %d, %s, %s)",
        GL_MONDIST, GL_SCREENWIDTH, GL_CALFILENAME, FUNDS_FILENAME);
    udp_send_to(buf, "MAC");
    return 0;
}

/*
 * Getting UDP commands
 */

int checkMsg()
{
    char inMsg[REXBUFF];
    static message_struct message;
    int i, itmp[20], messageAvailable, n_hyp_types;
    double dsingle, darray[768]; // SEND_*_TO_STREAM macros in the header depend on these names

    inMsg[0] = '\0';
    messageAvailable = udp_check(0);

    if ( messageAvailable ) {
        udp_read(inMsg, REXBUFF);
        prepareMsg(&message, inMsg);

        if ( !strcmp(message.header, "MACSTIMON") ) {
            ec_send_code_hi(STIMONCD);
            gl_stimon = 1;
        } else if ( !strcmp(message.header, "MACSTIMOFF") ) {
            ec_send_code_hi(STIMOFFCD);
            gl_stimon = 0;
        } else if ( !strcmp(message.header, "GABORREADY") ) {
            gl_mac_ready = 1;
        } else if ( !strcmp(message.header, GETSTIM_REQ) ) {
            hex2double(message.contents, (void*) &darray);
            memcpy(gl_stimulus, darray, message.size*sizeof(double));
            gl_plx_ready = 1;
        } else if ( !strcmp(message.header, "OOG") ) {
            gl_trialparams[TP_OOG].val.lval = 1;
        } else if ( !strcmp(message.header, "DISPLAYINIT") ) {
            gl_mac_ready = 1;
        } else if ( !strcmp(message.header, "ONLINEINIT") ) {
            gl_plx_ready = 1;
        } else if ( !strcmp(message.header, "PARSEDHEADER") ) {
            gl_plx_ready = 1;
        } else if ( !strcmp(message.header, "CONDIDXUPDATED") ) {
            gl_condidx_updated = 1;
        } else if ( strstr(message.header, "init_LMTF") ) {
            /* darray = [
            # stimuli,
            # of hyperparameter categories (n_hyp_types below),
            # of hyperparameters per category,
            hyperparameter values (categories were lexicographically ordered)
            ]*/
            hex2double(message.contents, (void*) &darray);
            gl_nstim = (int) darray[0];
            n_hyp_types = (int) darray[1];
            for (i = 0; i < n_hyp_types; ++i) itmp[i] = (int) darray[i+2];
            // send # stimuli as an int
            sendToPlexDataStream("int", (void*) &gl_nstim, 1, NSTIMCD, 0);
            // send # of hyperparameters per type (cov, mean, lik, ...others?)
            sendToPlexDataStream("int", (void*) &itmp, n_hyp_types, NHYPCD, 1);
            // send the hyperparameters as a bookended double array
            sendToPlexDataStream("double", (void*) &darray[n_hyp_types+2], message.size-n_hyp_types-4, HYPCD, 1);
            gl_nstim = (int) darray[0];
            gl_trials_left = gl_nstim*gl_trials_per_stim;
            if ((int)darray[message.size-2] != 0 || (int)darray[message.size-1] != 0) {
	            gl_stim_x = (int) darray[message.size-2]; /* last two numbers are RFx and RFy */
   		        gl_stim_y = (int) darray[message.size-1];
            }
            gl_plx_ready = 1;
        } else if ( !strcmp(message.header, FRAMERATE_REQ) ) {
            SEND_DVAL_TO_STREAM(FRAMERATECD);
            gl_framerate = hex2double(message.contents, NULL); 
        } else if ( !strcmp(message.header, PIXPERDEG_REQ) ) {
            SEND_DVAL_TO_STREAM(PIXPERDEGCD);
        } else if ( !strcmp(message.header, BKGNDRGB_REQ) ) {
            SEND_DARY_TO_STREAM(BKGNDRGBCD);
        } else if ( !strcmp(message.header, GAMMATABLE_REQ) ) {
            SEND_DARY_TO_STREAM(GAMMATABLECD);
        } else if ( !strcmp(message.header, FUNDAMENTALS_REQ) ) {
            SEND_DARY_TO_STREAM(FUNDSCD);
        } else if ( !strcmp(message.header, MONSPD_REQ) ) {
            SEND_DARY_TO_STREAM(MONSPDCD);
            ec_send_code_tagged(HEADERREADYCD, VALOFFSET+1);
            // HEADERREADYCD tells the online program that it's safe to
            // parse the header (i.e., REX has sent all relevant variables)
        }
    }
    return 0;
}

/*
 * Functions for dealing with the header and trial parameters
 */

int dropHeader()
{
    int i, nints = 2, ndoubles = 6;
    int nplatframes, nrampframes;
    int intvals[2];
    float nframestot; /* number of total frames can be fractional to avoid roundoff error */
    double dvals[] = {gl_theta, gl_sf, gl_phi, gl_gamma, gl_sigma, gl_nsigmas};
    double dcodes[] = {THETACD, SFCD, PHICD, GAMMACD, SIGMACD, NSIGMASCD};
    
/* For backward compatibility. Have to drop platframes and rampframes in header. */
    int intcodes[] = {PLATFRAMESCD, RAMPFRAMESCD};
  

    nframestot = (float)(gl_framerate*gl_flash_time/1000); /* Needs to be identical to slavecode */
    nrampframes = ceilf(nframestot/4);
    nplatframes = ceilf(nframestot-2*nrampframes);
   // dprintf("nplatframes %d, nrampframes %d",nplatframes,nrampframes);
    
    intvals[0] = nplatframes;
    intvals[1] = nrampframes;
    
    ec_send_code_tagged(PARADIGMIDENTCD , VALOFFSET+PARADIGMID);
    ec_send_code_tagged(FIXXCD          , VALOFFSET+gl_fp_x);
    ec_send_code_tagged(FIXYCD          , VALOFFSET+gl_fp_y);
    ec_send_code_tagged(FPSIZECD        , VALOFFSET+gl_fp_size);
    ec_send_code_tagged(EYEWINWCD       , VALOFFSET+gl_eyewin_w);
    ec_send_code_tagged(EYEWINHCD       , VALOFFSET+gl_eyewin_h);
    ec_send_code_tagged(TARGSIZECD      , VALOFFSET+gl_targ_size);
    ec_send_code_tagged(TARGWINWCD      , VALOFFSET+gl_targwin_w);
    ec_send_code_tagged(TARGWINHCD      , VALOFFSET+gl_targwin_h);
    ec_send_code_tagged(STIMXCD         , VALOFFSET+gl_stim_x);
    ec_send_code_tagged(STIMYCD         , VALOFFSET+gl_stim_y);
    ec_send_code_tagged(TRIALSPERBLOCKCD, VALOFFSET+gl_trials_per_block);
    ec_send_code_tagged(NBLOCKSCD       , VALOFFSET+gl_nblocks);

    for (i = 0; i < nints; ++i)
        sendToPlexDataStream("int", &intvals[i], 1, intcodes[i], 0);
    for (i = 0; i < ndoubles; ++i)
        sendToPlexDataStream("double", &dvals[i], 1, dcodes[i], 0);

    gl_first_trial = 0;
    return 0;
}

int initTrialParams()
{
    int i;

    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
        gl_trialparams[i].casttype = gl_trialparamtypes[i];
        gl_trialparams[i].nelements = 1;
        gl_trialparams[i].val.dval = 0;
    }

    gl_trialparams[TP_LCC].val.dval = gl_stimulus[0];
    gl_trialparams[TP_MCC].val.dval = gl_stimulus[1];
    gl_trialparams[TP_TF].val.dval = gl_stimulus[2];
    gl_trialparams[TP_STIMIDX].val.lval = gl_stimulus[3];

    return 0;
}

int dropTrialParams()
{
    int i;
    char typebuffer[8];

    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        if ( gl_trialparams[i].casttype == FLOAT ) {
            sprintf(typebuffer, "float");
        } else if ( gl_trialparams[i].casttype == DOUBLE ) {
            sprintf(typebuffer, "double");
        } else { // INT
            sprintf(typebuffer, "int");
        }
        sendToPlexDataStream(typebuffer, (void *) &gl_trialparams[i].val.lval,
            gl_trialparams[i].nelements, gl_trialparams[i].flagcodeval, 0);
    }
    return 0;
}

/*
 * REX setup stuff
 */

int resetTrial()
{
    udp_send_to("AllOff();", "MAC");
    updateEyewin(WINDFIXN, gl_fp_x, gl_fp_y, gl_eyewin_w, gl_eyewin_h); // fixation
    updateEyewin(WINDCORR, 1000, 1000, gl_targwin_w, gl_targwin_h); // correct targ offscreen
    updateEyewin(WINDWRNG, 1000, 1000, gl_targwin_w, gl_targwin_h); // incorrect targ offscreen
    gl_stimon = 0;
    gl_plx_ready = 0;
    gl_mac_ready = 0;
    gl_stim_location = 0;
    gl_correct_choice = 0;
    gl_condidx_updated = 0;
    dio_off(FTTLBYTE1);
    return 0;
}

void updateEyewin(int window, int loc_x, int loc_y, int size_w, int size_h)
{
    wd_src_pos(window, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_check(window, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);
    wd_pos(window, loc_x, loc_y);
    wd_siz(window, size_w, size_h);
    wd_cntrl(window, WD_ON);
}

// this function runs whenever you press Reset States (or start the paradigm)
void rinitf()
{
    PARADIGM_PRELUDE("LMTF");

    resetTrial();

    gl_done = 0;
    gl_trial_count = 0;
    gl_first_trial = 1;
    gl_is_human = gl_human_subject;
    gl_trials_per_block += gl_trials_per_block%2; // we want equal amounts of left and right hemifield stim presentations
    gl_trials_per_stim = gl_trials_per_block*gl_nblocks;
    gl_nstim = 0;
    gl_trials_left = 0;
    i_b->t_wrate = 1;
}

/*
 * Functions specific to this paradigm
 */

int initLMTF()
{
    char buf[128];
    gl_plx_ready = 0;
    sprintf(buf, INITLMTF_REQ, gl_min_tf, gl_max_tf, gl_trials_per_block, gl_nblocks);
    udp_send_to(buf, "PLEXON");
    return 0;
}

int getStimulus()
{
    gl_plx_ready = 0;
    udp_send_to(GETSTIM_REQ, "PLEXON");
    return 0;
}

int prepareGabor()
{
    char buf[256];
    gl_mac_ready = 0;
    sprintf(buf, "PrepareGabor(%d, %lf, %lf, %f, %f, %f, %f, %f)",
        gl_flash_time, gl_theta, gl_sf, gl_phi, gl_sigma, gl_gamma, /*TF*/gl_stimulus[2], gl_nsigmas);
        
    udp_send_to(buf, "MAC");
    return 0;
}

int showStim()
{
    char buf[256];
    int tmpx, tmpy;
    gl_stimon = 0;
    gl_stim_location = (int) gl_stimulus[4]; // 1 -> same as to gl_stim_x,y, 0 -> flip signs
    gl_trialparams[TP_ACTSTIMX].val.lval = tmpx = gl_stim_location ? gl_stim_x : -gl_stim_x;
    gl_trialparams[TP_ACTSTIMY].val.lval = tmpy = gl_stim_y;
    sprintf(buf, "ShowStim(%d, %d, [%.15f %.15f %.15f])", tmpx, tmpy,
        gl_stimulus[0], gl_stimulus[1], 0);
    udp_send_to(buf, "MAC");

    return 0;
}

int saveTargChoice(int is_correct)
{
    gl_trialparams[TP_CORRECTCD].val.lval = gl_correct_choice = is_correct;
    return 0;
}

int showTargs()
{
    char buf[64];
    double scale;
    int correct_targx, correct_targy, desired_distance;
    correct_targx = gl_stim_location ? gl_stim_x : -gl_stim_x;
    correct_targy = gl_stim_y;
    desired_distance = gl_is_human ? 171 : 20; // 20 deg tenths for monkeys, 171 for humans w/ button box
    if (gl_is_human) {
        correct_targx = gl_stim_location ? desired_distance: -desired_distance;
        correct_targy = 0;
    } else {
        scale = sqrt(1.0*desired_distance*desired_distance/(correct_targx*correct_targx + correct_targy*correct_targy));
        if (correct_targx > 0) {
	        correct_targx = max(my_round(scale*correct_targx),20); // Target cannot be < 2.0 deg from vertical meridian
        } else {
	        correct_targx = min(my_round(scale*correct_targx),-20); // Target cannot be < 2.0 deg from vertical meridian    
        }
        correct_targy = my_round(scale*correct_targy);
    }
    updateEyewin(WINDCORR, correct_targx, correct_targy, gl_targwin_w, gl_targwin_h);
    updateEyewin(WINDWRNG,-correct_targx, correct_targy, gl_targwin_w, gl_targwin_h);
    sprintf(buf, "ShowTargs(%d, %d, %d)", correct_targx, correct_targy, gl_targ_size);
    udp_send_to(buf, "MAC");
    return 0;
}

int update_trial_counter()
{
    gl_trials_left--;
    if (++gl_trial_count >= gl_trials_per_stim*gl_nstim)
        gl_done = 1;
    return 0;
}

// Incoming codes are queued at index lx (_l_oad) and are dequeued starting at index dx (_d_ump).
// So lx is strictly greater than dx when there's codes waiting to be dumped, and lx == dx
// when the queue is empty.
int lopri_clear()
{
    static int prev_load_idx = 0;
    if (prev_load_idx != i_b->pl_lolx && i_b->pl_lodx == i_b->pl_lolx) {
        prev_load_idx = i_b->pl_lolx;
        return 1;
    }
    return 0;
}

/*
 * Functions that process menu inputs
 */

int sanityChecks(int flag, MENU *mp, char *astr, VLIST *vlp, int *tvadd)
{
    int tmp;
    if (gl_human_subject < 0) gl_human_subject = 0;
    else if (gl_human_subject > 1) gl_human_subject = 1;

    if (gl_fp_size < 0) gl_fp_size = 2;
    if (gl_targ_size < 0) gl_targ_size = 2;
    if (gl_trials_per_block < 0) gl_trials_per_block = 0;
    if (gl_nblocks < 0) gl_nblocks = 0;
    if (gl_min_tf < 0) gl_min_tf = 0;
    if (gl_max_tf < 0) gl_max_tf = 0;
    if (gl_max_tf < gl_min_tf) {
        tmp = gl_min_tf;
        gl_min_tf = gl_max_tf;
        gl_max_tf = tmp;
    }

    return 0;
}

/*
 * Menu declarations
 */

VLIST menu_fixed[] =
{
    {"FP X", &gl_fp_x, NP, NP, 0, ME_DEC},
    {"FP Y", &gl_fp_y, NP, NP, 0, ME_DEC},
    {"FP size", &gl_fp_size, NP, sanityChecks, ME_AFT, ME_DEC},
    {"FP window width", &gl_eyewin_w, NP, NP, 0, ME_DEC},
    {"FP window height", &gl_eyewin_h, NP, NP, 0, ME_DEC},
    {"Targ size", &gl_targ_size, NP, sanityChecks, ME_AFT, ME_DEC},
    {"Targ window width", &gl_targwin_w, NP, NP, 0, ME_DEC},
    {"Targ window height", &gl_targwin_h, NP, NP, 0, ME_DEC},
    {"Minimum TF", &gl_min_tf, NP, sanityChecks, ME_AFT, ME_DEC},
    {"Maximum TF", &gl_max_tf, NP, sanityChecks, ME_AFT, ME_DEC},
    {"Human subject?", &gl_human_subject, NP, sanityChecks, ME_AFT, ME_DEC},
    {"Trials/block", &gl_trials_per_block, NP, sanityChecks, ME_AFT, ME_DEC},
    {"# blocks", &gl_nblocks, NP, sanityChecks, ME_AFT, ME_DEC},
    {NS},
};

VLIST menu_stim[] =
{
    {"Stim X", &gl_stim_x, NP, NP, 0, ME_DEC},
    {"Stim Y", &gl_stim_y, NP, NP, 0, ME_DEC},
    {"Stim dur (ms)", &gl_flash_time, NP, NP, 0, ME_DEC},
        {NS},
};

MENU umenus[] = {
    {"Fixed Params", &menu_fixed, NP, NP, 0, NP, NP},
    {"Stimulus Params", &menu_stim, NP, NP, 0, NP, NP},
    {NS},
};

VLIST state_vl[] = {{NS}};

char hm_sv_vl[] = "";

// User-supplied real-time variable table

RTVAR rtvars[] = {
#ifdef TESTING
    {"mac ready?", &gl_mac_ready},
    {"plex ready?", &gl_plx_ready},
    {"act stim loc same?", &gl_stim_location},
    {"correct choice?", &gl_correct_choice},
    {"max # trials/stim", &gl_trials_per_stim},
    {"# stim", &gl_nstim},
    {"trials completed", &gl_trial_count},
    {"done?", &gl_done},
#endif
    {"Trials remaining", &gl_trials_left},
    {"", 0},
};

/*
 * The state sets
 */

%%
id PARADIGMID
restart rinitf
main_set {
status ON
begin
first:
    to stopPlex
stopPlex:
    do dio_off(PLXREC)
    to closePort
closePort:
    do udpClose()
    to openPort
openPort:
    do udpOpen()
    to stopSlave
stopSlave:
    do stopSlavePrgm()
    time 5
    to stopOnline
stopOnline:
    do stopOnlinePrgm()
    time 5
    to startSlave
startSlave:
    do startSlavePrgm()
    time 200
    to initMac
initMac:
    do setupMac()
    time 5000
    to startOnline on 1 = gl_mac_ready
    to failSafely0
startOnline:
    do startOnlinePrgm()
    to waitForUser on 1 = gl_plx_ready /* wait for user to hit cancel or select text file */
waitForUser:
    to plxStart on -PSTOP & softswitch
getBkgndrgb:
    do getBkgndrgb()
    time 30
    to getFundamentals
getFundamentals:
    do getFundamentals()
    time 1000
    to getMonspd
getMonspd:
    do getMonspd()
    time 2000
    to initLMTF
initLMTF: /* get the algorithm started while we get the rest of the header */
    do initLMTF()
    time 200
    to getGammaTable
getGammaTable:
    do getGammaTable()
    time 3000
    to getFrameRate
getFrameRate:
    do getFrameRate()
    time 30
    to getPixPerDeg
getPixPerDeg:
    do getPixPerDeg()
    time 20
    to dropHeader on 1 = gl_plx_ready
 dropHeader:
    do dropHeader()
    time 100
    to pause0
pause0:
    to pauseLoop on 1 = gl_done
    to pause1 on +PSTOP & softswitch
    to getStimulus
getStimulus:
    do getStimulus()
    time 5000
    to prepareGabor on 1 = gl_plx_ready
    to pauseLoop
prepareGabor:
    do prepareGabor()
    time 5000
    to initTrial on 1 = gl_mac_ready
    to pauseLoop
pauseLoop:
    time 1000
    to pause0
pause1:
    to getStimulus on -PSTOP & softswitch
initTrial:
    do initTrialParams()
    time 200
    to pauseLoop on 1 = gl_done
    to plxStart
plxStart:
    do dio_on(PLXREC)
    time 400
    to getBkgndrgb on 1 = gl_first_trial
    to fpOn
fpOn:
    code FPONCD /* 1004 */
    do showFP()
    time 2500
    to preStimDel on -WD0_XY & eyeflag
    to bail
preStimDel:
    code FPACQCD /* 1032 */
    time 500 /* 300 + 200 (preStimDel & frameOn, resp.) to match DTspot */
    to bail on +WD0_XY & eyeflag
    to showStim
showStim:
    do showStim()
    time 5000
    to enforceFix on 1 = gl_stimon /* STIMONCD dropped */
    to bail
enforceFix:
    time 5000
    to bail on +WD0_XY & eyeflag
    to preTargDel on 0 = gl_stimon /* STIMOFFCD dropped */
    to bail
preTargDel:
    time 100
    rand 500
    to bail on +WD0_XY & eyeflag
    to targsOn
targsOn:
    code TARGONCD /* 1008 */
    do showTargs()
    time 10
    to fpOff
fpOff:
    code FPOFFCD /* 1006 */
    do hideFP()
    to reactTime
reactTime:
    time 700
    to choseCorrect on -WD1_XY & eyeflag
    to choseWrong on -WD2_XY & eyeflag
    to bail
choseCorrect:
    do saveTargChoice(1)
    to targWait
choseWrong:
    do saveTargChoice(0)
    to targWait
targWait:
    code SACMADCD /* 1026 */
    time 200
    to evalChoice
evalChoice:
    to evalCorrect on 1 = gl_correct_choice
    to evalIncorrect
evalCorrect:
    to rewOn on -WD1_XY & eyeflag
    to bail
evalIncorrect:
    to wrongChoice on -WD2_XY & eyeflag
    to bail
rewOn:
    code REWCD /* 1015 */
    do dio_on(REW)
    time 130
    to rewOff
rewOff:
    do dio_off(REW)
    to correctChoice
correctChoice:
    code CORRECTCD /* 1018 */
    do resetTrial()
    to incTrialCounts
wrongChoice:
    code ERRCD /* 1014 */
    do resetTrial()
    to incTrialCounts
incTrialCounts:
    do update_trial_counter()
    to dropTrialParams
bail:
    code ABORTCD /* 1013 */
    do resetTrial()
    to dropTrialParams
dropTrialParams:
    do dropTrialParams() /* this function modifies this state's wait time */
    to markTrialEnd on 1 % lopri_clear
markTrialEnd:
    code EOT /* 1016 */
    to plxStop on 1 = gl_condidx_updated /* need to handshake here!!! */
plxStop:
    do dio_off(PLXREC)
    time 50
    to pause0
failSafely0:
    time 1000
    to failSafely1
failSafely1:
    time 1000
    to failSafely0
abort list:
}

msg_set {
status ON
begin
msgZero:
   to msgFirst
msgFirst:
   do checkMsg()
   to msgFirst
abort list:
}
