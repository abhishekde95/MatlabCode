#include <stdlib.h>
#include "ldev.h"
#include "labcodes.h"
#include "UDPtransfer.h"
#include "rigconsts.h"
#include "ChromCat.h"

//#define TESTING

// "public" globals
float gl_thresh_mult_override = 0.0f;
float gl_tf = 0.0f;
int gl_fp_x = 0;
int gl_fp_y = 0;
int gl_rf_x = 0;
int gl_rf_y = 0;
int gl_eyewin_w  = 7;
int gl_eyewin_h = 7;
int gl_targwin_w = 13;
int gl_targwin_h = 13;
int gl_fp_size = 2;
int gl_targ_size = 2;
int gl_crutchstim_size = 6;
float gl_wrong_targ_contrast = 1.0f;
int gl_stim_dur_public = 500; // must be an even integer
int gl_nblocks = 9999;
int gl_condition_override = 0; // 1 -> always show L+M, 2 -> always show L-M
int gl_color_targets = 0;
int gl_use_gabors = 1; // when 0, show blue or yellow squares instead of Gabors
int gl_subject_id = 0; // 0-> sedna, 1-> greg, otherwise-> zack
int gl_streak_thresh_lo = 999;
int gl_streak_thresh_hi = 999;

// "private" globals
int gl_streak_thresh = 0; // show a "probe" stimulus (i.e., not L+M or L-M) once the current streak reaches this threshold; gated by gl_use_gabors
float gl_cc_at_LMthresh = 0.0f;
float gl_cc_at_LvMthresh = 0.0f;
int gl_LvMtarg_x = 50; // DO NOT change the sign of this--breaks the convention for the right-ward sacc trial parameter
int gl_LvMtarg_y = 0;
int gl_current_streak = 0;
int gl_stimon = 0;
int gl_stim_dur; // the true stimulation duration
int gl_done = 0;
int gl_mac_ready = 0;
int gl_trials_left = 0;
int gl_firsttrial = 1;
int gl_show_LvM = 0; // 0 -> show an L+M stimulus; L-M otherwise
int gl_nconds = 0;
int gl_condidx = 0; // index into gl_condarray
int gl_nblockscompleted = 0; // this counter increments to gl_nblocks
int gl_correct_choice = 0;
int *gl_condorder = NULL;
int **gl_condarray = NULL;

// human readable eye window names
#define WINDFIXN 0
#define WINDCORR 1
#define WINDWRNG 2

/*
 * Stuff for TrialParams array
 */

#define TP_LCC           0
#define TP_MCC           1
#define TP_CORRECT       2
#define TP_TF            3
#define TP_THRESHMULT    4
#define TP_CORRTARGX     5
#define TP_CORRTARGY     6
#define TP_RIGHTWARDSACC 7
#define NTRIALPARAMS     8

int gl_trialparamcodes[NTRIALPARAMS] = {LCCCD, MCCCD, CORRECTTRIALCD, TFCD, THRESHMULTCD, CORRTARGXCD, CORRTARGYCD, RIGHTWARDSACCCD};
int gl_trialparamtypes[NTRIALPARAMS] = {FLOAT, FLOAT, INT, FLOAT, FLOAT, INT, INT, INT};
ecode_struct gl_trialparams[NTRIALPARAMS];

/*
 * Functions for opening/closing UDP channels
 */

int create_udp()
{
    rigconstants(whichrig());
    udp_open_2con(REXIP, MACIP, PLEXIP, UDPPORT);
    return 0;
}

int destroy_udp()
{
    udp_close();
    return 0;
}

/*
 * Functions for sending requests via UDP
 */

const char RFXY_REQ[] = "sendToRex(udpCom, gl.bar.xy, 'double', message);";
const char M_REQ[] = "sendToRex(udpCom, gl.cal.M, 'double', message);";
const char BKGNDRGB_REQ[] = "sendToRex(udpCom, gl.bkgndrgb, 'double', message);";

#define DECL_UDPSEND_FUNC( funcname, reqstr, dest ) int funcname##() { udp_send_to(reqstr, dest); return 0; }
DECL_UDPSEND_FUNC( start_slave, "ThreePulseMonteSlave();", "MAC" )
DECL_UDPSEND_FUNC( stop_slave, "return;", "MAC" )
DECL_UDPSEND_FUNC( get_M, M_REQ, "MAC" )
DECL_UDPSEND_FUNC( get_bkgndrgb, BKGNDRGB_REQ, "MAC" )
DECL_UDPSEND_FUNC( hide_FP, "HideFP();", "MAC" )
DECL_UDPSEND_FUNC( hide_visual, "HideVis();", "MAC" )
#undef DECL_UDPSEND_FUNC

/*
 * Mac commands via UDP
 */

 int getRFxy()
{
    if ( !gl_rf_x && !gl_rf_y )
        udp_send_to(RFXY_REQ, "MAC");
    return 0;
}

int show_FP()
{
    char buf[256];
    update_eyewin(WINDFIXN, gl_fp_x, gl_fp_y, gl_eyewin_w, gl_eyewin_h);
    sprintf(buf, "ShowFP(%d, %d, %d)", gl_fp_x, gl_fp_y, gl_fp_size);
    udp_send_to(buf, "MAC");
    return 0;
}

int prepare_gabor()
{
    char buf[256];
    gl_mac_ready = 0;
    gl_trialparams[TP_TF].val.fval = gl_tf;
    sprintf(buf, "PrepareGabor(%d, %.15f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f)",
        gl_stim_dur, gl_tf, 0.0/*theta*/, 1.0/*sf*/, 0.0/*phi*/, 1.0/*gamma*/,
        0.15/*sigma*/, 2.0/*nsigmas*/);
    udp_send_to(buf, "MAC");
    return 0;
}

int show_visual()
{
    char buf[256];

    gl_stimon = 0;
    if (gl_use_gabors) {
        sprintf(buf, "ShowStim(%d, %d, [%.15f %.15f %.15f])", gl_rf_x, gl_rf_y, gl_trialparams[TP_LCC].val.fval, gl_trialparams[TP_MCC].val.fval, 0);
    } else if (gl_show_LvM) {
        sprintf(buf, "ShowVis(%d, %d, %d, [%d %d %d])", gl_rf_x, gl_rf_y, gl_crutchstim_size, 0, 0, 255);
    } else {
        sprintf(buf, "ShowVis(%d, %d, %d, [%d %d %d])", gl_rf_x, gl_rf_y, gl_crutchstim_size, 255, 255, 0);
    }
    udp_send_to(buf, "MAC");

    return 0;
}

int show_targs()
{
    char buf[64];
    int correct_targx, correct_targy;

    correct_targx = gl_show_LvM ? gl_LvMtarg_x : -gl_LvMtarg_x;
    correct_targy = gl_show_LvM ? gl_LvMtarg_y : -gl_LvMtarg_y;

    gl_trialparams[TP_CORRTARGX].val.lval = correct_targx;
    gl_trialparams[TP_CORRTARGY].val.lval = correct_targy;

    update_eyewin(WINDCORR, correct_targx, correct_targy, gl_targwin_w, gl_targwin_h);
    update_eyewin(WINDWRNG,-correct_targx,-correct_targy, gl_targwin_w, gl_targwin_h);
    sprintf(buf, "ShowTargs(%d, %d, %d, %.4f, %d)", correct_targx, correct_targy,
        gl_targ_size, gl_wrong_targ_contrast, gl_color_targets);
    udp_send_to(buf, "MAC");
    return 0;
}

const char FUNDS_FILENAME[] = "'T_cones_smj10.mat'";
int init_screen()
{
    char buf[256];
    gl_mac_ready = 0;
    sprintf(buf, "InitDisplay(%d, %d, %s, %s)",
        GL_MONDIST, GL_SCREENWIDTH, GL_CALFILENAME, FUNDS_FILENAME);
    udp_send_to(buf, "MAC");
    return 0;
}

/*
 * Getting UDP commands
 */

int check_msg()
{
    char inMsg[REXBUFF];
    static message_struct message;
    int messageAvailable;
    double dvals[16];

    inMsg[0] = '\0';
    messageAvailable = udp_check(0);

    if ( messageAvailable ) {
        udp_read(inMsg, REXBUFF);
        prepareMsg(&message, inMsg);

        if ( !strcmp(message.header, "DISPLAYINIT") ) {
            gl_mac_ready = 1;
        } else if ( !strcmp(message.header, "MACSTIMON") ) {
            ec_send_code_hi(STIMONCD);
            gl_stimon = 1;
        } else if ( !strcmp(message.header, "MACSTIMOFF") ) {
            ec_send_code_hi(STIMOFFCD);
            gl_stimon = 0;
        } else if ( !strcmp(message.header, "GABORREADY") ) {
            gl_mac_ready = 1;
        } else if ( !strcmp(message.header, RFXY_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            if ( message.size && (dvals[0] || dvals[1]) ) {
                gl_rf_x = my_round(10. * dvals[0]);
                gl_rf_y = my_round(10. * dvals[1]);
            }
        } else if ( !strcmp(message.header, M_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            sendToPlexDataStream("double", (void *) &dvals, message.size, MMATRIXCD, 1);
        } else if ( !strcmp(message.header, BKGNDRGB_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            sendToPlexDataStream("double", (void *) &dvals, message.size, BKGNDRGBCD, 1);
        }
    }

    return 0;
}

/*
 * Functions for dealing with the header and trial parameters
 */

int drop_header()
{
    ec_send_code_tagged(PARADIGMIDENTCD , VALOFFSET+PARADIGMID);
    ec_send_code_tagged(FIXXCD          , VALOFFSET+gl_fp_x);
    ec_send_code_tagged(FIXYCD          , VALOFFSET+gl_fp_y);
    ec_send_code_tagged(FPSIZECD        , VALOFFSET+gl_fp_size);
    ec_send_code_tagged(EYEWINWCD       , VALOFFSET+gl_eyewin_w);
    ec_send_code_tagged(EYEWINHCD       , VALOFFSET+gl_eyewin_h);
    ec_send_code_tagged(TARGXCD         , VALOFFSET+gl_LvMtarg_x);
    ec_send_code_tagged(TARGYCD         , VALOFFSET+gl_LvMtarg_y);
    ec_send_code_tagged(TARGSIZECD      , VALOFFSET+gl_targ_size);
    ec_send_code_tagged(TARGWINWCD      , VALOFFSET+gl_targwin_w);
    ec_send_code_tagged(TARGWINHCD      , VALOFFSET+gl_targwin_h);
    ec_send_code_tagged(STIMDURCD       , VALOFFSET+gl_stim_dur);
    ec_send_code_tagged(NBLOCKSCD       , VALOFFSET+gl_nblocks);

    gl_firsttrial = 0;
    return 0;
}

int init_trial_params()
{
    int i;
    float *this_condition, thresh_mults[] = {1.5, 1.75, 2};
    int n_thresh_mults = 3;

    // initialize all trial parameters to zero
    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
        gl_trialparams[i].casttype = gl_trialparamtypes[i];
        gl_trialparams[i].nelements = 1;
        gl_trialparams[i].val.dval = 0;
    }

    if (gl_current_streak == 0) // roll a new streak threshold
        gl_streak_thresh = gl_streak_thresh_lo + my_randint(gl_streak_thresh_hi-gl_streak_thresh_lo+1);

    this_condition = gl_condarray[gl_condorder[gl_condidx]];

    if (gl_use_gabors) {
        if (gl_thresh_mult_override > 0)
            gl_trialparams[TP_THRESHMULT].val.fval = gl_thresh_mult_override;
        else
            gl_trialparams[TP_THRESHMULT].val.fval = thresh_mults[my_randint(n_thresh_mults)];
    }

    gl_show_LvM = gl_condition_override ? (gl_condition_override - 1) : my_randint(2);

    if (gl_use_gabors) {
        if (!gl_condition_override && gl_current_streak >= gl_streak_thresh) { // show a probe stimulus
            printf("probe (#%d): %f %f\n", gl_condorder[gl_condidx], this_condition[TP_LCC], this_condition[TP_MCC]);
            gl_trialparams[TP_LCC].val.fval = this_condition[TP_LCC];
            gl_trialparams[TP_MCC].val.fval = this_condition[TP_MCC];
        } else if (gl_show_LvM) { // L-M
            gl_trialparams[TP_LCC].val.fval = gl_trialparams[TP_THRESHMULT].val.fval*gl_cc_at_LvMthresh;
            gl_trialparams[TP_MCC].val.fval = -gl_trialparams[TP_THRESHMULT].val.fval*gl_cc_at_LvMthresh;
        } else { // L+M
            gl_trialparams[TP_LCC].val.fval = gl_trialparams[TP_THRESHMULT].val.fval*gl_cc_at_LMthresh;
            gl_trialparams[TP_MCC].val.fval = gl_trialparams[TP_THRESHMULT].val.fval*gl_cc_at_LMthresh;
        }
    }

    set_times("visTime", gl_stim_dur, 0);

    return 0;
}

int drop_trial_params()
{
    int i, ncodes = 0;
    char typebuffer[8];
    float mspercode = 0.5f; // Conservative; actual value is closer to 0.33

    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        if ( gl_trialparams[i].casttype == FLOAT ) {
            sprintf(typebuffer, "float");
            ncodes += 5;
        } else if ( gl_trialparams[i].casttype == DOUBLE ) {
            sprintf(typebuffer, "double");
            ncodes += 9;
        } else { // INT
            sprintf(typebuffer, "int");
            ncodes += 2;
        }
        sendToPlexDataStream(typebuffer,
                            (void *) &gl_trialparams[i].val.lval,
                            gl_trialparams[i].nelements,
                            gl_trialparams[i].flagcodeval, 0);
    }
    set_times("dropTrialParams", my_ceil(ncodes*mspercode), 0);
    return 0;
}

/*
 * Functions specific to this paradigm
 */

int init_conds()
{
    free_memory();
    build_conds();
    init_conds_order();
    shuffle_conds();
    gl_trials_left = gl_nconds*gl_nblocks;
    return 0;
}

void build_conds()
{
    int i, nsamples = 10;
    float *tmpptr;
    float (*subject_data)[2];

    // These are 10 equally spaced color directions at 1.75x detection
    // threshold on an ellipse whose axes are the L+M and L-M detection thresholds.
    // The number at the end of the variable names denotes the temporal frequency.

    float sedna6[][2] = { // 1.72x detection threshold NOT 1.75 like the others!
        {0.109610288,-0.082053234},
        {0.119159436,-0.065065975},
        {0.125740173,-0.046898660},
        {0.129362469,-0.028141083},
        {0.130163091,-0.009309448},
        {0.128354238,0.009180076},
        {0.124176343,0.027012910},
        {0.117863052,0.043960646},
        {0.109620330,0.059857229},
        {0.099618040,0.074573131},
    };

    float sedna5[][2] = {
        {0.098543013,-0.073768376},
        {0.108605040,-0.059302838},
        {0.116801153,-0.043564578},
        {0.122921108,-0.026739850},
        {0.126743020,-0.009064840},
        {0.128049825,0.009158304},
        {0.126655190,0.027552150},
        {0.122436612,0.045666495},
        {0.115370037,0.062996806},
        {0.105557321,0.079019221},
    };

    float sedna3[][2] = {
        {0.074970210,-0.056121996},
        {0.084505124,-0.046143289},
        {0.094037642,-0.035074227},
        {0.103622133,-0.022541615},
        {0.113185210,-0.008095166},
        {0.122407934,0.008754788},
        {0.130518375,0.028392535},
        {0.136016438,0.050731508},
        {0.136571212,0.074573524},
        {0.129703170,0.097094577},
    };

    float sedna1[][2] = {
        {0.061203390,-0.045816284},
        {0.069660940,-0.038037751},
        {0.078769685,-0.029379574},
        {0.088935726,-0.019346783},
        {0.100698275,-0.007202083},
        {0.114782100,0.008209378},
        {0.132063989,0.028728762},
        {0.153024126,0.057075048},
        {0.175027865,0.095572445},
        {0.185062192,0.138535822},
    };

    float greg5[][2] = {
        {0.070484801,-0.052764261},
        {0.072331731,-0.039496113},
        {0.072032654,-0.026866791},
        {0.070322231,-0.015297665},
        {0.067716485,-0.004843179},
        {0.064530406,0.004615306},
        {0.060933620,0.013255298},
        {0.056998825,0.021259462},
        {0.052733749,0.028794806},
        {0.048099068,0.036006511},
    };

    float zack5[][2] = {
        {0.039454795,-0.029535490},
        {0.040076214,-0.021883268},
        {0.040467733,-0.015093684},
        {0.040698984,-0.008853522},
        {0.040806655,-0.002918550},
        {0.040806655,0.002918550},
        {0.040698984,0.008853522},
        {0.040467733,0.015093684},
        {0.040076214,0.021883268},
        {0.039454795,0.029535490},
    };

    // These are the L and M components for each subjects' detection threshold
    // in the L+M and L-M directions.
    if (gl_subject_id == -3) {
        gl_cc_at_LMthresh = 0.051157680;
        gl_cc_at_LvMthresh = 0.056541081;
        subject_data = sedna6;
        gl_tf = 6.0f;
    } else if (gl_subject_id == -2) {
        gl_cc_at_LMthresh = 0.093407327;
        gl_cc_at_LvMthresh = 0.030326606;
        subject_data = sedna1;
        gl_tf = 1.0f;
    } else if (gl_subject_id == -1) {
        gl_cc_at_LMthresh = 0.065465702;
        gl_cc_at_LvMthresh = 0.037319144;
        subject_data = sedna3;
        gl_tf = 3.0f;
    } else if (gl_subject_id == 0) {
        gl_cc_at_LMthresh = 0.053278453;
        gl_cc_at_LvMthresh = 0.049603738;
        subject_data = sedna5;
        gl_tf = 5.0f;
    } else if (gl_subject_id == 1) {
        gl_cc_at_LMthresh = 0.024580831;
        gl_cc_at_LvMthresh = 0.037524562;
        subject_data = greg5;
        gl_tf = 5.0f;
    } else {
        gl_cc_at_LMthresh = 0.021985642;
        gl_cc_at_LvMthresh = 0.021985642;
        subject_data = zack5;
        gl_tf = 5.0f;
    }

    gl_condarray = malloc(sizeof(float*) * nsamples);

    for ( i = 0; i < nsamples; ++i ) {
        tmpptr = malloc(2*sizeof(float));
        tmpptr[TP_LCC] = subject_data[i][0];
        tmpptr[TP_MCC] = subject_data[i][1];
        gl_condarray[gl_nconds++] = tmpptr;
    }

    gl_done = gl_nconds == 0;
}

void init_conds_order()
{
    int i;
    gl_condorder = malloc(sizeof(int) * gl_nconds);
    for (i = 0; i < gl_nconds; ++i)
        gl_condorder[i] = i;
}

// generate a random integer on [0,n-1] without modulo bias
static int my_randint(int n)
{
    int rnd, limit = RAND_MAX - (RAND_MAX % n);

    do { rnd = rand(); } while ( rnd >= limit );
    return rnd % n;
}

void shuffle_conds()
{
    int i, r, swap;

    for ( i = gl_nconds-1; i > 0; --i ) {
        r = my_randint(i+1);
        swap = gl_condorder[r];
        gl_condorder[r] = gl_condorder[i];
        gl_condorder[i] = swap;
    }
}

/*
    0    1 = is_correct
0   R    L

1   L    R
^
gl_show_LvM
*/

int reset_if_probe()
{
    if (gl_current_streak >= gl_streak_thresh) gl_current_streak = 0;
    return 0;
}

int evaluate_choice(int is_correct)
{
    gl_trialparams[TP_RIGHTWARDSACC].val.lval = !(is_correct^gl_show_LvM); // RIGHTWARD -> 1, LEFTWARD -> 0
    gl_trialparams[TP_CORRECT].val.lval = gl_correct_choice = is_correct || gl_current_streak >= gl_streak_thresh; // read as: correct choice OR it's a probe trial

    if (!gl_correct_choice) // reset streak when the subject made an incorrect choice (and it's not a probe trial)
        gl_current_streak = 0;
    else
        gl_current_streak++;

    return 0;
}

int reset_trial()
{
    dio_off(FTTLBYTE1);
    udp_send_to("AllOff();", "MAC");
    update_eyewin(WINDFIXN, gl_fp_x, gl_fp_y, gl_eyewin_w, gl_eyewin_h); // fixation
    update_eyewin(WINDCORR, 1000, 1000, gl_targwin_w, gl_targwin_h); // correct targ offscreen
    update_eyewin(WINDWRNG, 1000, 1000, gl_targwin_w, gl_targwin_h); // incorrect targ offscreen

    gl_show_LvM = 0;
    gl_correct_choice = 0;
    gl_mac_ready = 0;
    gl_stimon = 0;

    return 0;
}

int update_trial_counter()
{
    if ( !gl_use_gabors || gl_current_streak > gl_streak_thresh ) {
        gl_trials_left--;
        gl_current_streak = 0;
        if ( ++gl_condidx == gl_nconds ) {
            if ( ++gl_nblockscompleted >= gl_nblocks ) {
                gl_done = 1;
                dprintf("maximum number of trials reached\n");
            } else
                shuffle_conds();
            gl_condidx = 0;
        }
    }
    return 0;
}

int free_memory()
{
    int i;

    if ( gl_condarray ) {
        for ( i = 0; i < gl_nconds; ++i )
            free(gl_condarray[i]);
        free(gl_condarray); gl_condarray = NULL;
    }
    free(gl_condorder); gl_condorder = NULL;
    return 0;
}

/*
 * REX setup stuff
 */

// this function is called whenever you press Reset States
void rinitf()
{
    PARADIGM_PRELUDE("ChromCat");

    reset_trial();
    srand(time(0));

    gl_stim_dur = gl_stim_dur_public;
    gl_firsttrial = 1;
    gl_nconds = 0;
    gl_condidx = 0;
    gl_nblockscompleted = 0;
    gl_done = 0;
    gl_trials_left = 0;
    gl_current_streak = 0;

    i_b->t_wrate = 1;
}

void update_eyewin(int window, int loc_x, int loc_y, int size_w, int size_h)
{
    wd_src_pos(window, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_check(window, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);
    wd_pos(window, loc_x, loc_y);
    wd_siz(window, size_w, size_h);
    wd_cntrl(window, WD_ON);
}

/*
 * Functions that process menu inputs
 */

int sanitize(int flag, MENU *mp, char *astr, VLIST *vlp, int *tvadd)
{
    if (gl_stim_dur_public <= 0) gl_stim_dur_public = 300;
    gl_stim_dur_public += gl_stim_dur_public%2;

    if (gl_streak_thresh_lo < 0) gl_streak_thresh_lo = 0;
    if (gl_streak_thresh_hi < 0) gl_streak_thresh_hi = 0;
    if (gl_streak_thresh_lo > gl_streak_thresh_hi) gl_streak_thresh_hi = gl_streak_thresh_lo;

    if (gl_wrong_targ_contrast < 0) gl_wrong_targ_contrast = 0.0f;
    else if (gl_wrong_targ_contrast > 1) gl_wrong_targ_contrast = 1.0f;

    if (gl_condition_override < 0) gl_condition_override = 0;
    else if (gl_condition_override > 2) gl_condition_override = 2;

    if (gl_nblocks-gl_nblockscompleted <= 1)
        gl_trials_left = gl_nconds-gl_condidx;
    else
        gl_trials_left = gl_nconds*(gl_nblocks-gl_nblockscompleted)-gl_condidx;

    return 0;
}

/*
 * Menu declarations
 */

VLIST menu_fp[] =
{
    {"X", &gl_fp_x, NP, NP, 0, ME_DEC},
    {"Y", &gl_fp_y, NP, NP, 0, ME_DEC},
    {"Size", &gl_fp_size, NP, NP, 0, ME_DEC},
    {"Window width", &gl_eyewin_w, NP, NP, 0, ME_DEC},
    {"Window height", &gl_eyewin_h, NP, NP, 0, ME_DEC},
    {NS},
};

VLIST menu_targ[] =
{
    {"Size", &gl_targ_size, NP, NP, 0, ME_DEC},
    {"Crutch stim size", &gl_crutchstim_size, NP, NP, 0, ME_DEC},
    {"Window width", &gl_targwin_w, NP, NP, 0, ME_DEC},
    {"Window height", &gl_targwin_h, NP, NP, 0, ME_DEC},
    {"Wrong targ contrast [0-1]", &gl_wrong_targ_contrast, NP, sanitize, ME_AFT, ME_FLOAT},
    {NS},
};

VLIST menu_expt[] =
{
    {"RF X", &gl_rf_x, NP, NP, 0, ME_DEC},
    {"RF Y", &gl_rf_y, NP, NP, 0, ME_DEC},
    {"Duration", &gl_stim_dur_public,  NP, sanitize, ME_AFT, ME_DEC},
    {"# blocks", &gl_nblocks, NP, sanitize, ME_AFT, ME_DEC},
    {"Cond override (0-2)", &gl_condition_override, NP, sanitize, ME_AFT, ME_DEC},
    {"Color targets?", &gl_color_targets, NP, NP, 0, ME_DEC},
    {"Streak thresh lo", &gl_streak_thresh_lo, NP, sanitize, ME_AFT, ME_DEC},
    {"Streak thresh hi", &gl_streak_thresh_hi, NP, sanitize, ME_AFT, ME_DEC},
    {"Subject? (S=0,G=1,Z=2)", &gl_subject_id, NP, NP, 0, ME_DEC},
    {NS},
};

VLIST menu_gabor[] =
{
    {"Gabors?", &gl_use_gabors, NP, NP, 0, ME_DEC},
    {"Thresh mult override", &gl_thresh_mult_override, NP, NP, 0, ME_FLOAT},
    {"Temp freq", &gl_tf, NP, NP, 0, ME_FLOAT},
    {NS},
};

MENU umenus[] = {
    {"Fixation", &menu_fp, NP, NP, 0, NP, NP},
    {"Target", &menu_targ, NP, NP, 0, NP, NP},
    {"Experiment", &menu_expt, NP, NP, 0, NP, NP},
    {"Gabor", &menu_gabor, NP, NP, 0, NP, NP},
    {NS},
};

VLIST state_vl[] = {{NS}};

char hm_sv_vl[] = "";

RTVAR rtvars[] = {
#ifdef TESTING
    {"mac ready?", &gl_mac_ready},
    {"# conds", &gl_nconds},
    {"# blocks", &gl_nblocks},
    {"correct choice?", &gl_correct_choice},
    {"done?", &gl_done},
#endif
    {"Current streak", &gl_current_streak},
    {"Streak thresh", &gl_streak_thresh},
    {"# conds shown", &gl_condidx},
    {"# blocks completed", &gl_nblockscompleted},
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
    to closePort
closePort:
    do destroy_udp()
    to openPort
openPort:
    do create_udp()
    to stopSlave
stopSlave:
    do stop_slave()
    time 50
    to startSlave
startSlave:
    do start_slave()
    time 50
    to initMac
initMac:
    do init_screen()
    time 5000
    to getRFxy on 1 = gl_mac_ready
getRFxy:
    do getRFxy()
    time 10
    to initConds
initConds:
    do init_conds()
    time 50
    to waitForUser
waitForUser:
    to finished on 1 = gl_done
    to initTrial on -PSTOP & softswitch
initTrial:
    do init_trial_params()
    time 200
    to prepareGabor
prepareGabor:
    do prepare_gabor()
    time 5000
    to plxStart on 1 = gl_mac_ready
plxStart:
    do dio_on(PLXREC)
    time 400
    to dropHeader on 1 = gl_firsttrial
    to fpOn
dropHeader:
    do drop_header()
    time 300
    to getM
getM:
    do get_M()
    time 50
    to getBkgndrgb
getBkgndrgb:
    do get_bkgndrgb()
    time 20
    to fpOn
fpOn:
    code FPONCD
    do show_FP()
    time 2500
    to fixTime on -WD0_XY & eyeflag
    to bail
fixTime:
    code FPACQCD
    time 100
    to bail on +WD0_XY & eyeflag
    to showVisual
showVisual:
    do show_visual()
    to visTime
visTime: /* stimulus duration = gl_stim_dur */
    to probeCheckBail on +WD0_XY & eyeflag
    to hideVisual
hideVisual:
    do hide_visual()
    to preTargsDel
preTargsDel:
    time 200 rand 500
    to probeCheckBail on +WD0_XY & eyeflag
    to fpOff
fpOff:
    code FPOFFCD
    do hide_FP()
    to probeCheckBail on +WD0_XY & eyeflag
    to targsOn
targsOn:
    code TARGONCD
    do show_targs()
    to reactTime
reactTime:
    time 700
    to choseCorrect on -WD1_XY & eyeflag
    to choseWrong on -WD2_XY & eyeflag
    to probeCheckBail
probeCheckBail: /* reset streak if it was a probe stim */
    do reset_if_probe()
    to bail
choseCorrect:
    do evaluate_choice(1)
    to targWait
choseWrong:
    do evaluate_choice(0)
    to targWait
targWait:
    code SACMADCD
    time 200
    to evalChoice
evalChoice:
    to rewOn on 1 = gl_correct_choice
    to wrongChoice
rewOn:
    code REWCD
    do dio_on(REW)
    time 130
    to rewOff
rewOff:
    do dio_off(REW)
    to correctChoice
correctChoice:
    code CORRECTCD
    do reset_trial()
    to incTrialCounts
wrongChoice:
    code ERRCD
    do reset_trial()
    time 1700
    to incTrialCounts
incTrialCounts:
    do update_trial_counter()
    to dropTrialParams
bail:
    code ABORTCD
    do reset_trial()
    to dropTrialParams
dropTrialParams:
    do drop_trial_params()
    to markTrialEnd
markTrialEnd:
    code EOT
    time 15
    to plxStop
plxStop:
    do dio_off(PLXREC)
    time 350
    to waitForUser
finished:
    time 50
    to waitForUser
abort list:
}

msg_set {
status ON
begin
msgZero:
    to msgFirst
msgFirst:
    do check_msg()
    to msgFirst
abort list:
}
