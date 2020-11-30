/*
 * ThreePulseMonte: how finely can monkeys discriminate the relative timing of three electrical pulses in area V1?
 */

#include <stdlib.h>
#include "ldev.h"
#include "labcodes.h"
#include "UDPtransfer.h"
#include "rigconsts.h"
#include "ThreePulseMonte.h"

//#define TESTING

// "public" globals
float gl_thresh_mult_override = 0.0f;
float gl_tf = 5.0f;
int gl_fp_x = 0;
int gl_fp_y = 0;
int gl_rf_x = 0;
int gl_rf_y = 0;
int gl_eyewin_w  = 7;
int gl_eyewin_h = 7;
int gl_targwin_w = 13;
int gl_targwin_h = 13;
int gl_targ_x = 50; // the x coordinate of the target when the middle pulse is delayed; DON'T CHANGE THIS SIGN otherwise "rightward sacc" convention flips
int gl_targ_y = 0; // ditto for y
int gl_fp_size = 2;
int gl_targ_size = 2;
int gl_crutchstim_size = 6;
float gl_wrong_targ_contrast = 1.0f;
int gl_stim_dur_public = 500; // must be an even integer
int gl_npulse_offsets = 2; // how many interleaved pulse times in (0,dur/2) and (dur/2,dur)? e.g., 1 -> 1 pulse at dur/4 and 1 at 3*dur/4 (optionally dur/2)
int gl_opt_middle_pulse = 0; // whether to include a pulse at dur/2 (reward subject via coin flip)
int gl_nblocks = 9999;
int gl_coinflip_override = 1; // flip a coin every trial to pick a pulse offset from (0,dur/2) when tails (0) or (dur/2,dur) when heads (1)
int gl_condition_override = 0; // 1 -> always pick early 2nd pulse conditions, 2 -> pick delayed 2nd pulse conditions
int gl_color_targets = 0;
int gl_use_gabors = 1;
int gl_show_probe_public = 0;

// "private" globals
int gl_show_probe = 0;
float gl_L = 0.0f;
float gl_M = 0.0f;
int gl_stimon = 0;
int gl_stim_dur; // the true stimulation duration
int gl_done = 0;
int gl_mac_ready = 0;
int gl_trials_left = 0;
int gl_firsttrial = 1;
int gl_delay_pulse = 0; // delay the 2nd stim pulse? 0 (no) -> ||....|, 1 (yes) -> |....||
int gl_pulse_offset = 0;
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

#define TP_PULSE2OFFSET 0 // cond param
#define TP_CORRECTCD    1
#define TP_THRESHMULT   2
#define TP_LCC          3
#define TP_MCC          4
#define TP_CORRTARGX     5
#define TP_CORRTARGY     6
#define TP_RIGHTWARDSACCCD 7
#define NTRIALPARAMS    8
#define NCONDPARAMS     1

int gl_trialparamcodes[NTRIALPARAMS] = {PULSE2OFFSETCD, CORRECTTRIALCD, THRESHMULTCD, LCCCD, MCCCD, CORRTARGXCD, CORRTARGYCD, RIGHTWARDSACCCD};
int gl_trialparamtypes[NTRIALPARAMS] = {INT, INT, FLOAT, FLOAT, FLOAT, INT, INT, INT};
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

#define DECL_UDPSEND_FUNC( funcname, reqstr, dest ) int funcname##() { udp_send_to(reqstr, dest); return 0; }
DECL_UDPSEND_FUNC( start_slave, "ThreePulseMonteSlave();", "MAC" )
DECL_UDPSEND_FUNC( stop_slave, "return;", "MAC" )
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
    sprintf(buf, "PrepareGabor(%d, %.15f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f)",
        gl_stim_dur, gl_tf, 0.0/*theta*/, 1.0/*sf*/, 0.0/*phi*/, 1.0/*gamma*/,
        0.15/*sigma*/, 2.0/*nsigmas*/);
    udp_send_to(buf, "MAC");
    return 0;
}

int show_visual()
{
    char buf[256];
    float thresh_mult = gl_trialparams[TP_THRESHMULT].val.fval;

    gl_stimon = 0;    
    
    if (gl_use_gabors) {
        sprintf(buf, "ShowStim(%d, %d, [%.15f %.15f %.15f])", gl_rf_x, gl_rf_y, gl_trialparams[TP_LCC].val.fval, gl_trialparams[TP_MCC].val.fval, 0);
    } else if (gl_delay_pulse) {
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
    // gl_delay_pulse 1 -> don't flip the signs on the user-specified targ x,y; 0 -> flip 'em
    // i.e., what the user specifies in gl_targ_{x,y} (via menus) IS the location
    // of the correct target for a delayed 2nd pulse
    correct_targx = gl_delay_pulse ? gl_targ_x : -gl_targ_x;
    correct_targy = gl_delay_pulse ? gl_targ_y : -gl_targ_y;
    
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
    double dvals[2];

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
    ec_send_code_tagged(TARGXCD         , VALOFFSET+gl_targ_x);
    ec_send_code_tagged(TARGYCD         , VALOFFSET+gl_targ_y);
    ec_send_code_tagged(TARGSIZECD      , VALOFFSET+gl_targ_size);
    ec_send_code_tagged(TARGWINWCD      , VALOFFSET+gl_targwin_w);
    ec_send_code_tagged(TARGWINHCD      , VALOFFSET+gl_targwin_h);
    ec_send_code_tagged(STIMDURCD       , VALOFFSET+gl_stim_dur);
    ec_send_code_tagged(NPULSESCD       , VALOFFSET+gl_npulse_offsets);
    ec_send_code_tagged(OPTPULSECD      , VALOFFSET+gl_opt_middle_pulse);
    ec_send_code_tagged(NBLOCKSCD       , VALOFFSET+gl_nblocks);

    gl_firsttrial = 0;
    return 0;
}

int init_trial_params()
{
    int i, *this_condition;
    float thresh_mults[] = {1.5, 1.75, 2};
    int n_thresh_mults = 3;

    this_condition = gl_condarray[gl_condorder[gl_condidx]];
    gl_pulse_offset = this_condition[TP_PULSE2OFFSET];

    gl_show_probe = gl_show_probe_public;
    
    if (gl_thresh_mult_override > 0)
        gl_trialparams[TP_THRESHMULT].val.fval = gl_thresh_mult_override;
    else
        gl_trialparams[TP_THRESHMULT].val.fval = thresh_mults[my_randint(n_thresh_mults)];

    // if the 2nd pulse offset is precisely half of the stimulation duration, we flip a coin
    if ( !gl_condition_override && gl_opt_middle_pulse && (gl_pulse_offset == gl_stim_dur/2) ) {
        gl_delay_pulse = my_randint(2); // flip a coin
    } else if ( !gl_coinflip_override && !gl_condition_override ) {
        // set gl_delay_pulse to 1 if the 2nd pulse occurs in the (dur/2,dur) interval
        gl_delay_pulse = gl_pulse_offset > gl_stim_dur/2;
    } else { // coin flip override (i.e., ignore the block format)
        gl_delay_pulse = gl_condition_override ? (gl_condition_override - 1) : my_randint(2);
        if ( gl_delay_pulse ) { // pick a pulse offset from (dur/2,dur)
            do {
                gl_pulse_offset = gl_condarray[my_randint(gl_nconds)][TP_PULSE2OFFSET];
            } while ( gl_pulse_offset < gl_stim_dur/2 );
        } else { // or from (0,dur/2)
            do {
                gl_pulse_offset = gl_condarray[my_randint(gl_nconds)][TP_PULSE2OFFSET];
            } while ( gl_pulse_offset > gl_stim_dur/2 );
        }
    }
    
    // Sedna's threshold at 5 degrees and 5 Hz
    if (gl_show_probe && (gl_L > 0 || gl_M > 0)) {
        gl_trialparams[TP_LCC].val.fval = gl_L;
        gl_trialparams[TP_MCC].val.fval = gl_M;
    } else if (gl_delay_pulse) { // L-M
        gl_trialparams[TP_LCC].val.fval = gl_trialparams[TP_THRESHMULT].val.fval*0.049604;
        gl_trialparams[TP_MCC].val.fval = -gl_trialparams[TP_THRESHMULT].val.fval*0.049604;
    } else { // L+M
        gl_trialparams[TP_LCC].val.fval = gl_trialparams[TP_THRESHMULT].val.fval*0.053278;
        gl_trialparams[TP_MCC].val.fval = gl_trialparams[TP_THRESHMULT].val.fval*0.053278;
    }

    set_times("pulse1off", gl_pulse_offset, 0);
    set_times("pulse2off", gl_stim_dur-gl_pulse_offset, 0);

    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
        gl_trialparams[i].casttype = gl_trialparamtypes[i];
        gl_trialparams[i].nelements = 1;
    }

    // load in any trialparameters derived from the conditions array
    for ( i = 0; i < NCONDPARAMS; ++i ) {
        switch (gl_trialparamtypes[i]) {
            case INT:
                gl_trialparams[i].val.lval = (i == TP_PULSE2OFFSET) ? gl_pulse_offset : this_condition[i];
                break;
        }
    }

    // load in any trialparameters that are independent of the conditions array
    for ( i = NCONDPARAMS; i < NTRIALPARAMS; ++i ) {
        switch (gl_trialparamtypes[i]) {
            case INT:
                gl_trialparams[i].val.lval = -1;
                break;
        }
    }
    
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
    int i, *tmpptr;

    gl_condarray = malloc(sizeof(int*) * (2*gl_npulse_offsets+gl_opt_middle_pulse));

    for ( i = 1; i <= 2*gl_npulse_offsets+1; ++i ) {
        if (!gl_opt_middle_pulse && i == gl_npulse_offsets+1) // don't add a pulse at time duration/2
            continue;
        tmpptr = malloc(sizeof(int) * NTRIALPARAMS);
        // think linspace(0,duration,2*npulses+3) and remove first and last entry
        tmpptr[TP_PULSE2OFFSET] = (int) (1.0*i*gl_stim_dur/(2*gl_npulse_offsets+2));
        gl_condarray[gl_nconds++] = tmpptr;
    }
#ifdef TESTING
    dprintf("gl_nconds = %d, condarray has %d pointers\n",
        gl_nconds, 2*gl_npulse_offsets+gl_opt_middle_pulse);
#endif
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
    int limit = RAND_MAX - RAND_MAX % n;
    int rnd;

    do { rnd = random(); } while ( rnd >= limit );
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
#ifdef TESTING
    dprintf("ShuffleTrials completed:");
    for ( i = 0; i < gl_nconds; dprintf(" %d", gl_condorder[i++]) );
    dprintf("\n");
#endif
}

/*
    0    1 = is_correct
0   R    L

1   L    R
^
gl_delay_pulse
*/

int save_targ_choice(int is_correct)
{
    gl_trialparams[TP_RIGHTWARDSACCCD].val.lval = !(is_correct^gl_delay_pulse); // RIGHTWARD -> 1, LEFTWARD -> 0
    gl_trialparams[TP_CORRECTCD].val.lval = gl_correct_choice = is_correct || (gl_show_probe && (gl_L > 0 || gl_M > 0));
    return 0;
}

int ttl_on(int nth_pulse)
{
    dio_on(FTTL2);

    switch (nth_pulse) {
        case 1:
            return TPULSE1CD;
        case 2:
            return TPULSE2CD;
        default:
            return TPULSE3CD;
    }
}

int reset_trial()
{
    dio_off(FTTLBYTE1);
    udp_send_to("AllOff();", "MAC");
    update_eyewin(WINDFIXN, gl_fp_x, gl_fp_y, gl_eyewin_w, gl_eyewin_h); // fixation
    update_eyewin(WINDCORR, 1000, 1000, gl_targwin_w, gl_targwin_h); // correct targ offscreen
    update_eyewin(WINDWRNG, 1000, 1000, gl_targwin_w, gl_targwin_h); // incorrect targ offscreen

    gl_delay_pulse = 0;
    gl_pulse_offset = 0;
    gl_correct_choice = 0;
    gl_mac_ready = 0;
    gl_stimon = 0;

    return 0;
}

int update_trial_counter()
{
    gl_trials_left--;
    if ( ++gl_condidx == gl_nconds ) {
        if ( ++gl_nblockscompleted >= gl_nblocks ) {
            gl_done = 1;
            dprintf("maximum number of trials reached\n");
        } else
            shuffle_conds();

        gl_condidx = 0;
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
    PARADIGM_PRELUDE("ThreePulseMonte");

    reset_trial();
    srandom(time(0));

    gl_stim_dur = gl_stim_dur_public;
    gl_firsttrial = 1;
    gl_nconds = 0;
    gl_condidx = 0;
    gl_nblockscompleted = 0;
    gl_done = 0;
    gl_trials_left = 0;

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
    {"X",             &gl_fp_x,     NP, NP, 0, ME_DEC},
    {"Y",             &gl_fp_y,     NP, NP, 0, ME_DEC},
    {"Size",          &gl_fp_size,  NP, NP, 0, ME_DEC},
    {"Window width",  &gl_eyewin_w, NP, NP, 0, ME_DEC},
    {"Window height", &gl_eyewin_h, NP, NP, 0, ME_DEC},
    {NS},
};

VLIST menu_targ[] =
{
    {"Size",                 &gl_targ_size, NP, NP, 0, ME_DEC},
    {"Crutch stim size",     &gl_crutchstim_size, NP, NP, 0, ME_DEC},
    {"Window width",         &gl_targwin_w, NP, NP, 0, ME_DEC},
    {"Window height",        &gl_targwin_h, NP, NP, 0, ME_DEC},
    {"Wrong targ contrast [0-1]", &gl_wrong_targ_contrast, NP, sanitize, ME_AFT, ME_FLOAT},
    {NS},
};

VLIST menu_expt[] =
{
    {"RF X",            &gl_rf_x, NP, NP, 0, ME_DEC},
    {"RF Y",            &gl_rf_y, NP, NP, 0, ME_DEC},
    {"Duration",        &gl_stim_dur_public,  NP, sanitize, ME_AFT, ME_DEC},
    {"# pulse offsets", &gl_npulse_offsets,   NP, NP, 0, ME_DEC},
    {"Pulse at dur/2?", &gl_opt_middle_pulse, NP, NP, 0, ME_DEC},
    {"# blocks",        &gl_nblocks,          NP, sanitize, ME_AFT, ME_DEC},
    {"coinflip override", &gl_coinflip_override, NP, NP, 0, ME_DEC},
    {"cond override (0-2)", &gl_condition_override, NP, sanitize, ME_AFT, ME_DEC},
    {"color targets?",  &gl_color_targets, NP, NP, 0, ME_DEC},
    {NS},
};

VLIST menu_gabor[] =
{
    {"Gabors?", &gl_use_gabors, NP, NP, 0, ME_DEC},
    {"Thresh mult override", &gl_thresh_mult_override, NP, NP, 0, ME_FLOAT},
    {"Temp freq", &gl_tf, NP, NP, 0, ME_FLOAT},
    {"probe L", &gl_L, NP, NP, 0, ME_FLOAT},
    {"probe M", &gl_M, NP, NP, 0, ME_FLOAT},
    {"show probe?", &gl_show_probe_public, NP, NP, 0, ME_DEC},
    {NS},
};

MENU umenus[] = {
    {"Fixation",   &menu_fp,   NP, NP, 0, NP, NP},
    {"Target",     &menu_targ, NP, NP, 0, NP, NP},
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
    {"# blocks completed", &gl_nblockscompleted},
    {"correct choice?", &gl_correct_choice},
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
    to fpOn
fpOn:
    code FPONCD
    do show_FP()
    time 2500
    to fixTime on -WD0_XY & eyeflag
    to bail
fixTime:
    code FPACQCD
    time 100 /* have this time dynamically set based on gl_stim_dur? */
    to bail on +WD0_XY & eyeflag
    to showVisual
showVisual:
    do show_visual()
    to pulse1
pulse1:
    do ttl_on(1)
    to bail on +WD0_XY & eyeflag
    to pulse1off
pulse1off: /* check init_trial_params for setting this state's `time` */
    do dio_off(FTTL2)
    to bail on +WD0_XY & eyeflag
    to pulse2
pulse2:
    do ttl_on(2)
    to bail on +WD0_XY & eyeflag
    to pulse2off
pulse2off: /* check init_trial_params for setting this state's `time` */
    do dio_off(FTTL2)
    to bail on +WD0_XY & eyeflag
    to pulse3
pulse3:
    do ttl_on(3)
    to bail on +WD0_XY & eyeflag
    to pulse3off
pulse3off:
    do dio_off(FTTL2)
    to bail on +WD0_XY & eyeflag
    to hideVisual
hideVisual:
    do hide_visual()
    to preTargsDel
preTargsDel:
    time 200 rand 500
    to bail on +WD0_XY & eyeflag
    to fpOff
fpOff:
    code FPOFFCD
    do hide_FP()
    to bail on +WD0_XY & eyeflag
    to targsOn
targsOn:
    code TARGONCD
    do show_targs()
    to reactTime
reactTime:
    time 700
    to choseCorrect on -WD1_XY & eyeflag
    to choseWrong on -WD2_XY & eyeflag
    to bail
choseCorrect:
    do save_targ_choice(1)
    to targWait
choseWrong:
    do save_targ_choice(0)
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
	time 800
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
    time 150
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
