/*
 * IsoSamp - a REX paradigm for presenting gabor stimuli that lie on iso-detection surfaces.
 *
 * ZALB Oct 2011
 *
 * Advanced feature: in the "Advanced" user menu, the user can specify a subregion of the surface
 * to pick stimuli from. This is parameterized in cylindrical coordinates: r_lo <= 1 <= r_hi
 * (where 1 means at detection threshold in most cases) defines an interval to pick an r uniformly;
 * 0 <= theta_lo <= theta_hi < pi defines a wedge in the LM-plane; and z_lo < z_hi restricts the
 * z-coordinate (for example, z could be S cone contrast or temporal frequency).
 * Note that the r of each data point will be scaled by a scalar picked uniformly in [r_lo,r_hi].
 *
 * There is a "Behavior" mode (toggle it in a user menu) that will present the stimulus either in
 * the receptive field or with a flipped x-component. Then the subject reports which side the
 * stimulus appeared with a saccade to one of two targets.
 *
 * The "use previous stimuli" flag will signal the Online code to return the same (L,M,S,TF) tuples
 * as the previous run. This is useful if stimulus generation involves some sort of randomization,
 * but you wish to reproduce the exact stimuli set. This flag is ignored if the receptive field
 * changed, and a warning will appear on the Online side.
 */

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "ldev.h"
#include "labcodes.h"
#include "UDPtransfer.h"
#include "rigconsts.h"
#include "IsoSamp.h"

// This paradigm uses global structs to hold parameters. Most other paradigms use standalone
// global variables to hold information. Using structs over variables was mostly a cosmetic choice.
// You'll need to edit the struct definition in the header file in order to add or remove variables.
// In C, all globally scoped variables are zero'd; this means the default value for any of the
// struct members is 0. To supply non-zero defaults, modify the initGlobalParams function below.
params_public gPublic;
params_private gPrivate;
params_gabor gGabor;

// eye window handles
#define WINDFIXN 0
#define WINDCORR 1
#define WINDWRNG 2

/////////////////////////////////
// Stuff for TrialParams array //
/////////////////////////////////

// Indices for gl_trialparams
#define TP_CORRECTCD 0
#define TP_STIMXCD   1
#define TP_STIMYCD   2

#define NTRIALPARAMS 14

// must put the codes in INT, FLOAT, DOUBLE order (see initTrialParams())
int gl_trialparamcodes[NTRIALPARAMS] = {CORRECTTRIALCD, STIMXCD, STIMYCD, GABFLASHTIMECD,
                                        GABPHICD, GABGAMMACD, GABTFCD, GABNSIGMASCD, GABSIGMACD, GABTHETACD, GABSFCD,
                                        GABLCD, GABMCD, GABSCD};
int gl_trialparamtypes[NTRIALPARAMS] = {INT, INT, INT, INT,
                                        FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT, FLOAT,
                                        DOUBLE, DOUBLE, DOUBLE};
ecode_struct gl_trialparams[NTRIALPARAMS];

////////////////////////////////////////////////
// Functions for opening/closing UDP channels //
////////////////////////////////////////////////

int udpOpen()
{
    rigconstants(whichrig());
    udp_open_2con(REXIP, MACIP, PLEXIP, UDPPORT);
    return 0;
}

int udpClose()
{
    udp_close();  // Closing a closed port is better than opening an open one
    return 0;
}

////////////////////////////////////////////
// Functions for sending requests via UDP //
////////////////////////////////////////////

const char FUNDFILENAME[] = "'T_cones_smj10.mat'";
const char RFXY_REQ[] = "sendToRex(udpCom, gl.bar.xy, 'double', message);";
const char PREFSFORIENT_REQ[] = "sendToRex(udpCom, [gl.grating.prefsf gl.grating.preforient], 'double', message);";
const char NEXTSTIM_REQ[] = "sendToRex(udpCom, nextStim(), 'double', message);";
const char GENSTIMS_REQ[] = "sendToRex(udpCom, genStimSet(%d,%d,%lf,%lf,%lf,[%lf %lf %lf %lf %lf %lf],[%d %d],%d,%d), 'integer', message);";
const char MODELREV_REQ[] = "sendToRex(udpCom, gl.model_revision, 'integer', message);";
const char MODELSTR_REQ[] = "sendToRex(udpCom, gl.model{1}, 'integer', message);";
const char FRAMERATE_REQ[] = "sendToRex(udpCom, gl.framerate, 'double', message);";
const char PIXPERDEG_REQ[] = "sendToRex(udpCom, gl.pixperdeg, 'double', message);";
const char BKGNDRGB_REQ[] = "sendToRex(udpCom, gl.bkgndrgb, 'double', message);";
const char GAMMATABLE_REQ[] = "sendToRex(udpCom, gl.cal.gammaTable, 'double', message);";
const char FUNDAMENTALS_REQ[] = "sendToRex(udpCom, gl.cal.fundamentals, 'double', message);";
const char MONSPD_REQ[] = "sendToRex(udpCom, gl.cal.monSpd, 'double', message);";
const char MODELPARAMS_REQ[] = "sendToRex(udpCom, gl.modelout, 'double', message);";

#define DECL_UDPSEND_FUNC( funcname, reqstr, dest ) int funcname##() { udp_send_to(reqstr, dest); return 0; }
DECL_UDPSEND_FUNC( startSlavePrgm, "IsoSampSlave();", "MAC" )
DECL_UDPSEND_FUNC( stopSlavePrgm, "return;", "MAC" )
DECL_UDPSEND_FUNC( macAllOff, "AllOff();", "MAC" )
DECL_UDPSEND_FUNC( startOnlinePrgm, "IsoSampOnline();", "PLEXON" )
DECL_UDPSEND_FUNC( stopOnlinePrgm, "return;", "PLEXON" )
DECL_UDPSEND_FUNC( getModelRev, MODELREV_REQ, "PLEXON" )
DECL_UDPSEND_FUNC( getModelStr, MODELSTR_REQ, "PLEXON" )
DECL_UDPSEND_FUNC( getFrameRate, FRAMERATE_REQ, "MAC" )
DECL_UDPSEND_FUNC( getPixperDeg, PIXPERDEG_REQ, "MAC" )
DECL_UDPSEND_FUNC( getBkgndrgb, BKGNDRGB_REQ, "MAC" )
DECL_UDPSEND_FUNC( getGammaTable, GAMMATABLE_REQ, "MAC" )
DECL_UDPSEND_FUNC( getFundamentals, FUNDAMENTALS_REQ, "MAC" )
DECL_UDPSEND_FUNC( getMonSpd, MONSPD_REQ, "MAC" )
DECL_UDPSEND_FUNC( getModelParams, MODELPARAMS_REQ, "PLEXON" )
#undef DECL_UDPSEND_FUNC

///////////////////////////
//  Mac commands via UDP //
///////////////////////////

int getPrefSfOrient()
{
    static int done_once = 0;
    if (!done_once) {
        done_once = 1;
        udp_send_to(PREFSFORIENT_REQ, "MAC");
    }
    return 0;
}

int getRFxy()
{
    if ( !gPublic.rf_x && !gPublic.rf_y )
        udp_send_to(RFXY_REQ, "MAC");
    return 0;
}

int setupMac()
{
    char buf[256];
    gPrivate.mac_ready = 0;
    sprintf(buf, "InitDisplay(%d, %d, %s, %s)",
            GL_MONDIST, GL_SCREENWIDTH, GL_CALFILENAME, FUNDFILENAME);
    udp_send_to(buf, "MAC");
    return 0;
}

int showFP()
{
    char buf[256];
    sprintf(buf, "ShowFP(%d, %d, %d, %d, %d, %d)",
            gPublic.fp_x, gPublic.fp_y, gPublic.fp_size,
            gPublic.fp_rgb[0], gPublic.fp_rgb[1], gPublic.fp_rgb[2]);
    udp_send_to(buf, "MAC");
    gPrivate.fpon = 1;
    update_eyewin(WINDFIXN, gPublic.fp_x, gPublic.fp_y, gPublic.eyewin_w, gPublic.eyewin_h);
    return 0;
}

int hideFP()
{
    if ( gPrivate.fpon ) {
        udp_send_to("HideFP()", "MAC");
        gPrivate.fpon = 0;
    }
    return 0;
}

int showStim()
{
    char buf[256];
    int x = gl_trialparams[TP_STIMXCD].val.lval, y = gl_trialparams[TP_STIMYCD].val.lval;

    sprintf(buf, "ShowStim(%d, %d, %lf, %lf, %lf)",
        x, y, gGabor.lms[0], gGabor.lms[1], gGabor.lms[2]);
    udp_send_to(buf, "MAC");
    return 0;
}


int showTargs()
{
    char buf[64];
    double scale;
    int desired_distance = 20;
    // At this point, stim_x could have a flipped sign if the behavior flag is true
    int x = gl_trialparams[TP_STIMXCD].val.lval, y = gl_trialparams[TP_STIMYCD].val.lval;

    // Make sure the targets are `desired_distance` from the fixation point
    scale = sqrt(1.0*desired_distance*desired_distance/(x*x + y*y));
    x = my_round(scale*x); y = my_round(scale*y);

    update_eyewin(WINDCORR, x, y, gPublic.targwin_w, gPublic.targwin_h);
    update_eyewin(WINDWRNG,-x, y, gPublic.targwin_w, gPublic.targwin_h);
    sprintf(buf, "ShowTargs(%d, %d, %d)", x, y, gPublic.targ_size);
    udp_send_to(buf, "MAC");
    return 0;
}

// This call is the main initialization call to the Online code. Most of these parameters in the
// sprintf below are passed, in order, to the module script that generates stimuli. If the order
// changes here, then so does the order of `varargin` on the Online side. The contents of varargin
// in the module script are the arguments passed here beginning with the third argument to sprintf
// below.
int genStimSet()
{
    char buf[512];
    gPrivate.nstims_block = 0;
    sprintf(buf, GENSTIMS_REQ, gPublic.target_nstims, gPublic.prev_stims, gGabor.sf, gGabor.theta,
        gGabor.tf,gPublic.r[0], gPublic.r[1], gPublic.th[0], gPublic.th[1], gPublic.z[0],
        gPublic.z[1], gPublic.rf_x, gPublic.rf_y, gPublic.special_case, gPublic.blanks);
    udp_send_to(buf, "PLEXON");
    return 0;
}

// Get the next stimulus at the top of every trial
int getNextStim()
{
    gPrivate.plx_ready = 0;
    udp_send_to(NEXTSTIM_REQ, "PLEXON");
    return 0;
}

//////////////////////////
// Getting UDP commands //
//////////////////////////

int checkMsg()
{
    int i, tmpi, messageAvailable;
    long model_str[32], lmtf_bytes[] = {76, 77, 84, 70}; // "LMTF" as ASCII bytes
    char inMsg[REXBUFF]; // 8500
    static message_struct message;
    double tmpd, dvals[800];

    inMsg[0] = '\0';
    messageAvailable = udp_check(0);

    if ( messageAvailable ) {
        udp_read(inMsg, REXBUFF);
        prepareMsg(&message, inMsg);

        if ( !strcmp(message.header, "MACSTIMON") ) {
            gPrivate.stimon = 1;
        } else if ( !strcmp(message.header, "MACSTIMOFF") ) {
            gPrivate.stimon = 0;
        } else if ( !strcmp(message.header, "STIMIDXUPDATED") ) {
            gPrivate.plx_ready = 1;
        } else if ( !strcmp(message.header, "ONLINEINIT") ) {
            gPrivate.plx_ready = 1;
        } else if ( !strcmp(message.header, "PARSEDHEADER") ) {
            gPrivate.plx_ready = 1;
        } else if ( !strcmp(message.header, "DISPLAYINIT") ) {
            gPrivate.mac_ready = 1;
        } else if ( !strcmp(message.header, "GABORREADY") ) {
            gPrivate.mac_ready = 1;
        } else if ( !strcmp(message.header, RFXY_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            if ( message.size && (dvals[0] || dvals[1]) ) {
                gPublic.rf_x = my_round(10. * dvals[0]);
                gPublic.rf_y = my_round(10. * dvals[1]);
            }
        } else if ( !strcmp(message.header, PREFSFORIENT_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            gGabor.sf = (float) dvals[0]; gGabor.theta = (float) dvals[1];
        } else if ( !strcmp(message.header, NEXTSTIM_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            for ( i = 0; i < 3; ++i ) gGabor.lms[i] = dvals[i];
            gGabor.tf = (float) dvals[i];
            gPrivate.plx_ready = 1;
        } else if ( strstr(message.header, "genStimSet") ) {
            gPrivate.nstims_block = hex2long(message.contents, NULL);
            gPrivate.trials_remaining_block = gPrivate.nstims_block;
            dprintf("-> Generated %d stimuli!\n", gPrivate.nstims_block);
            if (!gPrivate.nstims_block) // something went wrong, gracefully halt the paradigm
                gPrivate.done = 1;
        } else if ( !strcmp(message.header, MODELREV_REQ) ) {
            tmpi = hex2long(message.contents, NULL);
            sendToPlexDataStream("long", (void *) &tmpi, message.size, MODELREVCD, 0);
        } else if ( !strcmp(message.header, MODELSTR_REQ) ) {
            hex2long(message.contents, (void *) &model_str);
            if ( !memcmp(model_str, lmtf_bytes, 4*sizeof(long))) {
                gGabor.sf = 1.0f; // always force the sf to 1 for LMTF
            }
            sendToPlexDataStream("long", (void *) &model_str, message.size, MODELSTRCD, 1);
        } else if ( !strcmp(message.header, BKGNDRGB_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            sendToPlexDataStream("double", (void *) &dvals, message.size, BKGNDRGBCD, 1);
        } else if ( !strcmp(message.header, FRAMERATE_REQ) ) {
            tmpd = hex2double(message.contents, NULL);
            sendToPlexDataStream("double", (void *) &tmpd, message.size, FRAMERATECD, 0);
        } else if ( !strcmp(message.header, PIXPERDEG_REQ) ) {
            tmpd = hex2double(message.contents, NULL);
            sendToPlexDataStream("double", (void *) &tmpd, message.size, PIXPERDEGCD, 0);
        } else if ( !strcmp(message.header, GAMMATABLE_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            sendToPlexDataStream("double", (void *) &dvals, message.size, GAMMATABLECD, 1);
        } else if ( !strcmp(message.header, FUNDAMENTALS_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            sendToPlexDataStream("double", (void *) &dvals, message.size, FUNDAMENTALSCD, 1);
        } else if ( !strcmp(message.header, MONSPD_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            sendToPlexDataStream("double", (void *) &dvals, message.size, MONSPDCD, 1);
            gPrivate.plx_ready = 0; // the Online program will ping back when the header is parsed fully
            ec_send_code_tagged(HEADERREADYCD, VALOFFSET+1);
            // HEADERREADYCD tells the online program that it's safe to
            // parse the header (i.e., REX has sent all relevant variables)
        } else if ( !strcmp(message.header, MODELPARAMS_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            sendToPlexDataStream("double", (void *) &dvals, message.size, MODELPARAMSCD, 1);
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////
// Functions for dealing with the header and trial parameters //
////////////////////////////////////////////////////////////////

int dropHeader()
{
    ec_send_code_tagged(PARADIGMIDENTCD , VALOFFSET + PARADIGMID);
    ec_send_code_tagged(FIXXCD          , VALOFFSET + gPublic.fp_x);
    ec_send_code_tagged(FIXYCD          , VALOFFSET + gPublic.fp_y);
    ec_send_code_tagged(FPSIZECD        , VALOFFSET + gPublic.fp_size);
    ec_send_code_tagged(FPRCD           , VALOFFSET + gPublic.fp_rgb[0]);
    ec_send_code_tagged(FPGCD           , VALOFFSET + gPublic.fp_rgb[1]);
    ec_send_code_tagged(FPBCD           , VALOFFSET + gPublic.fp_rgb[2]);
    ec_send_code_tagged(EYEWINXCD       , VALOFFSET + gPublic.eyewin_w);
    ec_send_code_tagged(EYEWINYCD       , VALOFFSET + gPublic.eyewin_h);
    ec_send_code_tagged(RFXCD           , VALOFFSET + gPublic.rf_x);
    ec_send_code_tagged(RFYCD           , VALOFFSET + gPublic.rf_y);
    ec_send_code_tagged(NUMSTIMSCD      , VALOFFSET + gPrivate.nstims_block);

    sendToPlexDataStream("float", (void *) &gPublic.r, 2, RSCALESCD, 1);
    sendToPlexDataStream("float", (void *) &gPublic.th, 2, THBOUNDSCD, 1);
    sendToPlexDataStream("float", (void *) &gPublic.z, 2, ZBOUNDSCD, 1);

    gPrivate.first_trial = 0;
    return 0;
}

// This function won't work unless the codes in gl_trialparam* are ordered in a INT, FLOAT, DOUBLE fashion
int initTrialParams()
{
// This way of loading in the trial params is finicky and stupid; take care when modifying this function
// You need to match up the identities of gl_trialparam* so that the code below populates gl_trialparams
// in the correct order. 
    int i, ni_vals = 4, nf_vals = 7;
    // if gPublic.behavior != 0, then flip a coin to see if we flip the sign on rf_x
    int i_vals[] = {0, gPublic.behavior ? (my_randint(2) ? gPublic.rf_x : -gPublic.rf_x) : gPublic.rf_x, gPublic.rf_y, gGabor.flash_time};
    float f_vals[] = {gGabor.phi, gGabor.gamma, gGabor.tf, gGabor.num_sigmas, gGabor.sigma, gGabor.theta, gGabor.sf};
    double d_vals[] = {gGabor.lms[0], gGabor.lms[1], gGabor.lms[2]};

    gPrivate.plx_ready = 0;

    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
        gl_trialparams[i].casttype = gl_trialparamtypes[i];
        gl_trialparams[i].nelements = 1;

        if ( gl_trialparamtypes[i] == INT )
            gl_trialparams[i].val.lval = i_vals[i];
        else if ( gl_trialparamtypes[i] == FLOAT )
            gl_trialparams[i].val.fval = f_vals[i-ni_vals];
        else if ( gl_trialparamtypes[i] == DOUBLE )
            gl_trialparams[i].val.dval = d_vals[i-ni_vals-nf_vals];
    }
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

        sendToPlexDataStream(typebuffer,
                            (void *) &gl_trialparams[i].val.lval,
                            gl_trialparams[i].nelements,
                            gl_trialparams[i].flagcodeval, 0);
    }

    return 0;
}

/////////////////////////
//   REX setup stuff   //
/////////////////////////

// Here's where you set up the default values for the struct members
void initGlobalParams()
{
    static int firstTime = 1;
    if ( firstTime ) { // The if branch is called once at the very beginning
        firstTime = 0;

        srand(time(NULL));

        gPublic.eyewin_w = 5; gPublic.eyewin_h = 5;
        gPublic.targwin_w = 13; gPublic.targwin_h = 13;
        gPublic.targ_size = 2;
        gPublic.fp_size = 2; gPublic.nblocks_max = 10;

        gPublic.target_nstims = 30;
        gPublic.blanks = 1;
        gPublic.r[0] = 1.0f; gPublic.r[1] = 1.5f;

        gGabor.flash_time = 666;
        gGabor.gamma = 1.0f; gGabor.tf = 3.0f;
        gGabor.sigma = 0.15f; gGabor.num_sigmas = 2.0f;
        gGabor.sf = 1.0f;

        // not listed here? defaults to zero
    } else {
        // don't clear out gPublic, gGabor
        macAllOff();
        memset(&gPrivate, 0, sizeof(gPrivate));
    }
    gPrivate.first_trial = 1;
}

void rinitf()
{
    PARADIGM_PRELUDE("IsoSamp");

    resetTrial();
    initGlobalParams();

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

/////////////////////////////////////////
// Functions specific to this paradigm //
/////////////////////////////////////////

int prepareGabor()
{
    char buf[256];
    gPrivate.mac_ready = 0;
    sprintf(buf, "PrepareGabor(%d, %lf, %lf, %f, %f, %f, %f, %f)",
        gGabor.flash_time, gGabor.theta, gGabor.sf,
            gGabor.phi, gGabor.sigma, gGabor.gamma, gGabor.tf, gGabor.num_sigmas);
    udp_send_to(buf, "MAC");
    return 0;
}

int genStimCheck()
{
    return gPublic.target_nstims > 0;
}

/*
 * Copyright (c) 2008, Damien Miller <djm@openbsd.org>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

// generate a random integer on [0,upper_bound) without modulo bias
// https://github.com/libressl-portable/openbsd/blob/master/src/lib/libc/crypt/arc4random_uniform.c
static unsigned int my_randint(unsigned int upper_bound)
{
    unsigned int r, min;

    if (upper_bound < 2)
        return 0;

    min = -upper_bound % upper_bound;

    for (;;) {
        r = rand();
        if (r >= min)
            break;
    }

    return r % upper_bound;
}

// We need to handshake here because the EOT code may not be read and the Online program won't
// increment the stimulus index. Without this the same stimulus will appear back-to-back (sometimes).
int updateTrialCounts()
{
    gPrivate.plx_ready = 0; // so we get an accurate handshake at the end of the state set

    if ( --gPrivate.trials_remaining_block <= 0 ) {
        if ( ++gPrivate.nblocks_done >= gPublic.nblocks_max ) {
            dprintf("[!] we're done here\n");
            gPrivate.done = 1;
        } else {
            gPrivate.trials_remaining_block = gPrivate.nstims_block;
        }
    }

    return REQSTIMUPDATECD; // cue the Online program to update the stim idx
}

int save_targ_choice(int is_correct)
{
    gl_trialparams[TP_CORRECTCD].val.lval = gPrivate.correct_choice = is_correct;
    return 0;
}

int resetTrial()
{
    macAllOff();

    update_eyewin(WINDFIXN, gPublic.fp_x, gPublic.fp_y, gPublic.eyewin_w, gPublic.eyewin_h); // fixation
    update_eyewin(WINDCORR, 1000, 1000, gPublic.targwin_w, gPublic.targwin_h); // correct targ offscreen
    update_eyewin(WINDWRNG, 1000, 1000, gPublic.targwin_w, gPublic.targwin_h); // incorrect targ offscreen

    gPrivate.stimon = 0;
    gPrivate.fpon = 0;
    gPrivate.plx_ready = 0;
    gPrivate.mac_ready = 0;
    gPrivate.correct_choice = 0;

    return 0;
}

////////////////////////
// Menu declarations  //
////////////////////////

VLIST menu_fp[] =
{
    {"X",             &gPublic.fp_x,      NP, NP, 0, ME_DEC},
    {"Y",             &gPublic.fp_y,      NP, NP, 0, ME_DEC},
    {"Size",          &gPublic.fp_size,   NP, NP, 0, ME_DEC},
    {"Red",           &gPublic.fp_rgb[0], NP, NP, 0, ME_DEC},
    {"Green",         &gPublic.fp_rgb[1], NP, NP, 0, ME_DEC},
    {"Blue",          &gPublic.fp_rgb[2], NP, NP, 0, ME_DEC},
    {"Window width",  &gPublic.eyewin_w,  NP, NP, 0, ME_DEC},
    {"Window height", &gPublic.eyewin_h,  NP, NP, 0, ME_DEC},
    {NS},
};

VLIST menu_gabor[] =
{
    {"X",              &gPublic.rf_x,        NP, NP, 0, ME_DEC},
    {"Y",              &gPublic.rf_y,        NP, NP, 0, ME_DEC},
    {"Theta",          &gGabor.theta,        NP, NP, 0, ME_FLOAT},
    {"SF",             &gGabor.sf,           NP, NP, 0, ME_FLOAT},
    {"Phi",            &gGabor.phi,          NP, NP, 0, ME_FLOAT},
    {"TF",             &gGabor.tf,           NP, NP, 0, ME_FLOAT},
    {"Gamma",          &gGabor.gamma,        NP, NP, 0, ME_FLOAT},
    {"Sigma",          &gGabor.sigma,        NP, NP, 0, ME_FLOAT},
    {"# sigmas",       &gGabor.num_sigmas,   NP, NP, 0, ME_FLOAT},
    {NS},
};

// At the top of this spot file is a detailed comment above these variables
VLIST menu_advanced[] =
{
    {"Lower r scale", &gPublic.r[0], NP, NP, 0, ME_FLOAT},
    {"Upper r scale", &gPublic.r[1], NP, NP, 0, ME_FLOAT},
    {"Lower theta bound", &gPublic.th[0], NP, NP, 0, ME_FLOAT},
    {"Upper theta bound", &gPublic.th[1], NP, NP, 0, ME_FLOAT},
    {"Lower z bound", &gPublic.z[0], NP, NP, 0, ME_FLOAT},
    {"Upper z bound", &gPublic.z[1], NP, NP, 0, ME_FLOAT},
    {NS},
};

VLIST menu_stimuli[] =
{
    {"Target # stims", &gPublic.target_nstims, NP, NP, 0, ME_DEC},
    {"Max # blocks", &gPublic.nblocks_max, NP, NP, 0, ME_DEC},
    {"Special case #", &gPublic.special_case, NP, NP, 0, ME_DEC},
    {"Blank trials?", &gPublic.blanks, NP, NP, 0, ME_DEC},
    {NS},
};

VLIST menu_behavior[] =
{
    {"Behavior?", &gPublic.behavior, NP, NP, 0, ME_DEC},
    {"Use prev. stims?", &gPublic.prev_stims, NP, NP, 0, ME_DEC},
    {"Target width", &gPublic.targwin_w, NP, NP, 0, ME_DEC},
    {"Target height", &gPublic.targwin_h, NP, NP, 0, ME_DEC},
    {"Target size", &gPublic.targ_size, NP, NP, 0, ME_DEC},
    {NS},
};

MENU umenus[] = {
    {"Fixation", &menu_fp,       NP, NP, 0, NP, NP},
    {"Behavior", &menu_behavior, NP, NP, 0, NP, NP},
    {"Gabor",    &menu_gabor,    NP, NP, 0, NP, NP},
    {"Stimuli",  &menu_stimuli,  NP, NP, 0, NP, NP},
    {"Advanced", &menu_advanced, NP, NP, 0, NP, NP},
    {NS},
};

VLIST state_vl[] = {{NS}};

char hm_sv_vl[] = "";

RTVAR rtvars[] = {
    {"Trials left", &gPrivate.trials_remaining_block},
    {"Blocks done", &gPrivate.nblocks_done},
    {"", 0},
};

////////////////////
// The state sets //
////////////////////

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
        to endOnline
    endOnline:
        do stopOnlinePrgm()
        time 5
        to startSlave
    startSlave:
        do startSlavePrgm()
        time 200
        to initMac
    initMac:
        do setupMac()
        to startOnline on 1 = gPrivate.mac_ready
    startOnline:
        do startOnlinePrgm()
        to getRFxy on 1 = gPrivate.plx_ready
    getRFxy:
        do getRFxy()
        time 10
        to getPrefSfOrient
    getPrefSfOrient:
        do getPrefSfOrient()
        time 10
        to plxStart on -PSTOP & softswitch /* don't drop stuff until the user is ready */
    pause0:
        to pauseLoop on 1 = gPrivate.done
        to pause1 on +PSTOP & softswitch
        to getNextStim
    pause1:
        do hideFP()
        to getNextStim on -PSTOP & softswitch
    pauseLoop:
        do hideFP()
        time 1000
        to pause0
    getNextStim:
        do getNextStim()
        to prepareGabor on 1 = gPrivate.plx_ready
    prepareGabor:
        do prepareGabor()
        to initTrial on 1 = gPrivate.mac_ready
    initTrial:
        do initTrialParams()
        time 200
        to plxStart
    plxStart:
        do dio_on(PLXREC)
        time 200
        to sendHeaderReqs on 1 = gPrivate.first_trial
        to fpOn on 0 = gPrivate.fpon
        to error on +WD0_XY & eyeflag
        to showStim
    sendHeaderReqs:
        to getPixperDeg
    getPixperDeg:
        do getPixperDeg()
        time 20
        to getBkgndrgb
    getBkgndrgb:
        do getBkgndrgb()
        time 20
        to getFrameRate
    getFrameRate:
        do getFrameRate()
        time 20
        to getFundamentals
    getFundamentals:
        do getFundamentals()
        time 1000
        to getMonSpd
    getMonSpd:
        do getMonSpd()
        time 2000
        to getGammaTable
    getGammaTable:
        do getGammaTable()
        time 3000
        to dropHeader
    dropHeader:
        do dropHeader()
        time 300
        to plexWait
    plexWait:
        time 1000
        to onlineParseHeader
    onlineParseHeader:
        time 5000
        to genStimWait on 1 = gPrivate.plx_ready
        to first
    genStimWait: /* wait for user to input a target # of stims via the menu */
        to genStimSet on 1 % genStimCheck
    genStimSet:
        do genStimSet()
        to getModelRev on 0 < gPrivate.nstims_block
        to pauseLoop on 1 = gPrivate.done
    getModelRev:
        do getModelRev()
        time 20
        to getModelStr
    getModelStr:
        do getModelStr()
        time 100
        to getModelParams
    getModelParams:
        do getModelParams()
        time 200
        to pause0
    fpOn:
        code FPONCD
        do showFP()
        time 2000
        to preStimDel on -WD0_XY & eyeflag
        to error
    preStimDel:
        code FPACQCD
        time 1500
        to showStim
    showStim:
        do showStim()
        to stimCode on 1 = gPrivate.stimon
    stimCode:
        code STIMONCD
        to error on +WD0_XY & eyeflag
        to postStimDel on 0 = gPrivate.stimon
    postStimDel:
        code STIMOFFCD
        time 100
        to error on +WD0_XY & eyeflag
        to behaviorFork
    behaviorFork:
        to error on +WD0_XY & eyeflag
        to fpOff on 0 < gPublic.behavior
        to rewFixOn
    rewFixOn:
        code REWCD
        do dio_on(REW)
        time 60
        to rewFixOff
    rewFixOff:
        do dio_off(REW)
        to incTrialCounts
    rewBehavOn:
        code REWCD
        do dio_on(REW)
        time 120
        to rewBehavOff
    rewBehavOff:
        do dio_off(REW)
        to resetTrial
    fpOff:
        code FPOFFCD
        do hideFP()
        time 100 rand 500
        to error on +WD0_XY & eyeflag
        to targsOn
    targsOn:
        code TARGONCD
        do showTargs()
        to reactTime
    reactTime:
        time 700
        to choseCorrect on -WD1_XY & eyeflag
        to choseWrong on -WD2_XY & eyeflag
        to error
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
        to evalCorrect on 1 = gPrivate.correct_choice
        to evalIncorrect
    evalCorrect:
        to rewBehavOn on -WD1_XY & eyeflag
        to error
    evalIncorrect:
        to resetTrial on -WD2_XY & eyeflag
        to error
    resetTrial:
        do resetTrial()
        to incTrialCounts
    incTrialCounts:
        do updateTrialCounts()
        to dropTrialParams on 1 = gPrivate.plx_ready /* wait until stimidx updates */
    error:
        code ABORTCD
        do resetTrial()
        time 20
        to dropTrialParams
    dropTrialParams:
        time 200
        do dropTrialParams()
        to markTrialEnd
    markTrialEnd:
        code EOT
        to plxStop
    plxStop:
        do dio_off(PLXREC)
        time 50
        to pause0
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
