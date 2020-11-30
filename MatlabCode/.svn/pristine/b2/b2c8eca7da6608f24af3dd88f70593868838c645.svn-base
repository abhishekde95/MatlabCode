/*
 * Paradigm for optically stimulating V1 during the presentation of a visual stimulus.
 * Does stimulation of the V1->SC axons change visual responses in the SC?
 *
 * In the _contrastresponse version, there is only the Fix condition with varying optical stimulation
 * intensities and varying contrast of a visual target.
 * There is no behavior required: we're looking for a pure electrophysiological effect.
 *
 * Original framework: GDLH 5/6/11
 * First production version: ZALB 6/7/11
 * Second production version: ZALB 10/26/11
 * Third production version: ZALB 02/08/12
 * _interfreq version: ZALB 01/29/13
 * _contrastresponse version: GDLH 2/16/13
 */

char gl_fundfilename [] = "'T_cones_smj10.mat'";

//#define TESTING

#include <stdlib.h>
#include <math.h>
#include "ldev.h"
#include "labcodes.h"
#include "UDPtransfer.h"
#include "rigconsts.h"
#include "FixStim.h"

#define MAXNCONDS 100

/********************/
/* Request strings   */
/********************/
const char RFXY_REQ[] = "sendToRex(udpCom, gl.bar.xy, 'double', message);";
const char GAMMA_REQ[] = "sendToRex(udpCom, gl.cal.gammaTable, 'double', message);";
const char FUNDAMENTALS_REQ[] = "sendToRex(udpCom, gl.cal.fundamentals, 'double', message);";
const char MONSPD_REQ[] = "sendToRex(udpCom, gl.cal.monSpd, 'double', message);";
const char BKGNDRGB_REQ[] = "sendToRex(udpCom, gl.bkgndrgb, 'double', message);";

//////////////////////////////////
// Global variable declarations //
//////////////////////////////////

// Parameters that the user can change
int gl_fpx = 0;
int gl_fpy = 0;
#ifdef TESTING
int gl_eyewinx = 150;
int gl_eyewiny = 150;
#else
int gl_eyewinx = 7;
int gl_eyewiny = 7;
#endif
int gl_fpsize = 2;
int gl_fprgb[] = {0, 0, 0};
int gl_ntrlspercond = 5;
int gl_rfx = 0;
int gl_rfy = 0;
int gl_max_laser_power = 255;
float gl_maxcontrast = 0.2;
int gl_ncontrasts = 4; // ncontrasts = gl_ncontrasts (one polarity only) + 1 (no visual stim)
float gl_ccmults[] = {1, 1, 1};

// Housekeeping globals
int gl_fpon = 0;
int gl_done = 0;
int gl_firsttrial = 1;
int gl_laser_power = 0;
int gl_nconds = 0;
int gl_condidx = 0; // index into gl_condarray
int gl_nblockscompleted = 0; // this counter increments to gl_ntrlspercond
int *gl_condorder;
double **gl_condarray;
// New variables
int gl_stimheight = 20;
int gl_stimwidth = 2;

/////////////////////////////////
// Stuff for TrialParams array //
/////////////////////////////////

// Indices for gl_trialparams and gl_condarray
#define TP_LASERPOWER 0
#define TP_CONTRAST   1
#define TP_HEIGHT     2
#define TP_WIDTH      3
#define TP_TARGSHOWN  4
#define TP_LCC        5
#define TP_MCC        6
#define TP_SCC        7
#define NTRIALPARAMS  8

int gl_trialparamcodes[NTRIALPARAMS] = {LASERPOWERCD, CONTRASTCD, STIMHEIGHTCD, STIMWIDTHCD,
                                         TARGSHOWNCD, LCCCD, MCCCD, SCCCD};
int gl_trialparamtypes[NTRIALPARAMS] = {INT, FLOAT, INT, INT, INT, FLOAT, FLOAT, FLOAT};
ecode_struct gl_trialparams[NTRIALPARAMS];

//////////////////////
//  Command strings //
//////////////////////

const char STARTSLAVE_ALERT[] = "FixStim;";
const char RETURN_ALERT[] = "return;";

////////////////////////////////////////////////
// Functions for opening/closing UDP channels //
////////////////////////////////////////////////

int udpopen()
{
    rigconstants(whichrig());
    udp_open_2con(REXIP, MACIP, PLEXIP, UDPPORT);
    return 0;
}

int udpclose()
{
    udp_close();
    return 0;
}

////////////////////////////////////////////
// Functions for sending requests via UDP //
////////////////////////////////////////////

int getBkgndrgb()
{
    udp_send_to(BKGNDRGB_REQ, "MAC");
    return 0;
}

int getGammaTable()
{
    udp_send_to(GAMMA_REQ, "MAC");
    return 0;
}

int getFundamentals()
{
    udp_send_to(FUNDAMENTALS_REQ, "MAC");
    return 0;
}

int getMonSpd()
{
    udp_send_to(MONSPD_REQ, "MAC");
    return 0;
}

int getRFxy()
{
    if ( !gl_rfx && !gl_rfy )
        udp_send_to(RFXY_REQ, "MAC");
    return 0;
}

///////////////////////////
//  Mac commands via UDP //
///////////////////////////

int SetupMac()
{
    char buf[256];
    sprintf(buf, "InitDisplay(%d, %d, %s, %s)",
            GL_MONDIST, GL_SCREENWIDTH, GL_CALFILENAME, gl_fundfilename);
    udp_send_to(buf, "MAC");
    return 0;
}

int StartSlavePrgm()
{
    udp_send_to(STARTSLAVE_ALERT, "MAC");
    return 0;
}

int EndSlavePrgm()
{
    udp_send_to(RETURN_ALERT, "MAC");
    return 0;
}

//////////////////////////
// Getting UDP commands //
//////////////////////////

int CheckMsg()
{
    char inMsg[REXBUFF];
    static message_struct message;
    int messageAvailable;
    double doublearray[790];

    inMsg[0] = '\0';
    messageAvailable = udp_check(0);

    if ( messageAvailable ) {
        udp_read(inMsg, REXBUFF);
        prepareMsg(&message, inMsg);

        if ( !strcmp(message.header, BKGNDRGB_REQ) ) {
            hex2double(message.contents, (void *)&doublearray);
            sendToPlexDataStream("double", (void *)&doublearray, message.size, BKGNDRGBCD, 1);
        } else if ( !strcmp(message.header, GAMMA_REQ) ) {
            hex2double(message.contents, (void *)&doublearray);
            sendToPlexDataStream("double", (void *)&doublearray, message.size, GAMMATABLECD, 1);
        } else if ( !strcmp(message.header, FUNDAMENTALS_REQ) ) {
            hex2double(message.contents, (void *)&doublearray);
            sendToPlexDataStream("double", (void *)&doublearray, message.size, FUNDAMENTALSCD, 1);
        } else if ( !strcmp(message.header, MONSPD_REQ) ) {
            hex2double(message.contents, (void *)&doublearray);
            sendToPlexDataStream("double", (void *)&doublearray, message.size, MONSPDCD, 1);
        } else if ( !strcmp(message.header, RFXY_REQ) ) {
            hex2double(message.contents, (void *)&doublearray);
            if ( message.size && doublearray[0] || doublearray[1] ) {
                gl_rfx = my_round(10 * doublearray[0]);
                gl_rfy = my_round(10 * doublearray[1]);
            }
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////
// Functions for dealing with the header and trial parameters //
////////////////////////////////////////////////////////////////

int DropHeader()
{
    ec_send_code_tagged(PARADIGMIDENTCD, VALOFFSET+PARADIGMID);
    ec_send_code_tagged(FIXXCD, VALOFFSET+gl_fpx);
    ec_send_code_tagged(FIXYCD, VALOFFSET+gl_fpy);
    ec_send_code_tagged(FPSIZECD, VALOFFSET+gl_fpsize);
    ec_send_code_tagged(FPRCD, VALOFFSET+gl_fprgb[0]);
    ec_send_code_tagged(FPGCD, VALOFFSET+gl_fprgb[1]);
    ec_send_code_tagged(FPBCD, VALOFFSET+gl_fprgb[2]);
    ec_send_code_tagged(EYEWINXCD, VALOFFSET+gl_eyewinx);
    ec_send_code_tagged(EYEWINYCD, VALOFFSET+gl_eyewiny);
    ec_send_code_tagged(RFXCD, VALOFFSET+gl_rfx);
    ec_send_code_tagged(RFYCD, VALOFFSET+gl_rfy);
    ec_send_code_tagged(NTRLSPERCONDCD, VALOFFSET+gl_ntrlspercond);

    gl_firsttrial = 0;
    dio_on(FTTL3); // turn on the TTL for the laser, control intensity with analog input
    return 0;
}

int InitTrialParams()
{
    int i;
    double *tmpptr;

    tmpptr = gl_condarray[gl_condorder[gl_condidx]];

    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
        gl_trialparams[i].casttype = gl_trialparamtypes[i];
        gl_trialparams[i].nelements = 1;
        switch ( gl_trialparamtypes[i] ) {
            case DOUBLE:
                gl_trialparams[i].val.dval = tmpptr[i];
                break;
            case FLOAT:
                gl_trialparams[i].val.fval = (float) tmpptr[i];
                break;
            case INT:
                gl_trialparams[i].val.lval = (int) tmpptr[i];
                break;
        }
    }
    gl_laser_power = gl_trialparams[TP_LASERPOWER].val.lval;
    return 0;
}

int DropTrialParams()
{
    int i, ncodes = 0;
    char typebuffer[8];
    float mspercode = 0.5f; // Conservative; actual value is closer to 0.33

    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        if ( gl_trialparams[i].casttype == FLOAT ) {
            sprintf(typebuffer, "float");
            ncodes = ncodes + 5;
        } else if ( gl_trialparams[i].casttype == DOUBLE ) {
            sprintf(typebuffer, "double");
            ncodes = ncodes + 9;
        } else { // INT
            sprintf(typebuffer, "int");
            ncodes = ncodes + 2;
        }
        sendToPlexDataStream(typebuffer,
            (void *) &gl_trialparams[i].val.lval,
            gl_trialparams[i].nelements,
            gl_trialparams[i].flagcodeval, 0);
    }
    set_times("dropTrialParams", my_ceil(ncodes * mspercode), 0);
    return 0;
}

/////////////////////////
//   REX setup stuff   //
/////////////////////////

void rinitf()
{
    PARADIGM_PRELUDE("FixStim_contrastresponse");

    udp_send_to("AllOff()", "MAC");

    // WIND0 is used for testing eye position
    wd_src_pos(WIND0, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_check(WIND0, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);
    wd_cntrl(WIND0, WD_ON); // don't use WD_OFF, just move the eye window around if required

    center_eyewin(WIND0, 1); // size and position the eye window

    dio_off(FTTLBYTE1); // turn off TTLs 1A - 8A (the first byte)

    // default global variables
    gl_firsttrial = 1;
    gl_nconds = 0;
    gl_condidx = 0;
    gl_nblockscompleted = 0;
    gl_done = 0;

    FreeMemory();
    srand(time(NULL));
}

int center_eyewin(int window, int center)
{
    if ( !center ) {
        wd_pos(window, 1000, 1000);
    }

    // reposition and resize if necessary
    switch ( window ) {
        case WIND0:
            wd_siz(WIND0, gl_eyewinx, gl_eyewiny);
            if ( center ) {
                wd_pos(window, gl_fpx, gl_fpy);
            }
            break;
    }
    return 0;
}

/////////////////////////////////////////
// Functions specific to this paradigm //
/////////////////////////////////////////

int ShowFP()
{
    char buf[256];
    sprintf(buf, "ShowFP(%d, %d, %d, %d, %d, %d)",
        gl_fpx, gl_fpy, gl_fpsize, gl_fprgb[0], gl_fprgb[1], gl_fprgb[2]);
    udp_send_to(buf, "MAC");
    gl_fpon = 1;
    center_eyewin(WIND0, 1);
    return 0;
}

int HideFP()
{
    char buf[256];
    sprintf(buf, "HideFP()");
    udp_send_to(buf, "MAC");
    gl_fpon = 0;
    return 0;
}

int ShowStimLMScc()
{
    float cc;
    char buf[256];
    cc = gl_trialparams[TP_CONTRAST].val.fval;
    sprintf(buf, "ShowLMScc(%d, %d, %d, %d, %f, %f, %f)", gl_rfx, gl_rfy,
        gl_stimheight, gl_stimwidth,
        gl_ccmults[0]*cc, gl_ccmults[1]*cc, gl_ccmults[2]*cc);
    udp_send_to(buf, "MAC");
    return 0;
}

int HideStimLMScc()
{
    char buf[256];
    sprintf(buf, "HideLMScc()");
    udp_send_to(buf, "MAC");
    return 0;
}

int ErrorTrial()
{
    udp_send_to("AllOff()", "MAC");
    SendLaserIntensity(NULL);
    center_eyewin(WIND0, 0);
    gl_fpon = 0;
    return 0;
}

int InitCondArray()
{
    FreeMemory(); // before we change the number of conditions!
    BuildCondArray();
    ShuffleTrials();
    return 0;
}

/* set gl_ncontrasts = 0 for no stimulus */
/* set gl_ncontrasts = 1 for maximum contrast stimulus */
/* set gl_ncontrasts >= 2 for log-spaced contrasts between 0 and maximum, inclusive */
void BuildCondArray()
{
    int i, j;
    double *tmpptr;

    // change the following two lines to amend the available frequencies.
    int laserpowers[] = {0, gl_max_laser_power};
    int num_powers = 2;

    gl_condarray = malloc(sizeof(double *) * MAXNCONDS); // Upper bound on number of conditions allowed

    if ( !gl_max_laser_power ) {
        laserpowers[1] = NULL;
        num_powers = 1;
    }

    for ( i = 0; i < num_powers; ++i ) {
        for ( j = 0; j < gl_ncontrasts; ++j ) {
            tmpptr = (double *) malloc(sizeof(double) * NTRIALPARAMS);
            tmpptr[TP_LASERPOWER] = (double) laserpowers[i];
            tmpptr[TP_HEIGHT] = (double) gl_stimheight;
            tmpptr[TP_WIDTH] = (double) gl_stimwidth;
            tmpptr[TP_LCC] = (double) gl_ccmults[0];
            tmpptr[TP_MCC] = (double) gl_ccmults[1];
            tmpptr[TP_SCC] = (double) gl_ccmults[2];

            if ( j == 0 ) {
                tmpptr[TP_CONTRAST] = 0;
            } else if ( gl_ncontrasts < 3 ) {
                tmpptr[TP_CONTRAST] = gl_maxcontrast;
            } else {
                tmpptr[TP_CONTRAST] = pow(10, log10(gl_maxcontrast/(gl_ncontrasts-1))
                    + (j-1) * log10(gl_ncontrasts-1) / (gl_ncontrasts-2));
            }

            tmpptr[TP_TARGSHOWN] = tmpptr[TP_CONTRAST] > 0;
            gl_condarray[gl_nconds++] = tmpptr;
        }
    }
#ifdef TESTING
    dprintf("Number of conditions = %d\n", gl_nconds);
#endif
    if ( !gl_nconds ) gl_done = 1;

    return 0;
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
unsigned int my_randint(unsigned int upper_bound)
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

void ShuffleTrials() // Knuth shuffle
{
    int i, j, tmp;

    gl_condorder = malloc(sizeof(int) * gl_nconds);

    for ( i = 0; i < gl_nconds; ++i )
        gl_condorder[i] = i;

    for ( i = gl_nconds-1; i > 0; --i ) {
        j = my_randint(i+1);
        tmp = gl_condorder[j];
        gl_condorder[j] = gl_condorder[i];
        gl_condorder[i] = tmp;
    }
#ifdef TESTING
    dprintf("ShuffleTrials completed:");
    for ( i = 0; i < gl_nconds; dprintf(" %d", gl_condorder[i++]) );
    dprintf("\n");
#endif
}

int CorrectTrial()
{
    if ( ++gl_condidx == gl_nconds ) {
        if ( ++gl_nblockscompleted >= gl_ntrlspercond ) {
            gl_done = 1;
            dprintf("maximum number of trials reached\n");
        } else {
            ShuffleTrials();
        }
        gl_condidx = 0;
    }
    return 0;
}

int FreeMemory()
{
    int i;

    if ( gl_condarray ) {
        for ( i = 0; i < gl_nconds; ++i )
            free(gl_condarray[i]);

        free(gl_condarray);
        gl_condarray = NULL;
        free(gl_condorder);
        gl_condorder = NULL;
    }
    return 0;
}

#define TIMEBETWEENPULSES 2 // in ms
int SendLaserIntensity(int * laserpower)
{
    static int codeDrops = 0;
    int counter = 0, intensitybits = 0;
    long prev_t = i_b->i_time, curr_t;

    if ( laserpower )
        intensitybits = (*laserpower << 1) + 513;
    else
        intensitybits = (0 << 1) + 513;

    do {
        curr_t = i_b->i_time;
        if ( curr_t - prev_t >= TIMEBETWEENPULSES ) {
            if ( intensitybits & 1 ) // send a '1'
                dio_on(FTTL4);
            else // send a '0'
                dio_off(FTTL4);

            prev_t = curr_t;
            intensitybits >>= 1; // nibble LSB off
            ++counter;
        }
    } while ( counter < 10 );
    dio_off(FTTL4);

    if ( laserpower && (*laserpower) ) {
        ec_send_code_hi(STIMONCD);
    } else if ( !laserpower ) {
        ec_send_code_hi(STIMOFFCD);
    }
    return 0;
}
#undef TIMEBETWEENPULSES

// Declaration of statelist menu

VLIST state_vl[] = {{NS}};

VLIST fp_vl[] =
{
    {"X",             &gl_fpx, NP, NP, 0, ME_DEC},
    {"Y",             &gl_fpy, NP, NP, 0, ME_DEC},
    {"Size",          &gl_fpsize, NP, NP, 0, ME_DEC},
    {"Red",           &gl_fprgb[0], NP, NP, 0, ME_DEC},
    {"Green",         &gl_fprgb[1], NP, NP, 0, ME_DEC},
    {"Blue",          &gl_fprgb[2], NP, NP, 0, ME_DEC},
    {"Window width",  &gl_eyewinx, NP, NP, 0, ME_DEC},
    {"Window height", &gl_eyewiny, NP, NP, 0, ME_DEC},
    {NS},
};

VLIST optstim_vl[] =
{
    {"peak power", &gl_max_laser_power , NP, NP, 0, ME_DEC},
/*    {"nlasersteps", &gl_nlasersteps , NP, NP, 0, ME_DEC}, */
    {NS},
};

VLIST visstim_vl[] =
{
    {"Maximum contrast", &gl_maxcontrast, NP, NP, 0, ME_FLOAT},
    {"# of contrasts",   &gl_ncontrasts, NP, NP, 0, ME_DEC},
    {"Height",           &gl_stimheight, NP, NP, 0, ME_DEC},
    {"Width",            &gl_stimwidth, NP, NP, 0, ME_DEC},
    {"Lcc multiplier",   &gl_ccmults[0], NP, NP, 0, ME_FLOAT},
    {"Mcc multiplier",   &gl_ccmults[1], NP, NP, 0, ME_FLOAT},
    {"Scc multiplier",   &gl_ccmults[2], NP, NP, 0, ME_FLOAT},

    {NS},
};

VLIST misc_vl[] =
{
    {"Max # blocks", &gl_ntrlspercond, NP, NP, 0, ME_DEC},
    {NS},
};

MENU umenus[] = {
    {"Fixation",     &fp_vl, NP, NP, 0, NP, NP},
    {"Visual Stim",  &visstim_vl, NP, NP, 0, NP, NP},
    {"Optical Stim", &optstim_vl, NP, NP, 0, NP, NP},
    {"Misc.",        &misc_vl, NP, NP, 0, NP, NP},
    {NS},
};

char hm_sv_vl[] = "";

// User-supplied real-time variable table

RTVAR rtvars[] = {};

%%
id PARADIGMID
restart rinitf
main_set {
status ON
begin
    first:
        to closePort
    closePort:
        do udpclose()
        to openPort
    openPort:
        do udpopen()
        to endSlave
    endSlave:
        do EndSlavePrgm()
        time 50
        to startSlave
    startSlave:
        do StartSlavePrgm()
        time 50
        to initMac
    initMac:
        do SetupMac()
        time 1000
        to getRFxy
    getRFxy:
        do getRFxy()
        time 20
        to initCondArray
    initCondArray:
        do InitCondArray()
        to pause0
    pause0:
        to finish1 on 1 = gl_done
        to pause1 on +PSTOP & softswitch
        to initTrial
    pause1:
        code FPOFFCD
        do HideFP()
        to initTrial on -PSTOP & softswitch
    initTrial:
        do InitTrialParams()
        to plxStart
    plxStart:
        do dio_on(PLXREC)
        time 400
        to dropHeader on 1 = gl_firsttrial
        to fpOn on 0 = gl_fpon
        to fpacqwait
    dropHeader:
        do DropHeader()
        time 300
        to getGammaTable
    getGammaTable:
        do getGammaTable()
        time 3000
        to getMonSpd
    getMonSpd:
        do getMonSpd()
        time 1600
        to getFundamentals
    getFundamentals:
        do getFundamentals()
        time 1500
        to getBkgndrgb
    getBkgndrgb:
        do getBkgndrgb()
        time 200
        to fpOn on 0 = gl_fpon
        to fpacqwait
    fpOn:
        code FPONCD
        do ShowFP()
        to fpacqwait
    fpacqwait:
        time 2000
        to fixTime on -WD0_XY & eyeflag
        to error
    fixTime:
        code FPACQCD
        time 200
        to error on +WD0_XY & eyeflag
        to visstimon
    visstimon:
        code TARGONCD
        do ShowStimLMScc()
        time 30
        to error on +WD0_XY & eyeflag
        to laserStimOn
   laserStimOn:
        do SendLaserIntensity(&gl_laser_power)
        time 150
        to error on +WD0_XY & eyeflag
        to visstimoff
    visstimoff:
        code TARGOFFCD
        do HideStimLMScc()
        time 30
        to error on +WD0_XY & eyeflag
        to laserStimOff
    laserStimOff:
        do SendLaserIntensity(NULL)
        to error on +WD0_XY & eyeflag
        to postStimHold
    postStimHold: /* maintain fixation until the timer expires */
        time 200
        to error on +WD0_XY & eyeflag
        to rewon
    fpOff:
        code FPOFFCD
        do HideFP()
        to rewon
    rewon:
        do dio_on(REW)
        time 40
        to rewoff
    rewoff:
        do dio_off(REW)
        to correctTrial
    correctTrial:
        do CorrectTrial()
        to dropTrialParams
    dropTrialParams:
        do DropTrialParams()
        time 100
        to markTrialEnd
    markTrialEnd:
        code EOT
        time 10
        to plxStop
    plxStop:
        do dio_off(PLXREC)
        time 140
        to pause0
    error:
        code ABORTCD
        do ErrorTrial()
        time 200
        to dropTrialParams
    finish1:
        do FreeMemory()
        time 50
        to finish2
    finish2:
        do HideFP()
        time 50
        to pause0
    abort list:
}

msg_set {
status ON
begin
    msgFirst:
        to msgSecond
    msgSecond:
        do CheckMsg()
        to msgSecond
    abort list:
}
