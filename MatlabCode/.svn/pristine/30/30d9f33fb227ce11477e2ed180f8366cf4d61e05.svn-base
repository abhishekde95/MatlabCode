/*
 * Paradigm for optically and/or electrically stimulating V1 just after the fixation point goes out.
 * Can we evoke saccades into the RFs of the stimulated neurons?
 *
 * In the _interfreq version, there is only the Fix condition with varying stimulation frequencies.
 * There is no behavior required: we're looking for a pure electrophysiological effect.
 *
 * Original framework: GDLH 5/6/11
 * First production version: ZALB 6/7/11
 * Second production version: ZALB 10/26/11
 * Third production version: ZALB 02/08/12
 * _interfreq version: ZALB 01/29/13
 */

//#define TESTING

#include <stdlib.h>
#include <math.h>
#include "ldev.h"
#include "labcodes.h"
#include "UDPtransfer.h"
#include "rigconsts.h"
#include "FixStim.h"

#define MAXNCONDS 100
#define NTASKS 1 // # of saccade/stimulation tasks

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
int gl_ntrlspercond = 20;
int gl_laserstimon = 0;
int gl_stimdur = 500;
int gl_laseron_dur = 2; // on pulse dur
int gl_laseroff_dur = 0; // off pulse dur
int gl_actualfreq[] = {0, 0}; // whole part and remainder
int gl_force_duty_cycle = 1; // force .5 duty cycle

// Housekeeping globals
int gl_done = 0;
int gl_firsttrial = 1;
int gl_nconds = 0;
int gl_condidx = 0; // index into gl_condarray
int gl_nblockscompleted = 0; // this counter increments to gl_ntrlspercond
int *gl_condorder;
int **gl_condarray; // this restricts the trial parameters to be ints only
int gl_use_modulator_public = 1;
int gl_use_modulator = 1;
int gl_freq_idx = 0; // the frequency index -- overwritten by the trial params

/////////////////////////////////
// Stuff for TrialParams array //
/////////////////////////////////

// Indices for gl_trialparams and gl_condarray
#define TP_LASERSTIMFREQ 0
#define NTRIALPARAMS     1

int gl_trialparamcodes[NTRIALPARAMS] = {LASERSTIMFREQCD};
int gl_trialparamtypes[NTRIALPARAMS] = {INT};
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
    udp_close();  // Closing a closed port is better than opening an open one
    return 0;
}

///////////////////////////
//  Mac commands via UDP //
///////////////////////////

int SetupMac()
{
    char buf[256];
    sprintf(buf, "InitDisplay(%d, %d, %s)",
            GL_MONDIST, GL_SCREENWIDTH, GL_CALFILENAME);
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
    ec_send_code_tagged(NTRLSPERCONDCD, VALOFFSET+gl_ntrlspercond);
    ec_send_code_tagged(MODULATORCD, VALOFFSET+gl_use_modulator);

    gl_firsttrial = 0;
    return 0;
}

int InitTrialParams()
{
    int i;
    int *tmpptr;

    tmpptr = gl_condarray[gl_condorder[gl_condidx]];

    if (gl_use_modulator) {
        gl_freq_idx = tmpptr[TP_LASERSTIMFREQ];
    } else {
        SetStimTimes(tmpptr[TP_LASERSTIMFREQ]);
        tmpptr[TP_LASERSTIMFREQ] = gl_actualfreq[0]; // just the integer part
    }

    for (i = 0; i < NTRIALPARAMS; ++i) {
        gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
        gl_trialparams[i].casttype = gl_trialparamtypes[i];
        gl_trialparams[i].nelements = 1;
        // assign the trial params if they are known, otherwise set reasonable defaults
        switch ( gl_trialparamtypes[i] ) {
        case DOUBLE:
            gl_trialparams[i].val.dval = tmpptr[i];
            break;
        case FLOAT:
            gl_trialparams[i].val.fval = tmpptr[i];
            break;
        case INT:
            gl_trialparams[i].val.lval = tmpptr[i];
            break;
        default:
            dprintf("InitTrialParams: Unknown trial param [%d] cast type %d\n", i, gl_trialparams[i].casttype);
        }
    }
    
    SetStimDur(); /* need to add 40 ms to gl_stimdur to get laser state dwell time GDLH 5/19/15 */
    return 0;
}

int DropTrialParams()
{
    int i, ncodes = 0;
    char typebuffer[8];
    float mspercode = 0.5f; // Conservative; actual value is closer to 0.33

    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        if ( gl_trialparams[i].casttype == DOUBLE ) {
            sprintf(typebuffer, "double");
            ncodes += 9;
        } else if ( gl_trialparams[i].casttype == FLOAT ) {
            sprintf(typebuffer, "float");
            ncodes += 5;
        } else if ( gl_trialparams[i].casttype == INT ) {
            sprintf(typebuffer, "int");
            ncodes += 2;
        } else {
            dprintf("DropTrialParams: Unknown trial param [%d] cast type %d\n", i, gl_trialparams[i].casttype);
            continue;
        }
        sendToPlexDataStream(typebuffer, (void*) &gl_trialparams[i].val.lval,
            gl_trialparams[i].nelements, gl_trialparams[i].flagcodeval, 0);
    }
    set_times("dropTrialParams", my_ceil(ncodes * mspercode), 0);
    return 0;
}

/////////////////////////
//   REX setup stuff   //
/////////////////////////

void rinitf()
{
    PARADIGM_PRELUDE("FixStim_interfreq");

    udp_send_to("AllOff()", "MAC");

    // WIND0 is used for testing eye position
    wd_src_pos(WIND0, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_check(WIND0, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);
    wd_cntrl(WIND0, WD_ON); // don't use WD_OFF, just move the eye window around

    center_eyewin(WIND0, 1); // size and position the eye window

    dio_off(FTTLBYTE1); // turn off TTLs 1A - 8A (the first byte)

    // default global variables
    gl_firsttrial = 1;
    gl_nconds = 0;
    gl_condidx = 0;
    gl_nblockscompleted = 0;
    gl_laserstimon = 0;
    gl_done = 0;

    gl_use_modulator = !!gl_use_modulator_public;

    set_times("laserStimOn", gl_stimdur, 0);
    set_times("stimThink", gl_stimdur, 0);
    set_times("setLaserModParams", gl_stimdur, 0);

    i_b->t_wrate = 1;

    FreeMemory();
    srand(time(NULL));
}

int center_eyewin(int window, int center)
{
    if (!center) {
        wd_pos(window, 1000, 1000);
    }

    // reposition and resize if necessary
    switch (window) {
    case WIND0:
        wd_siz(WIND0, gl_eyewinx, gl_eyewiny);
        if (center) {
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
    center_eyewin(WIND0, 1);
    return 0;
}

int HideFP()
{
    char buf[256];
    sprintf(buf, "HideFP()");
    udp_send_to(buf, "MAC");
    return 0;
}

int ErrorTrial()
{
    gl_laserstimon = 0;
    udp_send_to("AllOff()", "MAC");
    center_eyewin(WIND0, 0);
    if (gl_use_modulator) dio_off(FTTL3); // make sure the GATE is off
    return 0;
}

int SetStimState(int * stim, int val)
{
    *stim = val;
    return 0;
}

int SetStimDur(int flag, MENU *mp, char *astr, VLIST *vlp, int *tvadd)
{
    set_times("laserStimOn", gl_stimdur, 0);
    set_times("stimThink", gl_stimdur+40, 0); /* We wait in stimthink on 0 frequency trials */
    set_times("setLaserModParams", gl_stimdur+40, 0); /* We wait in setLaserModParams on > 0 frequency trials */
     /* +40 to account for the 40 ms required to send 20 bits at 2 ms/bit in setLaserModParams */
    return 0;
}

void SetStimTimes(int laserfreq)
{
    float freqtemp;

    if ( laserfreq > 0 ) {
        if (gl_force_duty_cycle) {
            gl_laseron_dur = gl_laseroff_dur = (int) (500.0f / laserfreq);
        } else {
            gl_laseroff_dur = (int) (1000.0f / laserfreq - gl_laseron_dur);
        }

        // the laser pulse cannot be on (or off) for less than 2 ms
        gl_laseron_dur = gl_laseron_dur < 2 ? 2 : gl_laseron_dur;
        gl_laseroff_dur = gl_laseroff_dur < 2 ? 2 : gl_laseroff_dur;

        // set the states' duration -- subtract 1 ms from each state to account for REX latency
        // the state times of 0 and 1 process the state in 1 ms (i.e., they're equivalent)
        set_times("laserPulseDurOn", gl_laseron_dur-1, 0);
        set_times("laserPulseDurOff", gl_laseroff_dur-1, 0);

        // here we populate some real-time variables with the actual frequency
        freqtemp = 1000.0f / (gl_laseron_dur + gl_laseroff_dur);
        gl_actualfreq[0] = (int) freqtemp;
        gl_actualfreq[1] = 100 * (freqtemp - gl_actualfreq[0]);
    } else {
        gl_actualfreq[0] = 0;
        gl_actualfreq[1] = 0;
    }
}

int InitCondArray()
{
    FreeMemory();  // before we change the number of conditions!
    BuildCondArray();
    ShuffleTrials();
    return 0;
}

void BuildCondArray()
{
    int i, j, *tmpptr, num_freqs, frequencies[32];

    // change the following two lines to amend the available frequencies.
    // INCLUDE 0 if you want a no stimulation trial!!
	frequencies[0] = 0;
    if (gl_use_modulator) {		
		num_freqs = 10;
		for (i = 1; i < num_freqs-1; ++i)
			frequencies[i] = 1 << i; // powers of 2
		frequencies[num_freqs-1] = 500;
    } else {
		frequencies[1] = 50;
		frequencies[2] = 100;
		frequencies[3] = 125;
		frequencies[4] = 250;
        num_freqs = 5;
    }

    // with the modulator, the sinusoid has frequency = .5*((1000^(1/255))^k), where k is the frequency index
    // command to modulator should be k = 255*ln(f/Fmin)/ln(Fmax/Fmin)
    // convert Hz to frequency index
    if (gl_use_modulator) {
        for (i = 0; i < num_freqs; i++) {
            if (frequencies[i] != 0) {
                frequencies[i] =  my_round(255*logf((float)frequencies[i]/0.5)/logf(500.0/0.5)); // 500 and 0.5 are max and min freq in Hz, respectively
            }
        }
    }

    gl_condarray = malloc(sizeof(int*) * MAXNCONDS); // Upper bound on number of conditions allowed

    for ( i = 0; i < NTASKS; ++i ) {
        for ( j = 0; j < num_freqs; ++j ) {
            tmpptr = malloc(sizeof(int) * NTRIALPARAMS);
            tmpptr[TP_LASERSTIMFREQ] = frequencies[j]; // later on these are updated when _not_ using the modulator
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
    memset(gl_condorder, 0, sizeof(int)*gl_nconds);

    for (i = 0; i < gl_nconds; ++i)
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

unsigned int bit_parity(unsigned int v)
{
    int p = 0;

    while (v) {
        p = !p;
        v = v & (v - 1);
    }
    return p;
}

#define TIMEBETWEENPULSES 2
int SendModulatorParams(int *pk)
{
    long prev_t = i_b->i_time, curr_t;
    unsigned int byte1, byte2, p;
    unsigned int word1, word2, message, counter = 0;

    /*
      put laser into sinusoidal power mode

      second 10-bit word:
        word2 = k << 1 + 1, where
          +1 is the start bit
          k frequency index

      calculate the bit parity of k, call it p

      first byte (with checksum bit (position 2) is zero):
        0x98 = 0b10011000 has a "3" in bits 6,5,4,3, and bit 7 is always set
        this always has odd parity so we flip the parity p and stick it in bit 2
        0x98 | !p << 2, which ensures the whole 20-bit message has even parity

      first 10-bit word:
        word1 = (0x98 | (!p << 2)) << 1 + 1

      whole 20-bit message:
        word2 << 10 + word1
        shift over the high 10-bit sequence to the higher bits
    */
    byte1 = 0x98; // b6,5,4,3 = 3 (base 10) and b7 = 1
    byte2 = *pk;
    p = bit_parity(byte2);
    word1 = ((byte1 | (!p << 2)) << 1) + 1;
    word2 = (byte2 << 1) + 1;
    message = (word2 << 10) + word1;

    do {
        curr_t = i_b->i_time;
        if ( curr_t - prev_t >= TIMEBETWEENPULSES ) {
            if ( message & 1 )
                dio_on(FTTL4);
            else // send a '0'
                dio_off(FTTL4);

            prev_t = curr_t;
            message >>= 1; // nibble LSB off
            ++counter;
        }
    } while ( counter < 20 );
    dio_off(FTTL4);
    dio_on(FTTL3); // turn on the gate -- laser should be on now
    return STIMONCD; // need to test the timing of this code rel to the laser
}
#undef TIMEBETWEENPULSES

int turnOffModulator()
{
    dio_off(FTTL3); // turn off the "gate" input turns off the laser
    return STIMOFFCD;
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
        gl_condarray = NULL;
    }
    return 0;
}

int UseModulator()
{
    return gl_use_modulator && gl_freq_idx > 0;
}

// Declaration of statelist menu

VLIST state_vl[] = {{NS}};

VLIST fp_vl[] =
{
    {"X",             &gl_fpx,      NP, NP, 0, ME_DEC},
    {"Y",             &gl_fpy,      NP, NP, 0, ME_DEC},
    {"Size",          &gl_fpsize,   NP, NP, 0, ME_DEC},
    {"Red",           &gl_fprgb[0], NP, NP, 0, ME_DEC},
    {"Green",         &gl_fprgb[1], NP, NP, 0, ME_DEC},
    {"Blue",          &gl_fprgb[2], NP, NP, 0, ME_DEC},
    {"Window width",  &gl_eyewinx,  NP, NP, 0, ME_DEC},
    {"Window height", &gl_eyewiny,  NP, NP, 0, ME_DEC},
    {NS},
};

VLIST stim_vl[] =
{
    {"Use modulator", &gl_use_modulator_public, NP, NP, 0, ME_DEC},
    {".5 duty cycle override", &gl_force_duty_cycle, NP, NP, 0, ME_DEC},
    {"Stim ON duration", &gl_laseron_dur, NP, NP, 0, ME_DEC},
    {"Stim duration",    &gl_stimdur,     NP, SetStimDur, ME_AFT, ME_DEC},
    {NS},
};

VLIST misc_vl[] =
{
    {"Max # blocks", &gl_ntrlspercond, NP, NP, 0, ME_DEC},
    {NS},
};

MENU umenus[] = {
    {"Fixation", &fp_vl,   NP, NP, 0, NP, NP},
    {"Stim",     &stim_vl, NP, NP, 0, NP, NP},
    {"Misc.",    &misc_vl, NP, NP, 0, NP, NP},
    {NS},
};

char hm_sv_vl[] = "";

// User-supplied real-time variable table

RTVAR rtvars[] = {
    {"Freq index",       &gl_freq_idx},
    {"ON duration",      &gl_laseron_dur},
    {"OFF duration",     &gl_laseroff_dur},
    {"Freq (int part)",  &gl_actualfreq[0]},
    {"Freq (remainder)", &gl_actualfreq[1]},
    {"", 0},
};

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
        to initCondArray
    initCondArray:
        do InitCondArray()
        to pause0
    pause0:
        to finished on 1 = gl_done
        to pause1 on +PSTOP & softswitch
        to initTrial
    pause1:
        to initTrial on -PSTOP & softswitch
    initTrial:
        do InitTrialParams()
        time 200
        to plxStart
    plxStart:
        do dio_on(PLXREC)
        time 400
        to dropHeader on 1 = gl_firsttrial
        to fpOn
    dropHeader:
        do DropHeader()
        time 300
        to fpOn
    fpOn:
        code FPONCD
        do ShowFP()
        time 2000
        to fixTime on -WD0_XY & eyeflag
        to error
    fixTime:
        code FPACQCD
        time 249
        to error on +WD0_XY & eyeflag
        to finished on 1 = gl_done
        to stimThink
    stimThink: /* goto stim state if frequency != 0, otherwise wait the stim duration (gl_stimdur) */
        to setLaserModParams on 1 % UseModulator
        to laserStimOn on 0 < gl_actualfreq[0]
        to error on +WD0_XY & eyeflag
        to postStimHold
    setLaserModParams: /* time of this state (laser stimulation duration) is via gl_stimdur */
        do SendModulatorParams(&gl_freq_idx)
        to error on +WD0_XY & eyeflag
        to turnOffModulator
    turnOffModulator:
        do turnOffModulator()
        to postStimHold
    laserStimOn:
        code STIMONCD
        /*time 500*/ /* this is the laser stimulation duration (now controled by user-defined variable) */
        do SetStimState(&gl_laserstimon,1)
        to error on +WD0_XY & eyeflag
        to laserStimOff
    laserStimOff:
        code STIMOFFCD
        do SetStimState(&gl_laserstimon,0)
        to postStimHold
    postStimHold: /* maintain fixation until the timer expires */
        time 249
        to error on +WD0_XY & eyeflag
        to fpOff
    fpOff:
        code FPOFFCD
        do HideFP()
        to rewOn
    rewOn:
        code REWCD
        do dio_on(REW)
        time 60 /* this is the reward duration */
        to rewOff
    rewOff:
        do dio_off(REW)
        to correctTrial
    correctTrial:
        time 500
        do CorrectTrial()
        to dropTrialParams
    dropTrialParams:
        do DropTrialParams()
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
    finished:
        do FreeMemory()
        time 100
        to pause0
    abort list:
}

laserstim_set {
status ON
begin
    stimZero:
        to stimFirst
    stimFirst:
        do dio_off(FTTL2)
        to laserOnPhase on 1 = gl_laserstimon
    laserOnPhase:
        do dio_on(FTTL2)
        to laserPulseDurOn
    laserPulseDurOn:
        to stimFirst on 0 = gl_laserstimon
        to laserOffPhase
    laserOffPhase:
        do dio_off(FTTL2)
        to laserPulseDurOff
    laserPulseDurOff:
        to stimFirst on 0 = gl_laserstimon
        to laserOnPhase
    abort list:
}
