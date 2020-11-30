/*
 * Paradigm for optically and/or electrically stimulating V1 just after the fixation point goes out.
 * Can we evoke saccades into the RFs of the stimulated neurons?
 *
 * Original framework: GDLH 5/6/11
 * First production version: ZALB 6/7/11
 * Second production version: ZALB 10/26/11
 * Implementing the TTL laser modulator: ZALB 04/13/15
 */

// #define TESTING

#include <stdlib.h>
#include <math.h>
#include "ldev.h"
#include "labcodes.h"
#include "UDPtransfer.h"
#include "rigconsts.h"
#include "FixStim.h"

#define MAXNCONDS 100
#define NTASKS 6 // # of saccade/stimulation tasks

//////////////////////////////////
// Global variable declarations //
//////////////////////////////////

// Parameters that the user can change
int gl_fpx = 0;
int gl_fpy = 0;
#ifdef TESTING
int gl_eyewinx = 12;
int gl_eyewiny = 12;
int gl_address_bits_override = -1;
int gl_second_byte_override = -1;
#else
int gl_eyewinx = 7;
int gl_eyewiny = 7;
#endif
int gl_targwinx = 14;
int gl_targwiny = 14;
int gl_fpsize = 2;
int gl_targsize = 4;
int gl_fprgb[] = {0, 0, 0};
int gl_targrgb[] = {0, 0, 0};
int gl_ntrlspercond = 20;
int gl_ttl_stimfreq = 1;
int gl_laser_power = 0;
float gl_sine_freq = 50.0f;
int gl_use_modulator_public = 1;
int gl_modulator_dc = 0;

#define WINDOPT 4
#define WINDTARG 6

// Stimulus parameters obtained via requests to Slave
int gl_rfx = 0;
int gl_rfy = 0;
int gl_optx = 0;
int gl_opty = 0;

// Housekeeping globals
int gl_done = 0;
int gl_firsttrial = 1;
int gl_curtasktype = -1;
int gl_curstimtype = 0;
int gl_ttlstimon = 0;
int gl_nconds = 0;
int gl_ntasks = NTASKS;
int gl_taskset[NTASKS] = {0, 0, 0, 1, 0, 0}; // all, bx, bx+ttl, bx+opt, fix+ttl, fix+opt
int gl_condidx = 0; // index into gl_condarray
int gl_nblockscompleted = 0; // this counter increments to gl_ntrlspercond
int gl_use_modulator = 0;
int gl_freq_idx = 170; // for the laser modulator; overwritten by user-specified frequency
int *gl_condorder = NULL;
int **gl_condarray = NULL;

// store certain event times
long gl_timer_begin, gl_timer_end, gl_fpoff_t, gl_targon_t;

/////////////////////////////////
//   Specifying the task sets  //
/////////////////////////////////

#define TASKS_FIXTTL (1<<0 | 1<<1)
#define TASKS_FIXOPT (1<<0 | 1<<2)
#define TASKS_BX     (1<<0 | 1<<3)
#define TASKS_BXTTL  TASKS_BX | (1<<1 | 1<<4)
#define TASKS_BXOPT  TASKS_BX | (1<<2 | 1<<5)
#define TASKS_ALL    TASKS_BXTTL | TASKS_BXOPT
#define NTASKSETS    6

int gl_task_masks[] = {TASKS_ALL, TASKS_BX, TASKS_BXTTL, TASKS_BXOPT, TASKS_FIXTTL, TASKS_FIXOPT};
int gl_ntasksperset[] = {6, 2, 4, 4, 2, 2};

/////////////////////////////////
// Stuff for TrialParams array //
/////////////////////////////////

// Indices for gl_trialparams and gl_condarray
#define TP_TARGSHOWN   0
#define TP_STIMTYPE    1
#define TP_STIMFREQ    2
#define TP_LASERPOWER  3
#define NTRIALPARAMS   4

int gl_trialparamcodes[NTRIALPARAMS] = {TARGSHOWNCD, STIMTYPECD, STIMFREQCD, LASERPOWERCD};
int gl_trialparamtypes[NTRIALPARAMS] = {INT, INT, INT, INT};
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

////////////////////////////////////////////
// Functions for sending requests via UDP //
////////////////////////////////////////////

const char RFXY_REQ[] = "sendToRex(udpCom, gl.bar.xy, 'double', message);";

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

//////////////////////////
// Getting UDP commands //
//////////////////////////

int CheckMsg()
{
    char inMsg[REXBUFF];
    static message_struct message;
    int messageAvailable;
    double dvals[2];

    inMsg[0] = '\0';
    messageAvailable = udp_check(0);

    if (messageAvailable)
    {
        udp_read(inMsg, REXBUFF);
        prepareMsg(&message, inMsg);

        if ( !strcmp(message.header, RFXY_REQ) ) {
            hex2double(message.contents, (void*) &dvals);
            if ( message.size && (dvals[0] || dvals[1]) ) {
                gl_rfx = my_round(10 * dvals[0]);
                gl_rfy = my_round(10 * dvals[1]);
                gl_optx = gl_rfx;
                gl_opty = gl_rfy;
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
    ec_send_code_tagged(TARGSIZECD, VALOFFSET+gl_targsize);
    ec_send_code_tagged(TARGRCD, VALOFFSET+gl_targrgb[0]);
    ec_send_code_tagged(TARGGCD, VALOFFSET+gl_targrgb[1]);
    ec_send_code_tagged(TARGBCD, VALOFFSET+gl_targrgb[2]);
    ec_send_code_tagged(TARGWINXCD, VALOFFSET+gl_targwinx);
    ec_send_code_tagged(TARGWINYCD, VALOFFSET+gl_targwiny);
    ec_send_code_tagged(RFXCD, VALOFFSET+gl_rfx);
    ec_send_code_tagged(RFYCD, VALOFFSET+gl_rfy);
    ec_send_code_tagged(NTRLSPERCONDCD, VALOFFSET+gl_ntrlspercond);

    gl_firsttrial = 0;
    if (!gl_use_modulator) dio_on(FTTL3); // turn on the laser trigger
    return 0;
}

int InitTrialParams()
{
    int i;
    int *tmpptr;

    gl_curtasktype = gl_condorder[gl_condidx]; // for use in the state set
    tmpptr = gl_condarray[gl_curtasktype];
    gl_curstimtype = tmpptr[TP_STIMTYPE]; // for use in the state set

    if (gl_modulator_dc && gl_use_modulator)
        tmpptr[TP_STIMFREQ] = 999;

    for (i = 0; i < NTRIALPARAMS; ++i) {
        gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
        gl_trialparams[i].casttype = gl_trialparamtypes[i];
        gl_trialparams[i].nelements = 1;
        switch (gl_trialparamtypes[i]) {
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
    PARADIGM_PRELUDE("FixStim");

    udp_send_to("AllOff()", "MAC");

    // WIND0 and WIND1 are used for testing eye position
    wd_src_pos(WIND0, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_pos(WIND1, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_check(WIND0, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);
    wd_src_check(WIND1, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);
    wd_pos(WIND0, gl_fpx, gl_fpy);
    wd_pos(WIND1, gl_rfx, gl_rfy);
    wd_siz(WIND0, gl_eyewinx, gl_eyewiny);
    wd_siz(WIND1, gl_targwinx, gl_targwiny);
    wd_cntrl(WIND0, WD_ON);
    wd_cntrl(WIND1, WD_ON);

    // WINDOPT and WINDTARG are for displaying when stimuli occur
    wd_src_pos(WINDOPT, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_pos(WINDTARG, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_pos(WINDOPT, 1000, 1000); // off screen
    wd_pos(WINDTARG, 1000, 1000); // off screen
    wd_siz(WINDOPT, gl_eyewinx, gl_eyewiny);
    wd_siz(WINDTARG, gl_targwinx, gl_targwiny);
    wd_cntrl(WINDTARG, WD_ON);
    wd_cntrl(WINDOPT, WD_ON);

    dio_off(FTTLBYTE1); // turn off TTLs 1A - 8A (the first byte)

    // default global variables
    gl_firsttrial = 1;
    gl_nconds = 0;
    gl_condidx = 0;
    gl_nblockscompleted = 0;
    gl_ttlstimon = 0;
    gl_curtasktype = -1;
    gl_curstimtype = 0;
    gl_done = 0;
    gl_ntasks = 0;

    i_b->t_wrate = 1;

    gl_timer_begin = 0;
    gl_timer_end = 0;
    gl_fpoff_t = 0;
    gl_targon_t = 0;

    gl_use_modulator = !!gl_use_modulator_public;

    SetStimTimes();

    FreeMemory();
    srand(time(NULL));
}

int update_eyewin(int window)
{
    switch (window) {
        case WIND0:
            wd_pos(WIND0, gl_fpx, gl_fpy);
            wd_siz(WIND0, gl_eyewinx, gl_eyewiny);
        case WIND1:
            wd_pos(WIND1, gl_rfx, gl_rfy);
            wd_siz(WIND1, gl_targwinx, gl_targwiny);
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
    update_eyewin(WIND0);
    return 0;
}

int HideFP()
{
    char buf[256];
    sprintf(buf, "HideFP()");
    udp_send_to(buf, "MAC");
    get_event_time(&gl_fpoff_t);
    return 0;
}

int ShowTarg()
{
    char buf[256];
    sprintf(buf, "ShowTarg(%d, %d, %d, %d, %d, %d)",
                gl_rfx, gl_rfy, gl_targsize, gl_targrgb[0], gl_targrgb[1], gl_targrgb[2]);
    udp_send_to(buf, "MAC");
    get_event_time(&gl_targon_t);
    update_eyewin(WIND1);

    wd_pos(WINDTARG, gl_rfx, gl_rfy); // move box on screen

    return 0;
}

int HideTarg()
{
    char buf[256];
    sprintf(buf, "HideTarg()");
    udp_send_to(buf, "MAC");
    wd_pos(WIND1, 1000, 1000); // off screen

    wd_pos(WINDTARG, 1000, 1000); // move box off screen

    return 0;
}

int ErrorTrial()
{
    gl_timer_begin = 0; gl_timer_end = 0;
    gl_fpoff_t = 0; gl_targon_t = 0;

    gl_ttlstimon = 0;

    if (gl_use_modulator) dio_off(FTTL3); // make sure the GATE is off

    udp_send_to("AllOff()", "MAC");

    wd_pos(WIND1, 1000, 1000); // off screen

    wd_pos(WINDOPT, 1000, 1000); // move box off screen
    wd_pos(WINDTARG, 1000, 1000); // move box off screen

    return 0;
}

int SetStimState(int * stim, int val)
{
    *stim = val;
    if ( *stim )
        wd_pos(WINDOPT, gl_optx, gl_opty);
    else
        wd_pos(WINDOPT, 1000, 1000);
    return 0;
}

void get_event_time(long * val)
{
    *val = i_b->i_time;
}

// credit to Shin!
static int timerCheck()
{
    return (i_b->i_time < gl_timer_end) ? 0 : 1;
}

// type- 0: A fixed delay, 1: A delay from a time point
#define DTIMER_XSTATES 0L /* test this for accuracy */
int timerSet(int type, long delay, long * from)
{
    if ( type == 1 && from )
        gl_timer_begin = *from;
    else
        gl_timer_begin = i_b->i_time;

    gl_timer_end = gl_timer_begin - DTIMER_XSTATES + delay;

    return 0;
}

void SetStimTimes()
{
    // assuming a .5 duty cycle
    float ttltemp = 500.0f / gl_ttl_stimfreq;
    // set the states' duration -- subtract 1 ms from each state to account for REX latency
    set_times("ttlPulseDurOn", ((int) ttltemp)-1, 0);
    set_times("ttlPulseDurOff", my_ceil(ttltemp)-1, 0);
}

int InitCondArray()
{
    FreeMemory(); // before we change the number of conditions!
    BuildCondArray();
    setupConditionOrder();
    ShuffleTrials();
    return 0;
}

void BuildCondArray()
{
    int i, *tmpptr;

    gl_condarray = malloc(sizeof(int*) * MAXNCONDS); // Upper bound on number of conditions allowed

    for ( i = 0; i < NTASKS; ++i ) {
        tmpptr = malloc(sizeof(int) * NTRIALPARAMS);
        tmpptr[TP_TARGSHOWN] = 0;
        tmpptr[TP_STIMTYPE] = 0;
        tmpptr[TP_STIMFREQ] = 0;
        tmpptr[TP_LASERPOWER] = 0;
        switch (i) {
            case 0: // blank
                break;
            case 1: // no target, ttl stim
                tmpptr[TP_STIMTYPE] = 1;
                tmpptr[TP_STIMFREQ] = gl_ttl_stimfreq;
                break;
            case 2: // no target, opt stim
                tmpptr[TP_STIMFREQ] = gl_use_modulator ? gl_freq_idx : 0;
                tmpptr[TP_STIMTYPE] = 2+gl_use_modulator;
                tmpptr[TP_LASERPOWER] = gl_use_modulator ? -1 : gl_laser_power;
                break;
            case 3: // target, no stim
                tmpptr[TP_TARGSHOWN] = 1;
                break;
            case 4: // target, ttl stim
                tmpptr[TP_TARGSHOWN] = 1;
                tmpptr[TP_STIMTYPE] = 1;
                tmpptr[TP_STIMFREQ] = gl_ttl_stimfreq;
                break;
            case 5: // target, opt stim
                tmpptr[TP_STIMFREQ] = gl_use_modulator ? gl_freq_idx : 0;
                tmpptr[TP_TARGSHOWN] = 1;
                tmpptr[TP_STIMTYPE] = 2+gl_use_modulator;
                tmpptr[TP_LASERPOWER] = gl_use_modulator ? -1 : gl_laser_power;
                break;
        }
        gl_condarray[gl_nconds++] = tmpptr; // index then increment
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

void setupConditionOrder()
{
    int i, taskmask, counter;

    for ( i = 0; i < NTASKSETS; ++i ) {
        if ( gl_taskset[i] ) {
            taskmask = gl_task_masks[i]; // get the bit mask for the task set
            gl_ntasks = gl_ntasksperset[i]; // how many tasks in this set?
            break;
        }
    }

    gl_condorder = malloc(sizeof(int) * gl_ntasks);
    memset(gl_condorder, 0, sizeof(int)*gl_ntasks);

    counter = 0;
    for (i = 0; i < NTASKS; ++i)
        if ( (taskmask & (1 << i)) != 0 ) // is task i part of the task set?
            gl_condorder[counter++] = i;
}

void ShuffleTrials()
{
    int i, j, tmp;

    for ( i = gl_ntasks-1; i > 0; --i ) {
        j = my_randint(i+1);
        tmp = gl_condorder[j];
        gl_condorder[j] = gl_condorder[i];
        gl_condorder[i] = tmp;
    }
#ifdef TESTING
    dprintf("ShuffleTrials completed:");
    for ( i = 0; i < gl_ntasks; dprintf(" %d", gl_condorder[i++]) );
    dprintf("\n");
#endif
}

int CorrectTrial()
{
    int *tp = gl_condarray[gl_curtasktype];

    gl_timer_begin = 0; gl_timer_end = 0;
    gl_fpoff_t = 0; gl_targon_t = 0;

    if ( gl_curtasktype == 1 )
        dprintf("### TTL stim (no target) trial\n");
    else if ( gl_curtasktype == 2 )
        dprintf("*** Optical stim (no target) trial\n");
#ifdef TESTING
    dprintf("block: %d, trial: #%d, task: #%d, params: (%d,%d,%d,%d)\n",
        gl_nblockscompleted, gl_condidx, gl_curtasktype, tp[TP_TARGSHOWN],
        tp[TP_STIMTYPE], tp[TP_STIMFREQ], tp[TP_LASERPOWER]);
#endif

    if ( ++gl_condidx == gl_ntasks ) {
        if ( ++gl_nblockscompleted >= gl_ntrlspercond ) {
            gl_done = 1;
            dprintf("maximum number of trials reached\n");
        }
        else
            ShuffleTrials();
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
int SendLaserIntensity(int* laserpower)
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
            if ( intensitybits & 1 )// send a '1'
                dio_on(FTTL4);
            else // send a '0'
                dio_off(FTTL4);

            prev_t = curr_t;
            intensitybits >>= 1; // nibble LSB off
            ++counter;
        }
    } while ( counter < 10 );
    dio_off(FTTL4);

    if ( codeDrops++ % 2 == 0 ) { // this is 3 ms before the laser turns on (CL473-050-O Crystal Laser) ZALB, GDLH 02/15/2013
        ec_send_code_hi(STIMONCD);
        wd_pos(WINDOPT, gl_optx, gl_opty); // move box on screen
    } else {
        ec_send_code_hi(STIMOFFCD);
        wd_pos(WINDOPT, 1000, 1000); // move box off screen
    }

    return 0;
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

int SendModulatorParams(int *pk)
{
    long prev_t = i_b->i_time, curr_t;
    unsigned int byte1, byte2, p;
    unsigned int word1, word2, message, counter = 0;

    /*
      Control the 15CE10 laser modulator:

      Construct two bytes and send them over TTL to change the various operating
      modes of the modulator.

      First byte with postions (MSB) 7 6 5 4 3 2 1 0 (LSB):
        bits 0, 1 are always zero
        bit 2 is a checksum flag maintains the even parity of the whole message
        bits 3-6 hold the "address" mode which is an integer 0..15
          * address == 2: continuous mode
          * address == 3: sinusoidal mode
          * this code doesn't implement any other mode
        bit 7 is always 1

      Second byte is a frequency index k = [0,255] for the sinusoidal mode. The
      freqency is related to the index k via f ~ .5*1.02746^k.

      Each byte is flanked by a start bit 1 (at the lower end) and a stop bit 0
      (at the higher end). In total, we send over two bytes, two start bits,
      and two stop bits for a total of a 20-bit message.

      The code below constructs byte1 and byte2, calculates the parity of this
      two-byte sequence, updates the checksum bit in byte1, tacks on the start
      bits to each 10-bit word (silently the stop bits are present), and sends
      the message over TTL.
    */

    if (gl_modulator_dc) { // DC mode
        // address mode 2: "continuous" laser while GATE is high
        byte1 = (2 << 3) + (1 << 7);
        byte2 = 0;
    } else {
        // address mode 3: sinusoidal amplitude
        byte1 = (3 << 3) + (1 << 7);
        // sinusoidal amplitude at frequency index k == *pk
        byte2 = *pk;
    }
#ifdef TESTING
    if (gl_address_bits_override != -1 && gl_second_byte_override != -1) {
        byte1 = (gl_address_bits_override << 3) + 128;
        byte2 = gl_second_byte_override;
    }
#endif
    p = bit_parity((byte2 << 8) + byte1); // parity of the bytes
    byte1 |= (p << 2); // put in checksum bit
    word1 = (byte1 << 1) + 1; // start bit
    word2 = (byte2 << 1) + 1; // start bit
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
    wd_pos(WINDOPT, gl_optx, gl_opty); // move box on screen
    dio_on(FTTL3); // turn on the gate -- laser should be on now
    return STIMONCD; // need to test the timing of this code rel to the laser
}
#undef TIMEBETWEENPULSES

int turnOffModulator()
{
    dio_off(FTTL3); // turn off the "gate" input turns off the laser
    wd_pos(WINDOPT, 1000, 1000);
    return STIMOFFCD;
}

int ProcessSineFreq(int flag, MENU *mp, char *astr, VLIST *vlp, int *tvadd)
{
    if ( gl_sine_freq < 0.5f || gl_sine_freq > 500.0f )
        gl_sine_freq = 50.0f;
#ifdef TESTING
    if (gl_address_bits_override != -1) {
        // don't let the user set some weird mode
        if (gl_address_bits_override > 9 || (gl_address_bits_override < 9 && gl_address_bits_override > 4) || gl_address_bits_override < 0)
            gl_address_bits_override = -1;
    }
    if (gl_second_byte_override > 255 || gl_second_byte_override < 0) gl_second_byte_override = -1;
#endif
    gl_freq_idx = my_round(255.0*log(1.0*gl_sine_freq/0.5)/log(500.0/0.5));
    if ( gl_use_modulator && gl_condarray ) {
        gl_condarray[2][TP_STIMFREQ] = gl_freq_idx;
        gl_condarray[5][TP_STIMFREQ] = gl_freq_idx;
    }
    return 0;
}

int ProcessLaserPower(int flag, MENU *mp, char *astr, VLIST *vlp, int *tvadd)
{
    if ( gl_laser_power < 0 || gl_laser_power > 255 )
        gl_laser_power = 0;
    if ( gl_condarray ) {
        gl_condarray[2][TP_LASERPOWER] = gl_laser_power;
        gl_condarray[5][TP_LASERPOWER] = gl_laser_power;
    }
    return 0;
}

int StimFreq(int flag, MENU *mp, char *astr, VLIST *vlp, int *tvadd)
{
    SetStimTimes();
    // change the condition array so the correct frequencies drop into the file
    if (gl_condarray) {
        gl_condarray[1][TP_STIMFREQ] = gl_ttl_stimfreq;
        gl_condarray[4][TP_STIMFREQ] = gl_ttl_stimfreq;
    }
    return 0;
}

int ClearTaskSet(int flag, MENU *mp, char *astr, VLIST *vlp, int *tvadd)
{
    int i;
    for ( i = 0; i < NTASKSETS; ++i ) gl_taskset[i] = 0;
    return 0;
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

VLIST targ_vl[] =
{
    {"X",             &gl_rfx,        NP, NP, 0, ME_DEC},
    {"Y",             &gl_rfy,        NP, NP, 0, ME_DEC},
    {"Opt X",         &gl_optx,       NP, NP, 0, ME_DEC},
    {"Opt Y",         &gl_opty,       NP, NP, 0, ME_DEC},
    {"Size",          &gl_targsize,   NP, NP, 0, ME_DEC},
    {"Red",           &gl_targrgb[0], NP, NP, 0, ME_DEC},
    {"Green",         &gl_targrgb[1], NP, NP, 0, ME_DEC},
    {"Blue",          &gl_targrgb[2], NP, NP, 0, ME_DEC},
    {"Window width",  &gl_targwinx,   NP, NP, 0, ME_DEC},
    {"Window height", &gl_targwiny,   NP, NP, 0, ME_DEC},
    {NS},
};

VLIST stim_vl[] =
{
    {"Use modulator", &gl_use_modulator_public, NP, NP, 0, ME_DEC},
    {"Modulator DC mode", &gl_modulator_dc, NP, NP, 0, ME_DEC},
    {"Laser power",    &gl_laser_power, NP, ProcessLaserPower, ME_AFT, ME_DEC},
    {"TTL frequency", &gl_ttl_stimfreq, NP, StimFreq,          ME_AFT, ME_DEC},
    {"Sine frequency",   &gl_sine_freq, NP, ProcessSineFreq,   ME_AFT, ME_FLOAT},
#ifdef TESTING
    {"Address bits", &gl_address_bits_override, NP, ProcessSineFreq, ME_AFT, ME_DEC},
    {"2nd byte", &gl_second_byte_override, NP, ProcessSineFreq, ME_AFT, ME_DEC},
#endif
    {NS},
};

VLIST tasks_vl[] =
{
    {"All tasks",        &gl_taskset[0],   NP, ClearTaskSet, ME_BEF, ME_DEC},
    {"Behavior",         &gl_taskset[1],   NP, ClearTaskSet, ME_BEF, ME_DEC},
    {"Behavior+TTL",     &gl_taskset[2],   NP, ClearTaskSet, ME_BEF, ME_DEC},
    {"Behavior+optical", &gl_taskset[3],   NP, ClearTaskSet, ME_BEF, ME_DEC},
    {"Fix+TTL",          &gl_taskset[4],   NP, ClearTaskSet, ME_BEF, ME_DEC},
    {"Fix+optical",      &gl_taskset[5],   NP, ClearTaskSet, ME_BEF, ME_DEC},
    {"Max # blocks",     &gl_ntrlspercond, NP, NP, 0, ME_DEC},
    {NS},
};

MENU umenus[] = {
    {"Fixation",    &fp_vl,    NP, NP, 0, NP, NP},
    {"Target",      &targ_vl,  NP, NP, 0, NP, NP},
    {"Stimulation", &stim_vl,  NP, NP, 0, NP, NP},
    {"Tasks",       &tasks_vl, NP, NP, 0, NP, NP},
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
        time 500 random 500
        to error on +WD0_XY & eyeflag
        to fpOff
    fpOff:
        code FPOFFCD
        do HideFP()
        to preTargDel
    preTargDel:
        do timerSet(1,100,&gl_fpoff_t)
        to rewEarlyOn on 3 > gl_curtasktype /* non-target tasks 0, 1, 2 */
        to showTarg on 1 % timerCheck
    showTarg:
        code TARGONCD
        do ShowTarg()
        to preTTLStimDel on 1 = gl_curstimtype
        to preOptStimDel on 2 = gl_curstimtype /* original laser interface */
        to preModStimDel on 3 = gl_curstimtype /* laser modulator */
        to rewHold /* tasktype 3 */
    rewEarlyOn:
        code REWCD
        do dio_on(REW)
        time 70
        to rewEarlyOff
    rewEarlyOff:
        do dio_off(REW)
        to waitState on 0 = gl_curtasktype
        to preTTLStimDel on 1 = gl_curstimtype
        to preOptStimDel on 2 = gl_curstimtype
        to preModStimDel
    preTTLStimDel:
        do timerSet(1,129,&gl_fpoff_t) /* forces a 130(-1 for stimThink) ms delay from gl_fpoff_t */
        to stimThink on 1 % timerCheck
    preOptStimDel:
        do timerSet(1,110,&gl_fpoff_t) /* need 130(-(2*10) send 10-bit string) ms between fp off and stim on */
        to setLaserIntensity on 1 % timerCheck
    preModStimDel:
        do timerSet(1,90,&gl_fpoff_t)  /* need 130(-(2*20) send 20-bit string) ms between fp off and stim on */
        to setLaserModParams on 1 % timerCheck
    stimThink:
        to error on +WD0_XY & eyeflag
        to finished on 1 = gl_done
        to ttlStimOn
    setLaserIntensity:
        do SendLaserIntensity(&gl_laser_power)
        time 100 /* this is the laser stimulation duration */
        to unsetLaserIntensity
    setLaserModParams:
        do SendModulatorParams(&gl_freq_idx) /* turn on laser drop STIMONCD */
        time 100 /* this is the laser stimulation duration */
        to turnOffModulator
    turnOffModulator:
        do turnOffModulator() /* turns off TTL3 and drops STIMOFFCD */
        time 30
        to rewHold
    unsetLaserIntensity:
        do SendLaserIntensity(NULL)
        time 30
        to rewHold
    ttlStimOn:
        code STIMONCD
        do SetStimState(&gl_ttlstimon,1)
        time 100 /* this is the TTL stimulation duration */
        to ttlStimOff
    ttlStimOff:
        code STIMOFFCD
        do SetStimState(&gl_ttlstimon,0)
        to rewHold
    rewHold:
        do timerSet(1,300,&gl_targon_t)
        to waitState on 3 > gl_curtasktype /* tasks 1, 2 */
        to rewLateOn on -WD1_XY & eyeflag
        to error on 1 % timerCheck /* no valid eye movement within 300 ms */
    waitState: /* 500 ms buffer time for neural response/saccades to phosphenes */
        do timerSet(1,630,&gl_fpoff_t)
        to correctTrial on 1 % timerCheck
    rewLateOn:
        code REWCD
        do dio_on(REW)
        time 70
        to rewLateOff
    rewLateOff:
        do dio_off(REW)
        to waitState
    error:
        code ABORTCD
        do ErrorTrial()
        time 200
        to dropTrialParams
    correctTrial:
        do CorrectTrial()
        to hideTarg on 2 < gl_curtasktype /* tasks 3-5, target trials */
        to dropTrialParams
    hideTarg:
        code TARGOFFCD
        do HideTarg()
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
    finished:
        do FreeMemory()
        time 100
        to pause0
    abort list:
}

ttlstim_set {
status ON
begin
    stimZero:
        to stimFirst
    stimFirst:
        do dio_off(FTTL2)
        to ttlOnPhase on 1 = gl_ttlstimon
    ttlOnPhase:
        do dio_on(FTTL2)
        to ttlPulseDurOn
    ttlPulseDurOn:
        to stimFirst on 0 = gl_ttlstimon
        to ttlOffPhase
    ttlOffPhase:
        do dio_off(FTTL2)
        to ttlPulseDurOff
    ttlPulseDurOff:
        to stimFirst on 0 = gl_ttlstimon
        to ttlOnPhase
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
