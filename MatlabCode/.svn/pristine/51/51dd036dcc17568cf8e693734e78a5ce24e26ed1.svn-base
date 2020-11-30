/*
 * A paradigm for running a white noise experiment.
 *
 * GDLH 11/28/07
 *
 * Adding the ability to stimulate subunits independently - ZALB Jan 2015
 */
#include <stdlib.h>
#include "ldev.h"
#include "labcodes.h"
#include "WNSubunit.h"
#include "UDPtransfer.h"
#include "rigconsts.h"

/**************************************/
/* Global variable declarations /*
/**************************************/

/* Which fundamentals to use.  Note: single quotes inside double quotes are important */
char gl_fundfilename [] = "'T_cones_smj10.mat'";

/* User changeable stuff */
int gl_fpx = 0;                       /* X position of fixation point */
int gl_fpy = 0;                       /* Y position of fixation point */
int gl_fpsize = 2;                    /* fixation point size in 0.1 deg units */
int gl_eyewinx = 5;                   /* Half width of eye window in 0.1 deg units */
int gl_eyewiny = 5;                   /* Half height of eye window in 0.1 deg units */
int gl_eyewinx_im = 3;                /* Half width of eye window during synth image trials */
int gl_eyewiny_im = 3;                /* Half height of eye window during synth image  trials */
int gl_rfx = 0;                       /* X position of RF in 0.1 deg units */
int gl_rfy = 0;                       /* Y position of RF in 0.1 deg units */
int gl_fprgb[] = {0, 0, 0};           /* RGB of the fixation point [0-255] */
int gl_mu_rgb[] = {0, 0, 0};          /* Mean of stimulus dist. relative to bkgnd (divided by 1000) */
int gl_mu_lms[] = {0, 0, 0};          /* Mean of stimulus dist. relative to bkgnd (divided by 1000) */
int gl_sigma_rgb[] = {130, 130, 130}; /* Standard deviation of stimulus dist. (divided by 1000)*/
// note: Changed from the default 150 in March 2013 (150s were out of gamut on both rigs) -- ZLB
int gl_sigma_lms[] = {0, 0, 0};       /* Standard deviation of stimulus dist. (divided by 1000)*/
int gl_gaussbounds[] = {5, 995};      /* probability cut offs for Gaussian (multiplied by 1000) */
long gl_seed = 1;                     /* Random number seed */
int gl_nstixperside_public = 10;      /* Number of stimulus elements (stixels) per side of white noise patch */
int gl_nstixperside;                  /* The "private" variable actually used internally */
int gl_npixperstix = 3;               /* Number of monitor pixels per side of each stixel */
int gl_nepochs = 30;                  /* Number of stimulus epochs */

int *gl_mask = NULL;
int gl_nmask_elem = 0;

/* Flags */
int gl_firsttrial = 1;                /* 1 =  We are on the first trial (header info should be dropped into Plexon) */
int gl_stimon = 0;                    /* Flag to indicate when the first screen('flip') occurs on the mac */
int gl_showimage = 0;                 /* 1 = Show a synthetic image on this trial  */
int gl_mac_ready = 0;                 /* the Mac sends a message when it has received the subunit mask */

/* Noise types */
int gl_noisetype = INTERLEAVE; // noise types defined in header

/**************************/
/* Stuff for TrialParams array */
/**************************/
/* Indices for the TrialParams */
#define TP_MU1        0
#define TP_MU2        1
#define TP_MU3        2
#define TP_SIGMA1     3
#define TP_SIGMA2     4
#define TP_SIGMA3     5
#define TP_BKGNDR     6
#define TP_BKGNDG     7
#define TP_BKGNDB     8
#define TP_NFRAMES    9
#define TP_SEED      10
#define TP_NOISETYPE 11
#define TP_IMAGE     12
#define NTRIALPARAMS 13

int gl_trialparamcodes[NTRIALPARAMS] = {MU1CD, MU2CD, MU3CD,
                                        SIGMA1CD, SIGMA2CD, SIGMA3CD,
                                        BKGNDRCD, BKGNDGCD, BKGNDBCD,
                                        NFRAMESCD, SEEDCD, NOISETYPECD,
                                        SYNTHIMAGECD};
int gl_trialparamtypes[NTRIALPARAMS] = {INT, INT, INT,
                                        INT, INT, INT,
                                        INT, INT, INT,
                                        INT, LONG, INT,
                                        LONGPOINTER};

ecode_struct gl_trialparams[NTRIALPARAMS];

/* Global variables for holding onto synthetic images and dropping them into the data file */
long gl_imagearray[1875];   /* assumes max image size is 25 x 25 (x 3) */
int gl_nimageelements = 0;

/*****************/
/* Request strings */
/*****************/
const char RFXY_REQ[] = "sendToRex(udpCom, gl.bar.xy, 'double', message);";
const char NFRAMES_REQ[] = "sendToRex(udpCom, gl.stim.framecounter, 'integer', message);";
const char FRAMERATE_REQ[] = "sendToRex(udpCom, gl.framerate, 'double', message);";
const char GAMMA_REQ[] = "sendToRex(udpCom, gl.cal.gammaTable, 'double', message);";
const char FUNDAMENTALS_REQ[] = "sendToRex(udpCom, gl.cal.fundamentals, 'double', message);";
const char MONSPD_REQ[] = "sendToRex(udpCom, gl.cal.monSpd, 'double', message);";
const char BKGNDRGB_REQ[] = "sendToRex(udpCom, gl.bkgndRGB, 'integer', message);";
const char SEED_REQ[] = "sendToRex(udpCom, gl.stim.seed, 'integer', message);";
const char IMAGE_REQ[] = "imagerequest();";
const char GABOR_REQ[] = "gaborrequest();";
const char MASK_REQ[] = "maskrequest();";

/********************/
/* Command strings   */
/********************/
const char IMAGE_ALERT[] = "imageReceiver(%d);"; /* Mac */
const char GABOR_ALERT[] = "gaborReceiver(%d);"; /* Mac */
const char MASK_ALERT[] = "maskReceiver(%d);"; /* Mac */
const char STARTSLAVE_ALERT[] = "WhiteNoise;"; /* Mac */
const char RETURN_ALERT[] = "return;";
const char CLEARSTATS_REQ[] = "ClearStats()"; /* Plexon */
const char STARTONLINE_ALERT[] = "WhiteNoiseOnline;"; /* Plexon */

/*********************/
/* Function definitions  */
/*********************/
int my_randint(int);

/*******************************/
/* Functions for dealing with UDP  */
/*******************************/

int
udpopen(void)
{
    rigconstants(whichrig());
    udp_open_2con(REXIP, MACIP, PLEXIP, UDPPORT);
    return 0;
}

int
udpclose(void)
{
    udp_close();  /* Closing an closed port is better than opening an open one. */
    return 0;
}

/***************************************/
/* Functions for sending requests via UDP  */
/***************************************/

int getFrameRate(void){
    udp_send_to(FRAMERATE_REQ, "MAC");
    return 0;
}

int getNumFrames(void){
    udp_send_to(NFRAMES_REQ, "MAC");
    return 0;
}

int getBkgndrgb(void){
    udp_send_to(BKGNDRGB_REQ, "MAC");
    return 0;
}

int getGammaTable(void){
    udp_send_to(GAMMA_REQ, "MAC");
    return 0;
}

int getFundamentals(void){
    udp_send_to(FUNDAMENTALS_REQ, "MAC");
    return 0;
}

int getMonSpd(void){
    udp_send_to(MONSPD_REQ, "MAC");
    return 0;
}

int getSeed(void){
    udp_send_to(SEED_REQ, "MAC");
    return 0;
}

int getMask()
{
    gl_mac_ready = 0;
    udp_send_to(MASK_REQ, "PLEXON");
    return 0;
}

/* Only send a request for a new image if there is not one currently waiting in the queue */
int getImage(void){

    if (!gl_nimageelements)
    {
        udp_send_to(IMAGE_REQ, "PLEXON");
    }
    return 0;
}

int getGabor(void){
    udp_send_to(GABOR_REQ, "PLEXON");
    return 0;
}

int prepImage(void){

    if (gl_nimageelements)
    {  /* This would be a better place to send the image to the mac, not in the listening loop!! */
        gl_showimage = 1;
        gl_trialparams[TP_IMAGE].val.lpval = (long *)&gl_imagearray;
        gl_trialparams[TP_IMAGE].nelements = gl_nimageelements;
        gl_trialparams[TP_NOISETYPE].val.lval = 0;


//      message.contents = &gl_imagearray; /* This won't work because gl_imagearry is longs not chars */
//      sprintf(buf, IMAGE_REQ);
//      message.header = &buf;
//      message -> size = gl_nimageelements;
//      sprintf(macfncall, IMAGE_ALERT, gl_nimageelements);
//      udpBulkTransfer(message, "MAC", outMsg);
    }
    return 0;
}

/* resizing the fixation window for synthetic image trials */
int shrinkeyewin(void)
{
    wd_siz(WIND0, gl_eyewinx_im, gl_eyewiny_im);
    return 0;
 }

/* Request the coordinates of the receptive field only if it hasn't already been set through the user menu */
int getRFxy(void) {
    if (gl_rfx == 0 & gl_rfy == 0)
    {
        udp_send_to(RFXY_REQ, "MAC");
    }
    return 0;
}

/********************************************/
/* Functions for dealing with incoming messages */
/********************************************/

int
checkmsg()
{
    char inMsg[REXBUFF];
    char outMsg[300];
    static message_struct message;
    int messageAvailable;
    double doublearray[790];
    double tmp;

    inMsg[0] = '\0';
    messageAvailable = udp_check(0);

    if (messageAvailable) {
        udp_read(inMsg, REXBUFF);
        prepareMsg(&message, inMsg);

        if (!strcmp(message.header, "MACSTIMON")) {
            gl_stimon = 1;
        } else if (!strcmp(message.header, "MASKRECV")) {
            gl_mac_ready = 1;
        } else if (!strcmp(message.header, "DISPLAYINIT")) {
            gl_mac_ready = 1;
        } else if (!strcmp(message.header, NFRAMES_REQ)) {
            gl_trialparams[TP_NFRAMES].val.lval  = (long)hex2long(message.contents, NULL);
        } else if (!strcmp(message.header, BKGNDRGB_REQ)) {
            hex2long(message.contents, (void *)&doublearray);
            gl_trialparams[TP_BKGNDR].val.lval = ((long *)&doublearray)[0];
            gl_trialparams[TP_BKGNDG].val.lval = ((long *)&doublearray)[1];
            gl_trialparams[TP_BKGNDB].val.lval = ((long *)&doublearray)[2];
        } else if (!strcmp(message.header, SEED_REQ)) {
            gl_seed = hex2long(message.contents, NULL);
        } else if (!strcmp(message.header, FRAMERATE_REQ)) {
            tmp = hex2double(message.contents, NULL);
            sendToPlexDataStream("double", (void *)&tmp, message.size, FRAMERATECD, 0);
        } else if (!strcmp(message.header, GAMMA_REQ)) {
            hex2double(message.contents, (void *)&doublearray);
            sendToPlexDataStream("double", (void *)&doublearray, message.size, GAMMATABLECD, 1);
        } else if (!strcmp(message.header, FUNDAMENTALS_REQ)) {
            hex2double(message.contents, (void *)&doublearray);
            sendToPlexDataStream("double", (void *)&doublearray, message.size, FUNDAMENTALSCD, 1);
        } else if (!strcmp(message.header, MONSPD_REQ)) {
            hex2double(message.contents, (void *)&doublearray);
            sendToPlexDataStream("double", (void *)&doublearray, message.size, MONSPDCD, 1);
        } else if (!strcmp(message.header, RFXY_REQ)) {
            hex2double(message.contents, (void *)&doublearray);
            if (message.size && (doublearray[0] || doublearray[1])) {
                gl_rfx = my_round(doublearray[0]*10);
                gl_rfy = my_round(doublearray[1]*10);
            }
        } else if (!strcmp(message.header, IMAGE_REQ)) {
            if (message.size > 0) {
                hex2long(message.contents, (long *)&gl_imagearray);
                gl_nimageelements = message.size;  /* This is the indication that there is an image ready to be shown */
                sprintf(outMsg, IMAGE_ALERT, gl_nimageelements);
                udpBulkTransfer(message, "MAC", outMsg);
            }
        } else if (!strcmp(message.header, GABOR_REQ)) {
            if (message.size > 0) {
                hex2double(message.contents, (void *)&doublearray);
                sprintf(outMsg, GABOR_ALERT, message.size);
                udpBulkTransfer(message, "MAC", outMsg);
            }
        } else if (!strcmp(message.header, MASK_REQ)) {
            if (message.size > 0) {
                free(gl_mask);
                gl_nmask_elem = message.size;
                gl_mask = malloc(sizeof(int)*gl_nmask_elem);
                hex2long(message.contents, gl_mask);
                sprintf(outMsg, MASK_ALERT, message.size);
                udpBulkTransfer(message, "MAC", outMsg);
            }
        }
    }
    return 0;
}

/*******************************************************/
/* Functions for dealing with the header and trial parameters */
/*******************************************************/

/* The first pair of codes that we drop are the magic "8999" code, which indicates that */
/* the next code is the paradigm identifier, and then the paradigm identifier. */

int
drophdr(void)
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
    ec_send_code_tagged(NPIXPERSTIXCD, VALOFFSET+gl_npixperstix);
    ec_send_code_tagged(NSTIXPERSIDECD, VALOFFSET+gl_nstixperside);
    ec_send_code_tagged(GAUSSLOCUTCD, VALOFFSET+gl_gaussbounds[0]);
    ec_send_code_tagged(GAUSSHICUTCD, VALOFFSET+gl_gaussbounds[1]);

    gl_firsttrial = 0;

    return 0;
}

int
inittrialparams(void)
{
    int i;

    for (i=0; i < NTRIALPARAMS; i++) {
        gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
        gl_trialparams[i].casttype = gl_trialparamtypes[i];
        gl_trialparams[i].nelements = 0;
        gl_trialparams[i].val.dval = 0;
        gl_trialparams[i].val.lpval = NULL;
    }
    gl_trialparams[TP_SEED].val.lval = gl_seed;  // Storing the seed away before it gets changed

    if (gl_noisetype != INTERLEAVE) {
        gl_trialparams[TP_NOISETYPE].val.lval = gl_noisetype;
    } else {
        gl_trialparams[TP_NOISETYPE].val.lval = (long) (my_randint(2))+1;  // interleave by default
        if (gl_sigma_lms[0] == 0 && gl_sigma_lms[1] == 0 && gl_sigma_lms[2] == 0) {
            gl_trialparams[TP_NOISETYPE].val.lval = GUN;  // Gun noise only
        }
        if (gl_sigma_rgb[0] == 0 && gl_sigma_rgb[1] == 0 && gl_sigma_rgb[2] == 0) {
            gl_trialparams[TP_NOISETYPE].val.lval = CONE;  // Cone noise only
        }
    }
    return 0;
}

/* Drop the trial parameters into the data file using sendToPlexDataStream */
/* By convention, trial parameters are never bookended */
int
droptrialparams(void)
{
    int i;
    int ncodes = 0;
    float mspercode = 0.5;   /* Conservative; actual value is closer to 0.33 */

    if (gl_nmask_elem) {
        sendToPlexDataStream("int", (void *) gl_mask, gl_nmask_elem, MASKCD, 1);
        ncodes = ncodes + gl_nmask_elem + 2;
    }

    for (i = 0; i < NTRIALPARAMS; i++)
    {
        if (gl_trialparams[i].casttype == INT)
        {
            sendToPlexDataStream("int", (void *)&gl_trialparams[i].val.lval, 1, gl_trialparams[i].flagcodeval, 0);
            ncodes = ncodes + 2;
        }
        if (gl_trialparams[i].casttype == LONG)
        {
            sendToPlexDataStream("long", (void *)&gl_trialparams[i].val.lval, 1, gl_trialparams[i].flagcodeval, 0);
            ncodes = ncodes + 5;
        }
        if (gl_trialparams[i].casttype == DOUBLE)
        {
            sendToPlexDataStream("double", (void *)&gl_trialparams[i].val.dval, 1, gl_trialparams[i].flagcodeval, 0);
            ncodes = ncodes + 9;
        }
        if (gl_trialparams[i].casttype == LONGPOINTER && gl_trialparams[i].val.lpval != NULL)
        {
            sendToPlexDataStream("long", gl_trialparams[i].val.lpval, gl_trialparams[i].nelements, gl_trialparams[i].flagcodeval, 1);
            ncodes = ncodes + 2 + 4*gl_trialparams[i].nelements;
        }
    }

    set_times("droptrialparams", my_ceil(ncodes*mspercode), 0);
    return 0;
}

/*
 * The way the logic works, we put gl_seed in gl_trialparams[TP_SEED] at the
 * begining of each trial.  Then, after the stimulus is displayed, we query the
 * Mac for the new seed and put that in gl_seed.  gl_trialparams[TP_SEED] is
 * unchanged so we can drop it into the data file.
 */

/*
* Initializations.
*/
void
rinitf(void)
{
    wd_disp(D_W_ALLCUR & ~D_W_JOY);     /* all cursors but joystick */
    wd_src_pos(WIND0, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_check(WIND0, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);
    wd_cntrl(WIND0, WD_ON);
    update_eyewin();
    gl_mac_ready = 0;
    gl_firsttrial = 1;      /* So reset_state command drops the header again */
    gl_stimon = 0;
    gl_showimage = 0;
    gl_nimageelements = 0;
    gl_nstixperside = gl_nstixperside_public;
    udp_send_to("AllOff()", "MAC");
    srand(time(NULL));
    return 0;
}

/************************/
/*    Plexon commands      */
/************************/

int
startonlineprgm(void)
{
    udp_send_to(STARTONLINE_ALERT, "PLEXON");
    return 0;
}

int
endonlineprgm(void)
{
    udp_send_to(RETURN_ALERT, "PLEXON");
    return 0;
}

int
initonlineprgm(void)
{
    udp_send_to(CLEARSTATS_REQ, "PLEXON");
    return 0;
}

/************************/
/* Mac plotting commands */
/************************/

int
startslaveprgm(void)
{
    udp_send_to(STARTSLAVE_ALERT, "MAC");
    return 0;
}

int
endslaveprgm(void)
{
    udp_send_to(RETURN_ALERT, "MAC");
    return 0;
}

int
setup(void)
{
    char buf[500];
    gl_mac_ready = 0;
    sprintf(buf, "InitDisplay(%d, %d, %s, %s)",
                  GL_MONDIST, GL_SCREENWIDTH, GL_CALFILENAME, gl_fundfilename);
    udp_send_to(buf, "MAC");
    return 0;
}

int
showFP(void)
{
    char buf[256];

    sprintf(buf, "ShowFP(%d, %d, %d, %d, %d, %d)",
            gl_fpx, gl_fpy, gl_fpsize, gl_fprgb[0], gl_fprgb[1], gl_fprgb[2]);
    udp_send_to(buf, "MAC");
    update_eyewin();
    return 0;
}

int
showStim(void)
{
    char buf[256];
    int mu[3];
    int sigma[3];

    gl_stimon = 0;  /* Initializing.  Will be set to '1' upon receipt of acknowledgment */

    switch (gl_trialparams[TP_NOISETYPE].val.lval) {
    case GUN:
        mu[0] = gl_mu_rgb[0];
        mu[1] = gl_mu_rgb[1];
        mu[2] = gl_mu_rgb[2];
        sigma[0] = gl_sigma_rgb[0];
        sigma[1] = gl_sigma_rgb[1];
        sigma[2] = gl_sigma_rgb[2];
        break;
    case RGCONE:
    case CONE:
        mu[0] = gl_mu_lms[0];
        mu[1] = gl_mu_lms[1];
        mu[2] = gl_mu_lms[2];
        sigma[0] = gl_sigma_lms[0];
        sigma[1] = gl_sigma_lms[1];
        sigma[2] = gl_sigma_lms[2];
        break;
    default:
        dprintf("Unknown noise type (%d); assuming gun noise\n", gl_trialparams[TP_NOISETYPE].val.lval);
        mu[0] = gl_mu_rgb[0];
        mu[1] = gl_mu_rgb[1];
        mu[2] = gl_mu_rgb[2];
        sigma[0] = gl_sigma_rgb[0];
        sigma[1] = gl_sigma_rgb[1];
        sigma[2] = gl_sigma_rgb[2];
    }

    gl_trialparams[TP_MU1].val.lval = mu[0];
    gl_trialparams[TP_MU2].val.lval = mu[1];
    gl_trialparams[TP_MU3].val.lval = mu[2];
    gl_trialparams[TP_SIGMA1].val.lval = sigma[0];
    gl_trialparams[TP_SIGMA2].val.lval = sigma[1];
    gl_trialparams[TP_SIGMA3].val.lval = sigma[2];
    sprintf(buf, "ShowStim(%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %ld)",
            gl_rfx, gl_rfy, gl_seed, gl_npixperstix, gl_nstixperside, mu[0], mu[1], mu[2],
            sigma[0], sigma[1], sigma[2], gl_gaussbounds[0], gl_gaussbounds[1],
            gl_trialparams[TP_NOISETYPE].val.lval);
    udp_send_to(buf, "MAC");
    return 0;
}

int
showImage(void)
{
    gl_stimon = 0;  /* Initializing.  Will be set to '1' upon receipt of acknowledgment */

    udp_send_to("ShowImage()", "MAC");
    return 0;
}


int
imageoff(void)
{
    gl_stimon = 0;

    udp_send_to("ImageOff()", "MAC");
    return 0;
}

int
successfulimagetrial(void)
{
    gl_showimage = 0;
    gl_nimageelements = 0;      /* This is the indication that we should ask for a new image */
    /* Don't clear out gl_imagearray yet - we need it around for when we drop the trial parameters */

    return 0;
}

int
checkStimDone(int flag)
{
    static int ntimes;

    if (flag == 0)
        ntimes = 0;
    if (++ntimes > gl_nepochs)
        gl_stimon = 0;

    return 0;
}

int
alloff(void)
{
    udp_send_to("AllOff()", "MAC");
    gl_stimon = 0;
    return 0;
}

int
giverew(void)
{
    dio_on(REW);
    return(REWCD);
}

int
update_eyewin()
{
    wd_src_pos(WIND0, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_pos(WIND0, gl_fpx, gl_fpy);
    wd_siz(WIND0, gl_eyewinx, gl_eyewiny);
    return 0;
}

// generate a random integer on [0,n-1] without modulo bias
// http://stackoverflow.com/questions/10984974/why-do-people-say-there-is-modulo-bias-when-using-a-random-number-generator#comment18071601_10984975
int my_randint(int n)
{
    int limit = RAND_MAX - ((RAND_MAX % n) + 1) % n;
    int rnd;

    do { rnd = rand(); } while ( rnd > limit );
    return rnd % n;
}

/*
* User accessible menu
*/
VLIST fp_vl[] =
{
    {"Fixation X",       &gl_fpx,        NP, NP, 0, ME_DEC},
    {"Fixation Y",       &gl_fpy,        NP, NP, 0, ME_DEC},
    {"Fixation size",    &gl_fpsize,     NP, NP, 0, ME_DEC},
    {"Fixation R",       &gl_fprgb[0],   NP, NP, 0, ME_DEC},
    {"Fixation G",       &gl_fprgb[1],   NP, NP, 0, ME_DEC},
    {"Fixation B",       &gl_fprgb[2],   NP, NP, 0, ME_DEC},
    {"Eyewin width",     &gl_eyewinx,    NP, NP, 0, ME_DEC},
    {"Eyewin height",    &gl_eyewiny,    NP, NP, 0, ME_DEC},
    {"Eyewin width im",  &gl_eyewinx_im, NP, NP, 0, ME_DEC},
    {"Eyewin height im", &gl_eyewiny_im, NP, NP, 0, ME_DEC},
    {NS},
};

VLIST stimstuff_vl[] =
{
    {"nstim_epochs",    &gl_nepochs,        NP, NP, 0, ME_DEC},
    {"xpos",            &gl_rfx,            NP, NP, 0, ME_DEC},
    {"ypos",            &gl_rfy,            NP, NP, 0, ME_DEC},
    {"nstixperside",    &gl_nstixperside_public,   NP, NP, 0, ME_DEC},
    {"npixperstix",     &gl_npixperstix,    NP, NP, 0, ME_DEC},
    {"gauss lo cutoff", &gl_gaussbounds[0], NP, NP, 0, ME_DEC},
    {"gauss hi cutoff", &gl_gaussbounds[1], NP, NP, 0, ME_DEC},
    {"random seed",     &gl_seed,           NP, NP, 0, ME_DEC},
    {"noise type",      &gl_noisetype,      NP, NP, 0, ME_DEC},
    {NS},
};

VLIST stimstats_vl[] =
{
    {"mu_red",      &gl_mu_rgb[0],    NP, NP, 0, ME_DEC},
    {"mu_green",    &gl_mu_rgb[1],    NP, NP, 0, ME_DEC},
    {"mu_blue",     &gl_mu_rgb[2],    NP, NP, 0, ME_DEC},
    {"sigma_red",   &gl_sigma_rgb[0], NP, NP, 0, ME_DEC},
    {"sigma_green", &gl_sigma_rgb[1], NP, NP, 0, ME_DEC},
    {"sigma_blue",  &gl_sigma_rgb[2], NP, NP, 0, ME_DEC},
    {"mu_L",        &gl_mu_lms[0],    NP, NP, 0, ME_DEC},
    {"mu_M",        &gl_mu_lms[1],    NP, NP, 0, ME_DEC},
    {"mu_S",        &gl_mu_lms[2],    NP, NP, 0, ME_DEC},
    {"sigma_L",     &gl_sigma_lms[0], NP, NP, 0, ME_DEC},
    {"sigma_M",     &gl_sigma_lms[1], NP, NP, 0, ME_DEC},
    {"sigma_S",     &gl_sigma_lms[2], NP, NP, 0, ME_DEC},
    {NS},
};

MENU umenus[] =
{
    {"FixpointStuff", &fp_vl       , NP, NP, 0, NP, NP},
    {"StimulusStuff", &stimstuff_vl, NP, NP, 0, NP, NP},
    {"StimulusStats", &stimstats_vl, NP, NP, 0, NP, NP},
    {NS}
};

/*
 * User-supplied real-time variable table
 */
RTVAR rtvars[] = {};

%%
id PARADIGMID
restart rinitf
main_set {
status ON
begin
first:
    to closeport
closeport:
    do udpclose()
    to openport
openport:
    do udpopen()
    to endslave
endslave:
    do endslaveprgm()
    to endonline
endonline:
    do endonlineprgm()
    time 60
    to startslave
startslave:
    do startslaveprgm()
    time 1000
    to initmac
initmac:
    do setup()
    time 5000
    to startonline on 1 = gl_mac_ready
    to failSafely0
startonline:
    do startonlineprgm()
    time 50
    to getRFxy
getRFxy:
    do getRFxy()
    time 10
    to pause0
pause0:
    to pause1 on +PSTOP & softswitch
    to inittrial
pause1:
    to inittrial on -PSTOP & softswitch
inittrial:
    do inittrialparams()
    to plxstart
plxstart:
    do dio_on(PLXREC)
    time 400
    to drophdr on 1 = gl_firsttrial
    to getImage
drophdr:
    do drophdr()
    time 300
    to getFrameRate
getFrameRate:
    do getFrameRate()
    time 10
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
    time 800
    to fpon
getImage:
    do getImage()  /* sending the image request */
    time 10
    to getMask
getMask:
    do getMask()
    time 5000
    to fpon on 1 = gl_mac_ready
    to error
fpon:
    code FPONCD
    do showFP()
    time 2000
    to prestimdel on -WD0_XY & eyeflag
    to error
prestimdel:
    time 450
    to error on +WD0_XY & eyeflag
    to prepimage
prepimage:
    do prepImage()  /* Do we have an image ready to show? */
    to prestimswitch
prestimswitch:
    to error on +WD0_XY & eyeflag
    to shrinkeyewin on 1 = gl_showimage
    to showStim
showStim:
    do showStim()
    to stimcode on 1 = gl_stimon   /* Mac has sent the "stimulus on" signal */
stimcode:
    code STIMONCD
    do checkStimDone(0)   /* Initializing reward counter */
    to stimpause
stimpause:
    time 810 random 0
    to alloff on +WD0_XY & eyeflag
    to checkdone
checkdone:
    do checkStimDone(1)    /* Are we done yet?  If yes, set gl_stimon = 0 */
    to rewon
shrinkeyewin:
    do shrinkeyewin()
    time 100   /* a brief grace period */
    to showImage
showImage:
    do showImage()
    to error on +WD0_XY & eyeflag
    to imagecode on 1 = gl_stimon
imagecode:
    code STIMONCD
    time 200
    to error on +WD0_XY & eyeflag
    to imageoff
imageoff:
    code STIMOFFCD
    do imageoff()
    to postimagedelay
postimagedelay:
    time 100
    to error on +WD0_XY & eyeflag
    to successfulimagetrial
successfulimagetrial:
    time 300
    do successfulimagetrial()
    to rewon
rewon:
    do giverew()
    time 45
    to rewoff
rewoff:
    do dio_off(REW)
    to alloff on  0 = gl_stimon
    to stimpause
error:             /* abort trial */
    code ABORTCD
    to alloff
alloff:
    code ALLOFFCD
    do alloff()
    to getnframes
getnframes:
    do getNumFrames()
    time 20     /* Need time to let the request go through */
    to getBkgndrgb
getBkgndrgb:
    do getBkgndrgb()
    time 20
    to getSeed
getSeed:
    do getSeed()
    time 20
    to droptrialparams
droptrialparams:
    do droptrialparams()
    to marktrialend
marktrialend:
    code EOT
    time 10
    to getGabor
getGabor:
    do getGabor()  /* sending the gabor parameters request */
    to plxstop
plxstop:
    do dio_off(PLXREC)
    time 140
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
msgFirst:
    to msgSecond
msgSecond:
    do checkmsg()
    to msgSecond
abort list:
}
