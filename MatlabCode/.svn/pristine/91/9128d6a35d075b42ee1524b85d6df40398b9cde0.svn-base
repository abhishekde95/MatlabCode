/*
 * This paradigm finds the subunits of a RF, and stimulates each subunit
 * with different combinations of L- and M-cone contrast.  Version 2.0
 */

#include <stdio.h>
#include "ldev.h"
#include "labcodes.h"
#include "UDPtransfer.h"
#include "rigconsts.h"
#include "GridLMSubunitDetection.h"

// Globals, define initial values in rinitf
params_public gPublic;
params_private gPrivate;

/////////////////////////////////
// Stuff for TrialParams array //
/////////////////////////////////
// Indices for gl_trialparams

#define TP_LCC          0
#define TP_MCC          1
#define TP_SCC          2
#define TP_NSTIXGRID    3
#define TP_DVAPERSTIX   4
#define TP_STIMDUR      5
#define TP_RFCORRECT	6
#define TP_CHOOSECORR   7
#define NTRIALPARAMS    8

#define PI 3.141592653589793
#define WINDFIXN 		0
#define WINDCORR		1
#define WINDWRNG		2

int gl_trialparamcodes[NTRIALPARAMS] = {LCCCD, MCCCD, SCCCD, NSTIXGRIDCD, DVAPERSTIXCD, STIMDURCD, RFCORRCD, CORRECTCD};
int gl_trialparamtypes[NTRIALPARAMS] =  {FLOAT, FLOAT, FLOAT, INT, FLOAT, FLOAT, INT, INT};
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

const char RFXY_REQ[] = "sendToRex(udpCom, gl.bar.xy, 'double', message);";
const char FUNDFILENAME[] = "'/Applications/Psychtoolbox/PsychColorimetricData/PsychColorimetricMatFiles/T_cones_smj10.mat'";
const char FRAMERATE_REQ[] = "sendToRex(udpCom, gl.framerate, 'double', message);";
const char PIXPERDEG_REQ[] = "sendToRex(udpCom, gl.pixperdeg, 'double', message);";
const char BKGNDRGB_REQ[] = "sendToRex(udpCom, gl.bkgndrgb, 'double', message);";
const char GAMMA_REQ[] = "sendToRex(udpCom, gl.cal.gammaTable, 'double', message);";
const char FUNDAMENTALS_REQ[] = "sendToRex(udpCom, gl.cal.fundamentals, 'double', message);";
const char MONSPD_REQ[] = "sendToRex(udpCom, gl.cal.monSpd, 'double', message);";
const char TRIALPARAMS_REQ[] = "sendToRex(udpCom, getStimParams(), 'double', message);";
const char XCOORDINATES_REQ[] = "sendToRex(udpCom, getX(), 'int', message);";
const char YCOORDINATES_REQ[] = "sendToRex(udpCom, getY(), 'int', message);";

#define DECL_UDPSEND_FUNC( funcname, reqstr, dest ) int funcname##() { udp_send_to(reqstr, dest); return 0; }
DECL_UDPSEND_FUNC( stopSlavePrgm, "return;", "MAC" );
DECL_UDPSEND_FUNC( macAllOff, "AllOff();", "MAC" );
DECL_UDPSEND_FUNC( stopOnlinePrgm, "return;", "PLEXON" );

#undef DECL_UDPSEND_FUNC


///////////////////////////
//  Mac commands via UDP //
///////////////////////////

int startSlavePrgm()
{
	gPrivate.macready = 0;
	udp_send_to("GLMSDetection_Slave();", "MAC");
	return 0;
}

int getFrameRate()
{
	gPrivate.macready = 0;
	udp_send_to(FRAMERATE_REQ, "MAC");
	return 0;
}

int getPixPerDeg()
{
	gPrivate.macready = 0;
	udp_send_to(PIXPERDEG_REQ, "MAC");
	return 0;
}

int getBkgndrgb()
{
	gPrivate.macready = 0;
	udp_send_to(BKGNDRGB_REQ, "MAC");
	return 0;
}

int getGamma()
{
	gPrivate.macready = 0;
	udp_send_to(GAMMA_REQ, "MAC");
	return 0;
}

int getFundamentals()
{
	gPrivate.macready = 0;
	udp_send_to(FUNDAMENTALS_REQ, "MAC");
	return 0;
}

int getMonspd()
{
	gPrivate.macready = 0;
	udp_send_to(MONSPD_REQ, "MAC");
	return 0;
}

int getRFxy()
{
	gPrivate.plexready = 0;
    if ( gPublic.rf_x == 0 && gPublic.rf_y == 0 ) udp_send_to(RFXY_REQ, "PLEX");
    dprintf ("Restarting paradigm...\n\n\n");
    return 0;
}

int setupMac()
{
    char buf[256];
	gPrivate.macready = 0;
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
    updateEyewin(WIND0, gPublic.fp_x, gPublic.fp_y,
            gPublic.fixwin_w, gPublic.fixwin_h);
    gPrivate.fpon = 1;
	//gPrivate.fpoff = 0;
    return 0;
}

int hideFP()
{
    //char buf[256];
    gPrivate.macready = 0;
    //gPrivate.plexready = 0;
    udp_send_to("HideFP()", "MAC");
    gPrivate.fpon = 0;
    //gPrivate.fpoff = 1;
    return 0;
}

int sendX()
{
    char argsbuf[1024], reqbuf[1100];
	gPrivate.macready = 0;
    sprintf(reqbuf, "sendX([%s])", format_grid_request(argsbuf, gPublic.grid_x));
    udp_send_to(reqbuf, "MAC");
    return 0;
}

int sendY()
{
    char argsbuf[1024], reqbuf[1100];
    gPrivate.macready = 0;
	sprintf(reqbuf, "sendY([%s])", format_grid_request(argsbuf, gPublic.grid_y));
    udp_send_to(reqbuf, "MAC");
    return 0;
}

char * format_grid_request(char * buffer, int * grid)
{
    int i;
    char argbuf[8];
    buffer[0] = 0;
    for ( i = 0; i < gPrivate.grid_len; ++i ) {
        sprintf(argbuf, "%s%d", (i==0?"":" "), grid[i]);
        strcat(buffer, argbuf);
    }
    return buffer;
}

int prepareStim()
{
    char buf[256];
	double tempX;
    gPrivate.macready = 0;
    gPrivate.plexready = 0;
	tempX = gPrivate.rfcorrect ? gPublic.rf_x : -gPublic.rf_x;
    sprintf(buf, "PrepareStim(%d, %f, %f, %.15f, %.15f, %.15f, %.15f, %.15f, %d)",
        gPublic.nstixgrid, gPublic.dvaperstix, gPublic.stimdur,
        gPublic.stim_lmscc[0], gPublic.stim_lmscc[1], gPublic.stim_lmscc[2],
		tempX, gPublic.rf_y, gPrivate.rfcorrect);
    udp_send_to(buf, "MAC");
    gPrivate.firsttrial = 0;
    return 0;
}

int displayStim()
{
    gPrivate.macready = 0;
	gPrivate.stimon = 1;
    udp_send_to("displayStim()", "MAC");
    return 0;
}

int showTargs()
{
    char buf[64];
    double scale;
	int correct_targx;

	correct_targx = gPrivate.rfcorrect ? gPublic.rf_x : -gPublic.rf_x;

    updateEyewin(WINDCORR, correct_targx, gPublic.rf_y, gPublic.targwin_x, gPublic.targwin_y);
    updateEyewin(WINDWRNG,-correct_targx, gPublic.rf_y, gPublic.targwin_x, gPublic.targwin_y);
    sprintf(buf, "ShowTargs(%d, %d, %d)", correct_targx, gPublic.rf_y, gPublic.targ_size);
    udp_send_to(buf, "MAC");
    return 0;
}


///////////////////////////////
//  Plexon commands via UDP  //
//////////////////////////////

int startOnlinePrgm()
{
	gPrivate.plexready = 0;
	udp_send_to("GridLMSubunitDetectionOnline();", "PLEXON");
	return 0;
}

int getX()
{
    //char buf[256];
    //sprintf(buf, XCOORDINATES_REQ);
    gPrivate.plexready = 0;
    udp_send_to(XCOORDINATES_REQ,"PLEXON");
    return 0;
}

int getY()
{
    //char buf[256];
    //sprintf(buf, YCOORDINATES_REQ);
    gPrivate.plexready = 0;
    udp_send_to(YCOORDINATES_REQ,"PLEXON");
    return 0;
}

int getStimParams()
{
    gPrivate.plexready = 0;
    udp_send_to(TRIALPARAMS_REQ,"PLEXON");
    return 0;
}

//////////////////////////
// Getting UDP commands //
//////////////////////////

int checkMsg()
{
    char inMsg[REXBUFF]; // 8500
	char buf[1024];
    static message_struct message;
    int messageAvailable;
    double doublearray[790];
    double tmpd;
    //long longarray[600];
    int i;

    inMsg[0] = '\0';
    messageAvailable = udp_check(0);

    if ( messageAvailable ) {
        udp_read(inMsg, REXBUFF);
        prepareMsg(&message, inMsg);
		
        if ( !strcmp(message.header, "macready") ) {
            gPrivate.macready = 1;
        } else if ( !strcmp(message.header, "plexready") ) {
            gPrivate.plexready = 1;
        } else if ( !strcmp(message.header, RFXY_REQ) ) {
            hex2double(message.contents, (void*) &doublearray);
            if ( message.size && (doublearray[0] || doublearray[1]) ) {
                gPublic.rf_x = (float) doublearray[0];
                gPublic.rf_y = (float) doublearray[1];
            }
        } else if ( !strcmp(message.header, FRAMERATE_REQ) ) {
            tmpd = hex2double(message.contents, NULL);
            sendToPlexDataStream("double", (void *)&tmpd, message.size, FRAMERATECD, 0);
			dprintf("Got Frame Rate\n");
        } else if ( !strcmp(message.header, BKGNDRGB_REQ) ) {
            hex2double(message.contents, (void *)&doublearray);
            sendToPlexDataStream("double", (void *)&doublearray, message.size, BKGNDRGBCD, 1);
			gPrivate.bkgndrgb[0] = doublearray[0];
			gPrivate.bkgndrgb[1] = doublearray[1];
			gPrivate.bkgndrgb[2] = doublearray[2];
            dprintf("Got Background RGB\n");
        } else if ( !strcmp(message.header, PIXPERDEG_REQ) ) {
            tmpd = hex2double(message.contents, NULL);
            sendToPlexDataStream("double", (void *)&tmpd, message.size, PIXPERDEGCD, 0);
            dprintf("Got Pixels Per Degree\n");
        } else if ( !strcmp(message.header, GAMMA_REQ) ) {
            hex2double(message.contents, (void *)&doublearray);
            sendToPlexDataStream("double", (void *)&doublearray, message.size, GAMMATABLECD, 1);
            dprintf("Got Gamma\n");
        } else if ( !strcmp(message.header, FUNDAMENTALS_REQ) ) {
            hex2double(message.contents, (void *)&doublearray);
            sendToPlexDataStream("double", (void *)&doublearray, message.size, FUNDAMENTALSCD, 1);
            dprintf("Got Fundamentals\n");
        } else if ( !strcmp(message.header, MONSPD_REQ) ) {
            hex2double(message.contents, (void *)&doublearray);
            sendToPlexDataStream("double", (void *)&doublearray, message.size, MONSPDCD, 1);
            dprintf("Got MonSPD\n");
        } else if ( !strcmp(message.header, TRIALPARAMS_REQ) ) {
            hex2double(message.contents, (void *)&doublearray);
			gPublic.stim_lmscc[0] = doublearray[0];
            gPublic.stim_lmscc[1] = doublearray[1];
            gPublic.stim_lmscc[2] = doublearray[2];
            gPublic.dvaperstix = doublearray[3];
            gPublic.stimdur = doublearray[4];
			gPublic.nstixgrid = (int) doublearray[5];
			gPrivate.rfcorrect = (int) doublearray[6];
			sprintf(buf, "%f %f %f %f %f %d %d\n", gPublic.stim_lmscc[0], gPublic.stim_lmscc[1], gPublic.stim_lmscc[2], gPublic.dvaperstix, gPublic.stimdur, gPublic.nstixgrid, gPrivate.rfcorrect);
            gPrivate.plexready = 1;
        } else if ( !strcmp(message.header, XCOORDINATES_REQ)) {
            hex2long(message.contents, (void *)&gPublic.grid_x);
            gPrivate.grid_len = message.size;
            gPrivate.plexready = 1;
        } else if ( !strcmp(message.header, YCOORDINATES_REQ)) {
            hex2long(message.contents, (void *)&gPublic.grid_y);
            gPrivate.getnewstim = 0;
            gPrivate.plexready = 1;
		} else {
			dprintf("junk header: %s", message.header);
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
    ec_send_code_tagged(EYEWINXCD       , VALOFFSET + gPublic.fixwin_w);
    ec_send_code_tagged(EYEWINYCD       , VALOFFSET + gPublic.fixwin_h);
	
    sendToPlexDataStream("float",(void *)&gPublic.rf_x, 1,RFXCD, 0);
    sendToPlexDataStream("float",(void *)&gPublic.rf_y, 1,RFYCD, 0);

	ec_send_code_tagged(HEADERCMPLTCD       , VALOFFSET + 1);
	
    gPrivate.firsttrial = 0;
	gPrivate.rexready = 1;
    return 0;
}

int initTrialParams()
{
    int i;
    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
        gl_trialparams[i].casttype = gl_trialparamtypes[i];
        gl_trialparams[i].nelements = 0;
        gl_trialparams[i].val.dval = 0;
        gl_trialparams[i].val.lpval = NULL;
    }
    return 0;
}

int dropTrialParams()
{
    int i, ncodes = 0;
    float mspercode = 0.7; // Conservative; actual value is closer to 0.33

    gPrivate.plexready = 0;

    gl_trialparams[TP_LCC].val.fval = gPublic.stim_lmscc[0];
    gl_trialparams[TP_MCC].val.fval = gPublic.stim_lmscc[1];
    gl_trialparams[TP_SCC].val.fval = gPublic.stim_lmscc[2];
    gl_trialparams[TP_NSTIXGRID].val.lval = gPublic.nstixgrid;
    gl_trialparams[TP_DVAPERSTIX].val.fval = gPublic.dvaperstix;
    gl_trialparams[TP_STIMDUR].val.fval = gPublic.stimdur;
	gl_trialparams[TP_RFCORRECT].val.lval = gPrivate.rfcorrect;

    sendToPlexDataStream("int", (void*) &gPublic.grid_x, gPrivate.grid_len, GRIDXCD, 1);
    sendToPlexDataStream("int", (void*) &gPublic.grid_y, gPrivate.grid_len, GRIDYCD, 1);
    ncodes = ncodes + gPrivate.grid_len*2 + 4;

    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        if ( gl_trialparams[i].casttype == INT ) {
            sendToPlexDataStream("int", (void *) &gl_trialparams[i].val.lval, 1, gl_trialparams[i].flagcodeval, 0);
            ncodes = ncodes + 2;
        }
		if (gl_trialparams[i].casttype == LONG)
        {
            sendToPlexDataStream("long", (void *)&gl_trialparams[i].val.lval, 1, gl_trialparams[i].flagcodeval, 0);
            ncodes = ncodes + 5;
        }
        if ( gl_trialparams[i].casttype == FLOAT ) {
            sendToPlexDataStream("long", (void *) &gl_trialparams[i].val.fval, 1, gl_trialparams[i].flagcodeval, 0);
            ncodes = ncodes + 5;
            // dprintf("dropping a float: %f,\n", gl_trialparams[i].val.fval);
        }
        if ( gl_trialparams[i].casttype == DOUBLE ) {
            sendToPlexDataStream("double", (void *) &gl_trialparams[i].val.dval, 1, gl_trialparams[i].flagcodeval, 0);
            ncodes = ncodes + 9;
            // dprintf("dropping a double: %f,\n", gl_trialparams[i].val.dval);
        }
    }
		
    set_times("dropTrialParams", my_ceil(ncodes*mspercode), 0);
    return 0;
}


/////////////////////////
//   REX setup stuff   //
/////////////////////////

void initGlobalParams() // call this AFTER freeing memory to any dynamically allocated members
{
    static int firstTime = 1;
    int i;
    if ( firstTime ) { // here's where you set up the defaults
        firstTime = 0;
        
		// public
        gPublic.fixwin_w = 5; 
		gPublic.fixwin_h = 5;
        gPublic.fp_size = 2;
		gPublic.targwin_x = 13;
		gPublic.targwin_y = 13;
		gPublic.targ_size = 2;
		gPublic.rf_x = 0; 
		gPublic.rf_y = 0;
        gPublic.stim_lmscc[0] = 0; 
		gPublic.stim_lmscc[1] = 0; 
		gPublic.stim_lmscc[2] = 0;
        gPublic.nstixgrid = 0; // In DVA
        gPublic.dvaperstix = .2;
        gPublic.stimdur = .2; // In seconds

        for (i = 0; i < 600; ++i) { // this size is defined in the header
            gPublic.grid_x[i] = 0;
            gPublic.grid_y[i] = 0;
        }

    } else {
	
        // don't clear out gPublic
        memset(&gPrivate, 0, sizeof(gPrivate));

    }
    macAllOff();
    gPrivate.firsttrial = 1;
    gPrivate.macready = 0;
    gPrivate.paradigmdone = 0;
    gPrivate.getnewstim = 0;
    gPrivate.fpon = 0;
	gPrivate.rexready = 0;
}


void rinitf()
{
    wd_disp(D_W_ALLCUR); // all cursors
    wd_src_pos(WIND0, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_check(WIND0, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);

    freeMemory();
    initGlobalParams(); // free memory BEFORE calling this
    updateEyewin(WIND0, gPublic.fp_x, gPublic.fp_y,
            gPublic.fixwin_w, gPublic.fixwin_h);
}

void updateEyewin(int window, int loc_x, int loc_y, int size_x, int size_y)
{
    wd_src_pos(window, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_check(window, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);
    wd_pos(window, loc_x, loc_y);
    wd_siz(window, size_x, size_y);
    wd_cntrl(window, WD_ON);
}

int freeMemory()
{
    // here's where dynamically allocated memory is freed
    return 0;
}

/////////////////
// Handshaking //
/////////////////

int plexreadyandfirsttrial()
{
    return gPrivate.plexready && gPrivate.firsttrial;
}

int macreadyandfpoff()
{
    return gPrivate.macready && !gPrivate.fpon;
}

int macreadyandfpon()
{
    return gPrivate.macready && gPrivate.fpon;
}

int setmacreadytozero()
{
	gPrivate.macready = 0;
	return 0;
}

int setplexreadytozero()
{
	gPrivate.plexready = 0;
	return 0;
}


//////////////////////
// Trial Management //
/////////////////////

int setUpStates()
{
	gPrivate.rexready = 0;
    gPrivate.plexready = 0;
    gPrivate.macready = 0;
	gPrivate.fpon = 0;
	gPrivate.stimon = 0;
	gPrivate.paradigmdone = 0;
	gPrivate.firsttrial = 1;
	gPrivate.getnewstim = 0;
	return 0;
}


int resetTrial()
{
    udp_send_to("AllOff();", "MAC");
    updateEyewin(WINDFIXN, gPublic.fp_x, gPublic.fp_y, gPublic.fixwin_w, gPublic.fixwin_h); // fixation
    updateEyewin(WINDCORR, 1000, 1000, gPublic.targwin_x, gPublic.targwin_y); // correct targ offscreen
    updateEyewin(WINDWRNG, 1000, 1000, gPublic.targwin_x, gPublic.targwin_y); // incorrect targ offscreen
    gPrivate.stimon = 0;
    gPrivate.plexready = 0;
    gPrivate.macready = 0;
	gPrivate.rexready = 0;
	gPrivate.fpon = 0;
	gPrivate.getnewstim = 1;
    dio_off(FTTLBYTE1);
    return 0;
}

int saveTargChoice(int is_correct)
{
    gl_trialparams[TP_CHOOSECORR].val.lval = is_correct;
    return 0;
}

int lopri_clear()
{
    static int prev_load_idx = 0;
    if (prev_load_idx != i_b->pl_lolx && i_b->pl_lodx == i_b->pl_lolx) {
        prev_load_idx = i_b->pl_lolx;
        return 1;
    }
    return 0;
}


///////////////////////////////
// Functions called by menus //
///////////////////////////////

// nothing here at the moment


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
    {"Window width",  &gPublic.fixwin_w,  NP, NP, 0, ME_DEC},
    {"Window height", &gPublic.fixwin_h,  NP, NP, 0, ME_DEC},
    {NS},
};

VLIST menu_stim[] =
{
    {"X Pos",             &gPublic.rf_x,          NP, NP, 0, ME_FLOAT},
    {"Y Pos",             &gPublic.rf_y,          NP, NP, 0, ME_FLOAT},
    {"DVA per Stixel",    &gPublic.dvaperstix,    NP, NP, 0, ME_FLOAT},
    {"Stimulus Duration", &gPublic.stimdur,       NP, NP, 0, ME_FLOAT},
    {NS},
};

MENU umenus[] = {
    {"Fixation", &menu_fp,   NP, NP, 0, NP, NP},
    {"Stimulus", &menu_stim, NP, NP, 0, NP, NP},
    {NS},
};

VLIST state_vl[] = {{NS}};

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
    to closePort
closePort:
    do udpClose()
    to openPort
openPort:
    do udpOpen()
	to stopPlex
stopPlex:
	time 140
	do dio_off(PLXREC)
    to stopSlave
stopSlave:
    do stopSlavePrgm()
    time 50
    to stopOnline
stopOnline:
    do stopOnlinePrgm()
    time 50
    to setUpStates
setUpStates:
	do setUpStates()
	to startSlave
startSlave:
    do startSlavePrgm()
    to startOnline on 1 = gPrivate.macready
startOnline:
    do startOnlinePrgm()
    to initMac on 1 = gPrivate.plexready
initMac:
    do setupMac()
    to getRFxy on 1 = gPrivate.macready
getRFxy:
    do getRFxy()
    to pause0 on 1 = gPrivate.plexready
pause0:
    to finished on 1 = gPrivate.paradigmdone
    to pause1 on +PSTOP & softswitch
    to initTrial
pause1:
    do hideFP()
    to initTrial on -PSTOP & softswitch
initTrial:
    do initTrialParams()
    to finished on 1 = gPrivate.paradigmdone
    to plxStart
plxStart:
    do dio_on(PLXREC)
    time 400
	to plxStart2
plxStart2:
    to getFrameRate on 1 = gPrivate.firsttrial
    to getStimParams on 1 = gPrivate.getnewstim
    to prepareStim
getFrameRate:
    do getFrameRate()
	to getPixPerDeg on 1 % lopri_clear
getPixPerDeg:
    do getPixPerDeg()
    to getBkgndrgb on 1 % lopri_clear
getBkgndrgb:
    do getBkgndrgb()
    to getGamma on 1 % lopri_clear
getGamma:
    do getGamma()
	to getFundamentals on 1 % lopri_clear
getFundamentals:
    do getFundamentals()
    to getMonspd on 1 % lopri_clear
getMonspd:
    do getMonspd()
    to dropHeader on 1 % lopri_clear
dropHeader:
    do dropHeader()
    to getStimParams on 1 % lopri_clear
getStimParams:
    do getStimParams()
	to getX on 1 = gPrivate.plexready
getX:
    do getX()
    to getY on 1 = gPrivate.plexready
getY:
    do getY()
    to sendX on 1 = gPrivate.plexready
sendX:
    do sendX()
    to sendY on 1 = gPrivate.macready
sendY:
    do sendY()
    to prepareStim on 1 = gPrivate.macready
prepareStim:
    do prepareStim() /* sets both plexready and macready to zero */
    to fpOn on 1 = gPrivate.macready
fpOn:
    code FPONCD
    do showFP()
    time 2000
    to fixTime on -WD0_XY & eyeflag
    to error
fixTime:
    code FPACQCD
	time 50
    to error on +WD0_XY & eyeflag
	to displayStim
displayStim:
	code STIMONCD
	do displayStim()
	to error on +WD0_XY & eyeflag
	to preTargDel on 0 = gPrivate.stimon 
preTargDel:
	code STIMOFFCD
	time 100
	rand 500
	to error on +WD0_XY & eyeflag
    to targsOn
targsOn:
	code TARGONCD
	do showTargs()
	time 10
	to fpOff
fpOff:
	code FPOFFCD
	do hideFP()
	to chooseWhich
chooseWhich:
	time 700
    to choseCorrect on -WD1_XY & eyeflag
    to choseWrong on -WD2_XY & eyeflag
    to error
choseCorrect:
	code SACMADCD
	do saveTargChoice(1)
	to evalChoice
choseWrong:
	code SACMADCD
	do saveTargChoice(0)
	to evalChoice
evalChoice:
	to rewOn on -WD1_XY & eyeflag
	to wrongChoice on -WD2_XY & eyeflag
	to error
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
    to dropTrialParams
wrongChoice:
    code ERRCD /* 1014 */
    do resetTrial()
    to dropTrialParams
error:
    code ABORTCD
    do resetTrial()
    time 200
    to dropTrialParams
dropTrialParams:
    do dropTrialParams()
    to markTrialEnd on 1 % lopri_clear
markTrialEnd:
    code EOT
    to plxStop on 1 % lopri_clear
plxStop:
    do dio_off(PLXREC)
    time 140
    to pause0
finished:
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
