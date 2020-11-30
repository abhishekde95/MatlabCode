/*
 * This paradigm finds the subunits of a RF, and stimulates each subunit
 * with different combinations of L- and M-cone contrast.  Version 2.0
 */

#include <stdio.h>
#include "ldev.h"
#include "labcodes.h"
#include "UDPtransfer.h"
#include "rigconsts.h"
#include "GridLMSubunitDN.h"

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
#define TP_NFRAMES		6
#define TP_EPOCH		7
#define TP_SEED			8
#define NTRIALPARAMS    9
#define PI 3.141592653589793

int gl_trialparamcodes[NTRIALPARAMS] = {LCCCD, MCCCD, SCCCD, NSTIXGRIDCD, DVAPERSTIXCD, STIMDURCD, NFRAMESCD, EPOCHCD, SEEDCD};
int gl_trialparamtypes[NTRIALPARAMS] =  {FLOAT, FLOAT, FLOAT, INT, FLOAT, FLOAT, INT, INT, LONG};
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
const char NFRAMES_REQ[] = "sendToRex(udpCom, gl.stim.framecounter, 'integer', message);";
const char SEED_REQ[] = "sendToRex(udpCom, gl.stim.seed, 'long', message);";


#define DECL_UDPSEND_FUNC( funcname, reqstr, dest ) int funcname##() { udp_send_to(reqstr, dest); return 0; }
DECL_UDPSEND_FUNC( startSlavePrgm, "GridLMSubunitDN_Slave();", "MAC" );
DECL_UDPSEND_FUNC( stopSlavePrgm, "return;", "MAC" );
DECL_UDPSEND_FUNC( macAllOff, "AllOff();", "MAC" );
DECL_UDPSEND_FUNC( startOnlinePrgm, "GridLMSubunitDNOnline();", "PLEXON" );
DECL_UDPSEND_FUNC( stopOnlinePrgm, "return;", "PLEXON" );
DECL_UDPSEND_FUNC( getFrameRate, FRAMERATE_REQ, "MAC" );
DECL_UDPSEND_FUNC( getPixPerDeg, PIXPERDEG_REQ, "MAC" );
DECL_UDPSEND_FUNC( getBkgndrgb, BKGNDRGB_REQ, "MAC" );
DECL_UDPSEND_FUNC( getGamma,GAMMA_REQ, "MAC" );
DECL_UDPSEND_FUNC( getFundamentals,FUNDAMENTALS_REQ, "MAC" );
DECL_UDPSEND_FUNC( getMonspd, MONSPD_REQ, "MAC" );

#undef DECL_UDPSEND_FUNC

///////////////////////////
//  Mac commands via UDP //
///////////////////////////

int getRFxy()
{
    if ( gPublic.rf_x == 0 && gPublic.rf_y == 0 ) udp_send_to(RFXY_REQ, "MAC");
    dprintf("Starting anew:\n\n\n");
    return 0;
}

int setupMac()
{
    char buf[256];
    sprintf(buf, "InitDisplay(%d, %d, %s, %s)",
        GL_MONDIST, GL_SCREENWIDTH, GL_CALFILENAME, FUNDFILENAME);
    udp_send_to(buf, "MAC");
    return 0;
}

//////////////////////////
// Getting UDP commands //
//////////////////////////

int getStimParams()
{
    char buf[256];
    gPrivate.plexdone = 0;
    sprintf(buf, TRIALPARAMS_REQ);
    udp_send_to(buf,"PLEXON");
    return 0;
}

int getX()
{
    char buf[256];

    sprintf(buf, XCOORDINATES_REQ);
    gPrivate.plexdone = 0;
    udp_send_to(buf,"PLEXON");
	
    return 0;
}

int getY()
{
    char buf[256];

    sprintf(buf, YCOORDINATES_REQ);
    gPrivate.plexdone = 0;
    udp_send_to(buf,"PLEXON");

    return 0;
}

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
		
		//dprintf("%s\n", message.header);

        if ( !strcmp(message.header, "macdone") ) {
            gPrivate.macdone = 1;
        } else if ( !strcmp(message.header, "plexdone") ) {
            gPrivate.plexdone = 1;
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
			gPublic.seed = (long) doublearray[6];
			gPrivate.epoch = (int) doublearray[7];
			sprintf(buf, "%f %f %f %f %f %d %ld %d\n", gPublic.stim_lmscc[0], gPublic.stim_lmscc[1], gPublic.stim_lmscc[2], gPublic.dvaperstix, gPublic.stimdur, gPublic.nstixgrid, gPublic.seed, gPrivate.epoch);
			//dprintf("%s",buf);
            gPrivate.plexdone = 1;
        } else if ( !strcmp(message.header, XCOORDINATES_REQ)) {
            hex2long(message.contents, (void *)&gPublic.grid_x);
            gPrivate.grid_len = message.size;
            gPrivate.plexdone = 1;
        } else if ( !strcmp(message.header, YCOORDINATES_REQ)) {
            hex2long(message.contents, (void *)&gPublic.grid_y);
            gPrivate.getnewstim = 0;
            gPrivate.plexdone = 1;
        } else if ( !strcmp(message.header, NFRAMES_REQ)) {
            gPublic.nframes = hex2long(message.contents, NULL);
			sprintf(buf, "%d\n", gPublic.nframes);
			//dprintf(buf);
		//} else if ( !strcmp(message.header, SEED_REQ)) {
		//	gPublic.seed = hex2long(message.contents, NULL);
		//	sprintf(buf, "%ld\n", gPublic.seed);
			//printf(buf);
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
    ec_send_code_tagged(EYEWINXCD       , VALOFFSET + gPublic.eyewin_x);
    ec_send_code_tagged(EYEWINYCD       , VALOFFSET + gPublic.eyewin_y);
	
    sendToPlexDataStream("float",(void *)&gPublic.rf_x, 1,RFXCD, 0);
    sendToPlexDataStream("float",(void *)&gPublic.rf_y, 1,RFYCD, 0);

	ec_send_code_tagged(HEADERCMPLTCD       , VALOFFSET + 1);
	
    gPrivate.firsttrial = 0;

    return 0;
}

// currently just sets the dummy trial parameter to 1
int initTrialParams()
{
    int i;

    // we need to make sure there aren't any buffer over- or underruns
/*  for (i = 0; i < 600; ++i) {
        gPublic.grid_x[i] = 0;
        gPublic.grid_y[i] = 0;
    } */

    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
        gl_trialparams[i].casttype = gl_trialparamtypes[i];
        gl_trialparams[i].nelements = 0;
        gl_trialparams[i].val.dval = 0;
        gl_trialparams[i].val.lpval = NULL;
        /* switch (gl_trialparamtypes[i]) {
            case DOUBLE:
                gl_trialparams[i].val.dval = 1;
                break;
            case FLOAT:
                gl_trialparams[i].val.fval = 1;
                break;
            case INT:
                gl_trialparams[i].val.lval = 1;
                break;
        }
 */
    }
    return 0;
}

int dropTrialParams()
{
    int i, ncodes = 0;
    float mspercode = 0.7; // Conservative; actual value is closer to 0.33
    //dprintf("Dropping Trial Parameters\n");

    gPrivate.plexdone = 0;

    // no better place to put these assignments
    gl_trialparams[TP_LCC].val.fval = gPublic.stim_lmscc[0];
    gl_trialparams[TP_MCC].val.fval = gPublic.stim_lmscc[1];
    gl_trialparams[TP_SCC].val.fval = gPublic.stim_lmscc[2];
    gl_trialparams[TP_NSTIXGRID].val.lval = gPublic.nstixgrid;
    gl_trialparams[TP_DVAPERSTIX].val.fval = gPublic.dvaperstix;
    gl_trialparams[TP_STIMDUR].val.fval = gPublic.stimdur;
	gl_trialparams[TP_NFRAMES].val.lval = gPublic.nframes;
	gl_trialparams[TP_SEED].val.lval = gPublic.seed;
	gl_trialparams[TP_EPOCH].val.lval = gPrivate.epoch;

    sendToPlexDataStream("int", (void*) &gPublic.grid_x, gPrivate.grid_len, GRIDXCD, 1);
    sendToPlexDataStream("int", (void*) &gPublic.grid_y, gPrivate.grid_len, GRIDYCD, 1);
	//sendToPlexDataStream("float", (void*) &gPrivate.bkgndrgb, 3, BKGNDRGBCD, 1);
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
        gPublic.eyewin_x = 6; gPublic.eyewin_y = 6;
        gPublic.fp_size = 2;
        gPublic.stim_lmscc[0] = 0; gPublic.stim_lmscc[1] = 0; gPublic.stim_lmscc[2] = 0;
        gPublic.nstixgrid = 0; // In DVA
        gPublic.dvaperstix = .2;
        gPublic.stimdur = 1; // In seconds
		gPublic.rf_x = 0; gPublic.rf_y = 0;
		gPublic.nframes = 0;
		gPublic.seed = 0;
		gPrivate.epoch = 0;

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
    //dprintf("In initGlobalParams: gPrivate.firsttrial = %d\n",gPrivate.firsttrial);
    gPrivate.macdone = 0;
    gPrivate.paradigmdone = 0;
    gPrivate.getnewstim = 0;
    gPrivate.fpoff = 1;
    gPrivate.fpon = 0;
}


void rinitf()
{
    wd_disp(D_W_ALLCUR); // all cursors
    wd_src_pos(WIND0, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_check(WIND0, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);

    freeMemory();
    initGlobalParams(); // free memory BEFORE calling this
    updateEyewin(WIND0);
    //dprintf("In rinitf; gPrivate.firsttrial = %d\n",gPrivate.firsttrial);
    //srand(time(NULL));
}

void updateEWHelper(int window, int loc_x, int loc_y, int size_x, int size_y)
{
    wd_src_pos(window, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_pos(window, loc_x, loc_y);
    wd_siz(window, size_x, size_y);
    wd_cntrl(window, WD_ON);
}

void updateEyewin(int window)
{
    if ( window == WIND0 ) {
        updateEWHelper(window, gPublic.fp_x, gPublic.fp_y,
                gPublic.eyewin_x, gPublic.eyewin_y);
    }

}

/////////////////////////////////////////
// Functions specific to this paradigm //
/////////////////////////////////////////

int showFP()
{
    char buf[256];
    gPrivate.macdone = 0; // This resets the "handshake" bw Mac and Rex that occurs in PrepareStim
    sprintf(buf, "ShowFP(%d, %d, %d, %d, %d, %d)",
        gPublic.fp_x, gPublic.fp_y, gPublic.fp_size,
        gPublic.fp_rgb[0], gPublic.fp_rgb[1], gPublic.fp_rgb[2]);
    udp_send_to(buf, "MAC");
    updateEyewin(WIND0);
    gPrivate.fpon = 1;
    gPrivate.fpoff = 0;

    return 0;
}

int hideFP()
{
    char buf[256];
    gPrivate.macdone = 0;
    gPrivate.plexdone = 0;
    sprintf(buf, "HideFP()");
    udp_send_to(buf, "MAC");
    gPrivate.fpon = 0;
    gPrivate.fpoff = 1;

    return 0;
}

int plexdoneandfirsttrial()
{
    return gPrivate.plexdone && gPrivate.firsttrial;
}

int macdoneandfpoff()
{
    return gPrivate.macdone && gPrivate.fpoff;
}

int macdoneandfpon()
{
    return gPrivate.macdone && gPrivate.fpon;
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

int sendX()
{
    char argsbuf[1024], reqbuf[1100];
	gPrivate.macdone = 0;
    sprintf(reqbuf, "sendX([%s])", format_grid_request(argsbuf, gPublic.grid_x));
    udp_send_to(reqbuf, "MAC");
    //dprintf("Arg into sendX: [%s]\n", argsbuf);
    return 0;
}

int sendY()
{
    char argsbuf[1024], reqbuf[1100];
    gPrivate.macdone = 0;
	sprintf(reqbuf, "sendY([%s])", format_grid_request(argsbuf, gPublic.grid_y));
    udp_send_to(reqbuf, "MAC");
    //dprintf("Arg into sendY: [%s]\n", argsbuf);
    return 0;
}

int prepareStim()
{
    char buf[256];
    gPrivate.macdone = 0;
    gPrivate.plexdone = 0;

    sprintf(buf, "PrepareStim(%d, %f, %f, %.15f, %.15f, %.15f, %d, %d)",
        gPublic.nstixgrid, gPublic.dvaperstix, gPublic.stimdur,
        gPublic.stim_lmscc[0], gPublic.stim_lmscc[1], gPublic.stim_lmscc[2],
		gPublic.seed, gPrivate.epoch);

    //dprintf("Arg into prepareStim: %s\n", buf);
		
    udp_send_to(buf, "MAC");

    gPrivate.firsttrial = 0;
    return 0;
}

int displayStim()
{
    char buf[256];
    gPrivate.macdone = 0;
    sprintf(buf, "displayStim()");
    udp_send_to(buf, "MAC");
    return 0;
}

int errorTrial()
{
    macAllOff();
    gPrivate.macdone = 0;
    gPrivate.plexdone = 0;
    gPrivate.fpon = 0;
    gPrivate.fpoff = 1;
    return 0;
}


int correctTrial()
{
    // here's where you increment the trial counter and
    // check if you've done enough trials for that block
    // or if you're done with the experiment
	if ( gPrivate.epoch == 1 ) {
        macAllOff();
    }
	gPrivate.macdone = 0;
	gPrivate.plexdone = 0;
	gPrivate.fpon = 0;
	gPrivate.fpoff = 1;
	gPrivate.getnewstim = 1;

	
    return 0;
}

int freeMemory()
{
    // here's where dynamically allocated memory is freed
    return 0;
}


int getNumFrames()
{
    udp_send_to(NFRAMES_REQ, "MAC");
    return 0;
}

int getSeed()
{
    udp_send_to(SEED_REQ, "MAC");
    return 0;
}

int setMacDoneToZero()
{
	gPrivate.macdone = 0;
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
    {"Window width",  &gPublic.eyewin_x,  NP, NP, 0, ME_DEC},
    {"Window height", &gPublic.eyewin_y,  NP, NP, 0, ME_DEC},
    {NS},
};

VLIST menu_stim[] =
{
    {"X Pos",             &gPublic.rf_x,          NP, NP, 0, ME_FLOAT},
    {"Y Pos",             &gPublic.rf_y,          NP, NP, 0, ME_FLOAT},
    {"DVA per Stixel",    &gPublic.dvaperstix,    NP, NP, 0, ME_FLOAT},
    {"Stimulus Duration", &gPublic.stimdur,       NP, NP, 0, ME_FLOAT},
    {"L(cc)",             &gPublic.stim_lmscc[0], NP, NP, 0, ME_FLOAT},
    {"M(cc)",             &gPublic.stim_lmscc[1], NP, NP, 0, ME_FLOAT},
    {"S(cc)",             &gPublic.stim_lmscc[2], NP, NP, 0, ME_FLOAT},

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
    to endOnline
endOnline:
    do stopOnlinePrgm()
    time 50
    to startSlave
startSlave:
    do startSlavePrgm()
    time 50
    to startOnline
startOnline:
    do startOnlinePrgm()
    time 1000
    to initMac
initMac:
    do setupMac()
    time 1000
    to getRFxy
getRFxy:
    do getRFxy()
    time 10
    to pause0
pause0:
    to finished on 1 = gPrivate.paradigmdone
    to pause1 on +PSTOP & softswitch
    to initTrial
pause1:
    do hideFP()
    to initTrial on -PSTOP & softswitch
initTrial:
    do initTrialParams()
    time 200
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
    time 20
    to getPixPerDeg
getPixPerDeg:
    do getPixPerDeg()
    time 20
    to getBkgndrgb
getBkgndrgb:
    do getBkgndrgb()
    time 10
    to getGamma
getGamma:
    do getGamma()
    time 2700
	to getFundamentals
getFundamentals:
    do getFundamentals()
    time 600
    to getMonspd
getMonspd:
    do getMonspd()
    time 1800
    to dropHeader
dropHeader:
    do dropHeader()
    time 300
    to getStimParams
getStimParams:
    do getStimParams()
	to getX on 1 = gPrivate.plexdone
getX:
    do getX()
    to getY on 1 = gPrivate.plexdone
getY:
    do getY()
    to sendX on 1 = gPrivate.plexdone
sendX:
    do sendX()
    to sendY on 1 = gPrivate.macdone
sendY:
    do sendY()
    to prepareStim on 1 = gPrivate.macdone
prepareStim:
    do prepareStim() /* sets plexdone to zero */
    to fpOn on 1 % macdoneandfpoff
    to fixTime on 1 % macdoneandfpon
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
	to waitforbsln
waitforbsln:
    to displayStim1 on 1 = gPrivate.epoch
	to displayStim2 on 2 = gPrivate.epoch
displayStim1:
	code STIMONCD
	time 100
	do displayStim()
	to error on +WD0_XY & eyeflag
	to stimPause
stimPause:
	time 1250
    to dropStimOffCd on +WD0_XY & eyeflag
	to rewOn
dropStimOffCd:
	code STIMOFFCD	
	time 10
	to correctTrial
displayStim2:
    do displayStim()
	to dropStimOnCode on 1 = gPrivate.macdone
	to error on +WD0_XY & eyeflag
dropStimOnCode:
	code STIMONCD
	do setMacDoneToZero()
	to waitForStimOff
waitForStimOff:
    to poststimdel on 1 = gPrivate.macdone
    to error on +WD0_XY & eyeflag
poststimdel:
    code STIMOFFCD
    time 200
    to error on +WD0_XY & eyeflag
    to rewOn
rewOn:
    code REWCD
    do dio_on(REW)
    time 30
    to rewOff
rewOff:
    do dio_off(REW)	
	to stimPause on 1 = gPrivate.epoch    
	to correctTrial on 2 = gPrivate.epoch
error:
    code ABORTCD
    do errorTrial()
    time 200
    to dropTrialParams
correctTrial:
    do correctTrial()
    to getnframes
getnframes:
    do getNumFrames()
    time 20 
    to dropTrialParams
dropTrialParams:
    do dropTrialParams()
    to markTrialEnd
markTrialEnd:
    code EOT
    time 20
    to waitForFR
waitForFR:
	time 15000
    to plxStop on 1 = gPrivate.plexdone
    to plxStop
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
