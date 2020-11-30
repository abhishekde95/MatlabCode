/*
 * A paradigm for presenting drifting gratings for an initial characterization of recorded neurons.
 *
 * GDLH 3/1/08
 *
 * Adding code to run "Charlie's protocol" (contrast-response functions in L-M and S directions.
 * I'm going to try bailing out of Grating.m on the slave and going into NeuroThresh.m during
 * the experiment.  This will be the first time a paradigm calls more than one slave program.
 *
 * GDLH 4/5/10
 *
 */

#include <stdlib.h>
#include "ldev.h"
#include "labcodes.h"
#include "Grating.h"
#include "UDPtransfer.h"
#include "rigconsts.h"

/***********************************************************************/
#define PI    3.14159
#define ORIENTPROTOCOL 1 /* In the PROTOCOL field from Plexon */
#define SFPROTOCOL 		 2
#define FLASHPROTOCOL	 3
#define COLORPROTOCOL	 4
#define SIZEPROTOCOL	 5
#define MODRATPROTOCOL	 6
#define GABORPROTOCOL	 7

#define OPTOSTIMTYPE 2  /* In the STIMTYTPE field from Plexon */

/********************/
/* Command strings   */
/********************/
const char STARTONLINE_ALERT[] = "GratingOnline;";
const char STARTSLAVE_ALERT1[] = "Grating;";
const char STARTSLAVE_ALERT2[] = "NeuroThresh;";
const char RETURN_ALERT[] = "return;";

/********************/
/* Request strings   */
/********************/
const char RFXY_REQ[] = "sendToRex(udpCom, gl.bar.xy, 'double', message);";
const char NFRAMES_REQ[] = "sendToRex(udpCom, gl.framecounter, 'integer', message);";
const char FRAMERATE_REQ[] = "sendToRex(udpCom, gl.framerate, 'double', message);";
const char GAMMA_REQ[] = "sendToRex(udpCom, gl.cal.gammaTable, 'double', message);";
const char FUNDAMENTALS_REQ[] = "sendToRex(udpCom, gl.cal.fundamentals, 'double', message);";
const char MONSPD_REQ[] = "sendToRex(udpCom, gl.cal.monSpd, 'double', message);";
const char BKGNDRGB_REQ[] = "sendToRex(udpCom, gl.bkgndrgb, 'double', message);";
const char GRATINGPARAMS_REQ[] = "sendToRex(udpCom, gl.gratingparams, 'double', message);";
const char PREFCOLOR1_REQ[] = "sendToRex(udpCom, gl.prefisolumcolors, 'double', message);";
const char PREFCOLOR2_REQ[] = "sendToRex(udpCom, gl.prefccdir, 'double', message);";
/**************************************/
/* Global variable declarations /*
/**************************************/
/* Which fundamentals to use.  Note: single quotes inside double quotes are important */
char gl_fundfilename [] = "'T_cones_smj10.mat'";

int gl_counter = 0;
int gl_fpx = 0;               /* X position of fixation point */
int gl_fpy = 0;               /* Y position of fixation point */
int gl_fpsize = 2;			/* fixation point size in 0.1 deg units */
int gl_eyewinx = 3;     /* Half width of eye window in 0.1 deg units */
int gl_eyewiny = 3;     /* Half height of eye window in 0.1 deg units */
int gl_fprgb[] = {0, 0, 0};
double gl_framerate = 0;
/*int gl_mondist = 100;*/

/* Default grating parameters */
int gl_rfx = 0;
int gl_rfy = 0;
int gl_diam = 100;  /* deg x100 */
int gl_sf = 200;    /* cycles/deg x100 */
int gl_tf = 300;    /* cycles/sec x100 */
int gl_orientation = 0;  /* radians x100 */
int gl_phase = 0;        /* radians x100 */
int gl_conecontrast[] = {40, 40, 40};    /* cone contrast x100 */
int gl_ncycles = 3;    /* number of cycles to show */
int gl_nframesflash = 10;
int gl_stimdur = 100;		/*stimulus duration in 100ths of seconds - ncycles is ignored in GratingOnline.m (YES's version at least) (YES) */
int gl_optoflag =0;        /*binary flag for optical stim. 0 is off; 1 is on (YES) */

int gl_firsttrial = 1;		/* So reset_state command drops the header again */
int gl_switchslaveprgm = 0;   /* Which slave program? 0 = stay, 1 = invoke Grating.m, 2 = invoke NeuroThresh.m */
int gl_currentslaveprgm = 0;   /* Which slave program currently active ? 1 = Grating.m, 2 = NeuroThresh.m */

int gl_done = 0;		/* Have we finished?  Set to '1' when we receive an all '0's message from plexon */
int gl_stimon = 0;
int gl_protocol = 1;  /* Used to force a protocol - don't use to check the current protocol */
int gl_ntrlspercond = 10;  /* Set to '0' to do trials until an ANOVA on the spike counts is significant */   				//YES: was 3, changed to 10//
int gl_getnewparams = 1;
int gl_flash = 0;
int gl_optostim = 0;
int gl_paramSetupComplete = 0;  /*indicates that the trial params have been sent from plexon */

/* Eye position compensation stuff */
int gl_nEPsamplestoaverage = 10;
int gl_EPintegrated = 0;
int gl_EP[2] = {0, 0};

/* Stuff for Charlie's protocol */
float gl_nsigmas = 3.0;
float gl_ncyclesplat = 1.0;
float gl_ncyclesramp = 0.5;
float gl_sigma = 0.4;   /* deg */
int gl_ncontrasts = 7;
float gl_maxcontscalefactor = 2.2;
float gl_mincontscalefactor = 0.25;
int gl_ngaborreps = 8;

/**************************/
/* Stuff for TrialParams array */
/**************************/
/* Indices for the TrialParams */
#define TP_DIAM					0
#define TP_SF						1
#define TP_TF						2
#define TP_LCONT				3
#define TP_MCONT				4
#define TP_SCONT				5
#define TP_ORIENT				6
#define TP_PHASE				7
#define TP_NCYCLES				8
#define TP_STIMTYPE			9
#define TP_PROTOCOL			10
#define TP_NFRAMES			11
#define NTRIALPARAMS  		12

int gl_trialparamcodes[NTRIALPARAMS] = {DIAMCD, SFCD, TFCD,
																LCONTCD, MCONTCD, SCONTCD,
																ORIENTCD, PHASECD, NCYCLESCD,
																 STIMTYPECD, PROTOCOLCD, NFRAMESCD,};
int gl_trialparamtypes[NTRIALPARAMS] = {DOUBLE, DOUBLE, DOUBLE,
																DOUBLE, DOUBLE, DOUBLE,
																DOUBLE, DOUBLE, DOUBLE,
																INT, INT, INT};
ecode_struct gl_trialparams[NTRIALPARAMS];
/* Due to clunky logic, the parameters returned by Plexon in response to sendtorex(gl.gratingparams) have to
 be the first n trial parameters. */

 /*********************/
 /* Function definitions  */
 /*********************/

/*******************************************/
/* Functions for opening/closing UDP channels  */
/*******************************************/

int
udpopen(void)
{
    udp_open_2con(REXIP, MACIP, PLEXIP, UDPPORT);
    return(0);
}

int
udpclose(void)
{
    udp_close();  /* Closing an closed port is better than opening an open one. */
    return(0);
}


/***************************************/
/* Functions for sending requests via UDP  */
/***************************************/

int getFrameRate(void){
	udp_send_to(FRAMERATE_REQ, "MAC");
	return(0);
}

int getNumFrames(void){
	udp_send_to(NFRAMES_REQ, "MAC");
	return(0);
}

int getBkgndrgb(void){
	udp_send_to(BKGNDRGB_REQ, "MAC");
	return(0);
}

int getGammaTable(void){
	udp_send_to(GAMMA_REQ, "MAC");
	return(0);
}

int getFundamentals(void){
	udp_send_to(FUNDAMENTALS_REQ, "MAC");
	return(0);
}

int getMonSpd(void){
	udp_send_to(MONSPD_REQ, "MAC");
	return(0);
}

int getRFxy(void) {
	if (gl_rfx == 0 & gl_rfy == 0)
	{
		udp_send_to(RFXY_REQ, "MAC");
	}
	return(0);
}

int getGratingParams(void){

	if (gl_getnewparams)
	{
		gl_paramSetupComplete = 0;  /* A flag that's set in the listening loop */
		udp_send_to(GRATINGPARAMS_REQ, "PLEXON");
	}
	return(0);
}

int prefColorRequest(void){
	udp_send_to(PREFCOLOR1_REQ, "PLEXON");
	udp_send_to(PREFCOLOR2_REQ, "PLEXON");
	return(0);
}

/* Function called when the user changes the protocol field */
int changeprotocol(void){
	char buf[256];

     if (!gl_firsttrial)
     {
        gl_done = 0;
     	sprintf(buf, "SetUpTrials(%d)", gl_protocol);
		udp_send_to(buf, "PLEXON");
		if (gl_protocol == GABORPROTOCOL & gl_currentslaveprgm == 1)
		{
			gl_switchslaveprgm = 2;
		}
		if (gl_protocol != GABORPROTOCOL & gl_currentslaveprgm == 2)
		{
			gl_switchslaveprgm = 1;
		}
	}
	return(0);
}

int changengratingtrlspercond(void){
	char buf[256];

     if (!gl_firsttrial)
     {
     	sprintf(buf, "gl.ntrlspercond = %d;", gl_ntrlspercond);
		udp_send_to(buf, "PLEXON");
	}
	return(0);
}

int changengabortrlspercond(void){
	char buf[256];

     if (!gl_firsttrial)
     {
     	sprintf(buf, "gl.NGABORREPS = %d;", gl_ngaborreps);
		udp_send_to(buf, "PLEXON");
	}
	return(0);
}
/***************************/
/* Plexon commands via UDP */
/***************************/

int
startonlineprgm(void)
{
     udp_send_to(STARTONLINE_ALERT, "PLEXON");
     return(0);
}

int
endonlineprgm(void)
{
     udp_send_to(RETURN_ALERT, "PLEXON");
     return(0);
}

/*************************/
/* Mac commands via UDP */
/*************************/

int
startslaveprgm1(void)
{
     rigconstants(whichrig());
     udp_send_to(STARTSLAVE_ALERT1, "MAC");
     gl_currentslaveprgm = 1;
   	 gl_switchslaveprgm = 0;
     return(0);
}

int
startslaveprgm2(void)
{

     udp_send_to(STARTSLAVE_ALERT2, "MAC");
     gl_currentslaveprgm = 2;
  	 gl_switchslaveprgm = 0;
     return(0);
}

int
endslaveprgm(void)
{
     udp_send_to(RETURN_ALERT, "MAC");
     return(0);
}

int launchgrating(void)
{
	 endslaveprgm();
	 startslaveprgm1();
	 return(0);
}

int launchneurothresh(void)
{
	 endslaveprgm();
	 startslaveprgm2();
	 return(0);
}

int
preparegabor(void)
{
    char buf[256];

	int nframestot = my_round((gl_ncyclesplat+2*gl_ncyclesramp)/((float)gl_tf/100)*(float)gl_framerate);
	int nframesramp =  my_round((gl_ncyclesramp/((float)gl_tf/100))*(float)gl_framerate);
	int nframesplateau = nframestot-2*nframesramp;
	int gamma = 1;

//	dprintf("%d %d %d\n", nframestot, nframesplateau, nframesramp);

	sprintf(buf, "PrepareGabor(%d, %d, %.7f, %.7f, %d, %.7f, %d, %.7f, %.7f)",
                      nframesplateau, nframesramp,
                      gl_trialparams[TP_ORIENT].val.dval,  gl_trialparams[TP_SF].val.dval,
                      gl_phase, gl_sigma, gamma, (float)gl_tf/100, gl_nsigmas);

    udp_send_to(buf, "MAC");
    showFP();

	return(0);
}

/************************/
/* Getting UDP commands */
/************************/

int
checkmsg()
{
   	char inMsg[REXBUFF];
   	char buf[100];
	static message_struct message;
	int messageAvailable;
	double doublearray[790];
	double tmp;
	int i;

	inMsg[0] = '\0';
	messageAvailable = udp_check(0);

    if (messageAvailable)
    {
         udp_read(inMsg, REXBUFF);
         prepareMsg(&message, inMsg);
         // dprintf("Number of udp_reads: %d\n", ++gl_counter);

         if (!strcmp(message.header, NFRAMES_REQ))
		 {
		    gl_trialparams[TP_NFRAMES].val.lval  = (long)hex2long(message.contents, NULL);
		 }
		  if (!strcmp(message.header, BKGNDRGB_REQ))
		 {
		    hex2double(message.contents, (void *)&doublearray);
		    sendToPlexDataStream("double", (void *)&doublearray, message.size, BKGNDRGBCD, 1);
		 }
		if (!strcmp(message.header, FRAMERATE_REQ))
		{
		    gl_framerate = hex2double(message.contents, NULL);
		    sendToPlexDataStream("double", (void *)&gl_framerate, message.size, FRAMERATECD, 0);
		}

		 if (!strcmp(message.header, GAMMA_REQ))
		 {
		 	  hex2double(message.contents, (void *)&doublearray);
			  sendToPlexDataStream("double", (void *)&doublearray, message.size, GAMMATABLECD, 1);
		 }

		 if (!strcmp(message.header, FUNDAMENTALS_REQ))
		 {
		 	  hex2double(message.contents, (void *)&doublearray);
	 	  	  sendToPlexDataStream("double", (void *)&doublearray, message.size, FUNDAMENTALSCD, 1);
		 }

		 if (!strcmp(message.header, MONSPD_REQ))
		 {
		     hex2double(message.contents, (void *)&doublearray);  //CAH .contents
		     sendToPlexDataStream("double", (void *)&doublearray, message.size, MONSPDCD, 1);
		 }
		  if (!strcmp(message.header, RFXY_REQ))
		 {
		     hex2double(message.contents, (void *)&doublearray);
		     if (message.size && doublearray[0] || doublearray[1])
		     {
					gl_rfx = my_round(doublearray[0]*10);
					gl_rfy = my_round(doublearray[1]*10);
			 }
		 }

		 if (!strcmp(message.header, GRATINGPARAMS_REQ) && (message.size > 0))
		{
			hex2double(message.contents, (void *)&doublearray);
			for (i=0; i < message.size; i++)
			{
				gl_trialparams[i].nelements = 1;
				if (gl_trialparamtypes[i] == DOUBLE)
					gl_trialparams[i].val.dval = doublearray[i];
				if (gl_trialparamtypes[i] == INT)
					gl_trialparams[i].val.lval = (long)doublearray[i];
			}
			if (gl_trialparams[TP_TF].val.dval == 0)
			{
				gl_trialparams[TP_NCYCLES].val.dval	 = (double)gl_nframesflash;
			}
			if (gl_trialparams[TP_PROTOCOL].val.lval == FLASHPROTOCOL)
				gl_flash = 1;
	     	else
				gl_flash = 0;
			if (gl_trialparams[TP_STIMTYPE].val.lval == OPTOSTIMTYPE)
				gl_optostim = 1;
	     	else 
				gl_optostim = 0;
			
			if (gl_trialparams[TP_PROTOCOL].val.lval == GABORPROTOCOL & gl_currentslaveprgm != 2)
			{
				   	gl_switchslaveprgm = 2;
			}

			/* Just for completeness sake */
			if (gl_trialparams[TP_PROTOCOL].val.lval != GABORPROTOCOL & gl_currentslaveprgm != 1)
			{
				   	gl_switchslaveprgm = 1;
			}

			if (gl_trialparams[TP_PROTOCOL].val.lval == 0)  /* if protocol is zero we must be done */
				{
					sprintf(buf, "gl.grating.prefdiam = %.7f; gl.grating.preforient = %.7f; gl.grating.prefsf = %.7f;",
							gl_trialparams[TP_DIAM].val.dval, gl_trialparams[TP_ORIENT].val.dval, gl_trialparams[TP_SF].val.dval);
					udp_send_to(buf,  "MAC");
					gl_done = 1;
				}
		}

		if(!strcmp(message.header, PREFCOLOR1_REQ)){
			hex2double(message.contents, (void *)&doublearray);
			sprintf(buf, "gl.grating.prefIsolum = [%.7f, %.7f, %.7f, %.7f, %.7f, %.7f];",
								doublearray[0], doublearray[1], doublearray[2], doublearray[3], doublearray[4], doublearray[5]);
		    udp_send_to(buf,  "MAC");
		}
		if(!strcmp(message.header, PREFCOLOR2_REQ)){
			hex2double(message.contents, (void *)&doublearray);
			sprintf(buf, "gl.grating.prefccdir = [%.7f, %.7f, %.7f];",
								doublearray[0], doublearray[1], doublearray[2]);
		    udp_send_to(buf,  "MAC");
		}
		if (!strcmp(message.header, "paramSetupComplete"))
			gl_paramSetupComplete = 1;
		if (!strcmp(message.header, "MACSTIMON"))
        	gl_stimon = 1;
        if (!strcmp(message.header, "MACSTIMOFF"))
        	gl_stimon = 0;
    }
     return(0);
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
    ec_send_code_tagged(INITPROTOCOLCD, VALOFFSET+gl_protocol);
    ec_send_code_tagged(NTRLSPERCONDCD, VALOFFSET+gl_ntrlspercond);
    ec_send_code_tagged(INITDIAMCD, VALOFFSET+gl_diam);
    ec_send_code_tagged(INITSFCD, VALOFFSET+gl_sf);
	ec_send_code_tagged(INITTFCD, VALOFFSET+gl_tf);
	ec_send_code_tagged(INITORIENTCD, VALOFFSET+gl_orientation);
	ec_send_code_tagged(INITPHASECD, VALOFFSET+gl_phase);
	ec_send_code_tagged(INITLCONTCD, VALOFFSET+gl_conecontrast[0]);
	ec_send_code_tagged(INITMCONTCD, VALOFFSET+gl_conecontrast[1]);
	ec_send_code_tagged(INITSCONTCD, VALOFFSET+gl_conecontrast[2]);
	ec_send_code_tagged(INITNCYCLESCD, VALOFFSET+gl_ncycles);
	ec_send_code_tagged(NCONTRASTSCD, VALOFFSET+gl_ncontrasts);
	ec_send_code_tagged(NGABORREPSCD, VALOFFSET+gl_ngaborreps);
    /*ec_send_code_tagged(MONDISTCD, VALOFFSET+gl_mondist);*/
	ec_send_code_tagged(STIMDURCD, VALOFFSET+gl_stimdur);
	ec_send_code_tagged(OPTOFLAGCD, VALOFFSET+gl_optoflag);

    sendToPlexDataStream("float", (void *)&gl_nsigmas, 1, NSIGMASCD, 0);
    sendToPlexDataStream("float", (void *)&gl_ncyclesplat, 1, NCYCLESPLATCD, 0);
    sendToPlexDataStream("float", (void *)&gl_ncyclesramp, 1,NCYCLESRAMPCD, 0);
    sendToPlexDataStream("float", (void *)&gl_maxcontscalefactor, 1, MAXCONTSCALEFACTCD, 0);
    sendToPlexDataStream("float", (void *)&gl_mincontscalefactor, 1, MINCONTSCALEFACTCD, 0);

    gl_firsttrial = 0;

    return(0);
}

int
inittrialparams(void)
{
	int i;

	if (gl_getnewparams)
	{
		for (i=0; i < NTRIALPARAMS; i++)
		{
			gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
			gl_trialparams[i].casttype = gl_trialparamtypes[i];
			gl_trialparams[i].nelements = 0;
			gl_trialparams[i].val.dval = 0;
			gl_trialparams[i].val.lpval = NULL;
		}
	}
	gl_EPintegrated = 0;
	gl_EP[0] = 0;
    gl_EP[1] = 0;
    return(0);
}

/* Drop the trial parameters into the data file using sendToPlexDataStream */
/* By convention, trial parameters are never bookended */
int
droptrialparams(void)
{
      int i;
      int ncodes = 0;
      float mspercode = 0.5;   /* Conservative; actual value is closer to 0.33 */

      for (i = 0; i < NTRIALPARAMS; i++)
      {
    		if (gl_trialparams[i].casttype == INT)
    		{
				sendToPlexDataStream("int", (void *)&gl_trialparams[i].val.lval, 1, gl_trialparams[i].flagcodeval, 0);
				ncodes = ncodes + 2;
      		}
     		if (gl_trialparams[i].casttype == DOUBLE)
      		{
     			sendToPlexDataStream("double", (void *)&gl_trialparams[i].val.dval, 1, gl_trialparams[i].flagcodeval, 0);
     			ncodes = ncodes + 9;
     		}
     }

     set_times("droptrialparams", my_ceil(ncodes*mspercode), 0);
     return(0);
}

void
rinitf(void)
{
    /*
     * Initializations.
     */
    udp_send_to("AllOff()", "MAC");
    wd_disp(D_W_ALLCUR);	    /* all cursors */
    wd_src_pos(WIND0, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_check(WIND0, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);
    wd_cntrl(WIND0, WD_ON);
    gl_firsttrial = 1;		/* So reset_state command drops the header again */
	gl_done = 0;
	gl_stimon = 0;
	gl_getnewparams = 1;
	gl_flash = 0;
	gl_optostim = 0;
	gl_EPintegrated = 0;
	gl_EP[0] = 0;
    gl_EP[1] = 0;

	dio_off(FTTL2);
    update_eyewin_big();
    return(0);
}

int
setup(void)
{
     char buf[350];

     /*printf("Remember to move the monitor to %d cm away from eye!\n",gl_mondist);*/
     sprintf(buf, "InitDisplay(%d, %d, %s, %s)",
                      GL_MONDIST, GL_SCREENWIDTH, GL_CALFILENAME, gl_fundfilename);
                
     udp_send_to(buf, "MAC");
     return(0);
}

int
showFP(void)
{
     char buf[256];

     sprintf(buf, "ShowFP(%d, %d, %d, %d, %d, %d)",
     			gl_fpx, gl_fpy, gl_fpsize, gl_fprgb[0], gl_fprgb[1], gl_fprgb[2]);
     udp_send_to(buf, "MAC");
     update_eyewin_big();
     return(0);
}

int
hideFP(void)
{
     char buf[256];

     sprintf(buf, "HideFP()");
     udp_send_to(buf, "MAC");
     return(0);
}

int
showStim(void)
{
     char buf[256];



	if (gl_trialparams[TP_STIMTYPE].val.lval == 1)
	{
     	sprintf(buf, "ShowStim(%d, %d, %.7f, %.7f, %.7f)",
     				gl_rfx, gl_rfy,
     				gl_trialparams[TP_LCONT].val.dval, gl_trialparams[TP_MCONT].val.dval, gl_trialparams[TP_SCONT].val.dval);
	} else {
	     sprintf(buf, "ShowGrating(%d, %d, %.7f, %.7f, %.7f, %.7f, %.7f, %.7f, %.7f, %.7f, %.7f)",
     			gl_rfx, gl_rfy,
     			gl_trialparams[TP_DIAM].val.dval,  gl_trialparams[TP_SF].val.dval, gl_trialparams[TP_TF].val.dval,
     			gl_trialparams[TP_LCONT].val.dval,  gl_trialparams[TP_MCONT].val.dval, gl_trialparams[TP_SCONT].val.dval,
     			gl_trialparams[TP_ORIENT].val.dval,  gl_trialparams[TP_PHASE].val.dval, gl_trialparams[TP_NCYCLES].val.dval);
    }

    udp_send_to(buf, "MAC");
	return(0);
}

int
errortrial(void)
{
    udp_send_to("AllOff()", "MAC");
    gl_stimon = 0;
	gl_getnewparams = 0;
    return(0);
}

int
markgoodtrial(void)
{
     gl_getnewparams = 1;   // Marking a good trial
}

int
update_eyewin_big()
{
    wd_src_pos(WIND0, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_pos(WIND0, gl_fpx, gl_fpy);
    wd_siz(WIND0, 2*gl_eyewinx, 2*gl_eyewiny);
    return(0);
}

int
update_eyewin_small()
{
    wd_src_pos(WIND0, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_pos(WIND0, gl_fpx, gl_fpy);
    wd_siz(WIND0, gl_eyewinx, gl_eyewiny);
    return(0);
}

int
integrateEP(void)
{
	static int counter = 0;

	gl_EP[0] = gl_EP[0]+eyeh;
	gl_EP[1] = gl_EP[1]+eyev;
	counter ++;
	if (counter == gl_nEPsamplestoaverage)
	{
		counter = 0;
		gl_EPintegrated = 1;
	}
}

int
sendEP()
{
	char buf[256];

	sprintf(buf, "eyepos(%f, %f)", (float)gl_EP[0]/gl_nEPsamplestoaverage, (float)gl_EP[1]/gl_nEPsamplestoaverage);
    udp_send_to(buf, "MAC");
}

/*
* User accessible menu
*/
VLIST fp_vl[] =
{
    
    {"Fixation_xpos",	&(gl_fpx),		NP,	NP,	0,	ME_DEC},
    {"Fixation_ypos",	&(gl_fpy),	NP,	NP,	0,	ME_DEC},
    {"Fixation_size",	&(gl_fpsize),		NP,	NP,	0,	ME_DEC},
    {"Fixation_red",	&(gl_fprgb[0]),		NP,	NP,	0,	ME_DEC},
    {"Fixation_green",	&(gl_fprgb[1]),		NP,	NP,	0,	ME_DEC},
    {"Fixation_blue",	&(gl_fprgb[2]),		NP,	NP,	0,	ME_DEC},
    {"Eyewin_width",	&(gl_eyewinx),	NP,	NP,	0,	ME_DEC},
    {"Eyewin_height",	&(gl_eyewiny),	NP,	NP,	0,	ME_DEC},
	{NS},
};

VLIST basic_vl[] =
{
    /*{"Monitor distance",    &(gl_mondist),  NP, NP, 0, ME_DEC},*/
    {"xpos",	&(gl_rfx),		NP,	NP,	0,	ME_DEC},
    {"ypos",	&(gl_rfy),	NP,	NP,	0,	ME_DEC},
	{"protocol", &(gl_protocol), NP, changeprotocol, ME_AFT, ME_DEC},
	{"ntrials per cond", &(gl_ntrlspercond), NP, changengratingtrlspercond, ME_AFT, ME_DEC},
	{"stimulus duration", &(gl_stimdur), NP, NP, 0, ME_DEC},
	{"optical stim on", &(gl_optoflag), NP, NP, 0, ME_DEC},
	{NS},
};

VLIST grating_vl[] =
{
    {"diam",	&(gl_diam),		NP,	NP,	0,	ME_DEC},
    {"sf",	&(gl_sf),	NP,	NP,	0,	ME_DEC},
	{"tf", &(gl_tf), NP, NP, 0, ME_DEC},
	{"orient", &(gl_orientation), NP, NP, 0, ME_DEC},
	{"phase", &(gl_phase), NP, NP, 0, ME_DEC},
	{"L-cone contrast", &(gl_conecontrast[0]), NP, NP, 0, ME_DEC},
	{"M-cone contrast", &(gl_conecontrast[1]), NP, NP, 0, ME_DEC},
	{"S-cone contrast", &(gl_conecontrast[2]), NP, NP, 0, ME_DEC},
	{"ncycles", &(gl_ncycles), NP, NP, 0, ME_DEC},
	{"nflashframes", &(gl_nframesflash), NP, NP, 0, ME_DEC},
	{NS},
};

VLIST gabor_vl[] =
{
    {"sigma",	&(gl_sigma),		NP,	NP,	0,	ME_FLOAT},
    {"nsigmas",	&(gl_nsigmas),		NP,	NP,	0,	ME_FLOAT},
    {"ncycles_plateau",	&(gl_ncyclesplat),		NP,	NP,	0,	ME_FLOAT},
    {"ncycles_ramp",	&(gl_ncyclesramp),		NP,	NP,	0,	ME_FLOAT},
    {"ncontrasts",	&(gl_ncontrasts),		NP,	NP,	0,	ME_DEC},
    {"nrepeats",	&(gl_ngaborreps),		NP,	changengabortrlspercond,	ME_AFT,	ME_DEC},
    {"max_contrast_scale_factor",	&(gl_maxcontscalefactor),		NP,	NP,	0,	ME_FLOAT},
    {"min_contrast_scale_factor",	&(gl_mincontscalefactor),		NP,	NP,	0,	ME_FLOAT},
	{NS},
};

MENU umenus[] =
{
	{"FixpointStuff",	&fp_vl,	NP,	NP,	0,	NP,	NP},
	{"BasicStuff",	&basic_vl,	NP,	NP,	0,	NP,	NP},
	{"GratingDefaults",	&grating_vl,	NP,	NP,	0,	NP,	NP},
	{"GaborDefaults",	&gabor_vl,	NP,	NP,	0,	NP,	NP},

	{NS}
};

/*
 * User-supplied real-time variable table
 */
RTVAR rtvars[] = {
};

%%
id 700
restart rinitf
main_set {
status ON
begin	first:
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
    	time 50
    	to startslave
     startslave:
     	do startslaveprgm1()
  	    to startonline
  	startonline:
  	 	do startonlineprgm()
     	time 500
     	to getRFxy
     getRFxy:
    	do getRFxy()
    	time 10
        to initmac
     initmac:
        do setup()
        to pause0
	pause0:
		to pause1 on +PSTOP & softswitch
		to inittrial
	pause1:
		to inittrial on -PSTOP & softswitch
    inittrial:
		do inittrialparams()
		to finished on 1 = gl_done
		to plxstart
    plxstart:
        do dio_on(PLXREC)
        time 400
        to drophdr on 1 = gl_firsttrial
        to getGratingParams
	drophdr:
        do drophdr()
        time 300
        to getFrameRate
    getFrameRate:
    	do getFrameRate()
    	time 10
		to getbkgndrgb
	getbkgndrgb:
		do getBkgndrgb()
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
        time 1500
     	to getGratingParams
	getGratingParams:
     	do getGratingParams()
		to requestPrefColor on 1= gl_done
		to fpon on 1 = gl_paramSetupComplete
	fpon:
		code FPONCD
		do showFP()
		time 2000
		to prestimdel on -WD0_XY & eyeflag
		to error
	prestimdel:
		code FPACQCD
		time 300
		to error on +WD0_XY & eyeflag
		to finished on 1 = gl_done
    	to branch
    branch:
    	to shrinkeyewin on 1 = gl_flash
    	to launchneurothresh on 2 = gl_switchslaveprgm
    	to launchgrating on 1 = gl_switchslaveprgm
    	to showStim
    shrinkeyewin:
    	do update_eyewin_small()
		time 50
		to error on +WD0_XY & eyeflag
		to epcomp
    epcomp:
    	do integrateEP()
    	to error on +WD0_XY & eyeflag
    	to epcomp2
	epcomp2:
    	to sendEP on 1 = gl_EPintegrated
    	to epcomp
    sendEP:
    	do sendEP()
    	to showStim
   launchgrating:
    	do launchgrating()
    	time 500
    	to showStim
    launchneurothresh:
    	do launchneurothresh()
    	time 200
    	to preparegabor
    preparegabor:
    	do preparegabor()
    	time 50
    	to showStim
	showStim:
		do showStim()
		to stimcode on 1 = gl_stimon   /* Mac has sent the "stimulus on" signal */
	stimcode:
		code STIMONCD
		to laseron on 1 = gl_optostim
		to error on +WD0_XY & eyeflag
		to laseroff on 0 = gl_stimon
	laseron:
		do dio_on(FTTL2)
		to error on +WD0_XY & eyeflag
		to laseroff on 0 = gl_stimon
   laseroff:
		do dio_off(FTTL2)
		to poststimdelay   
   poststimdelay:
   		code STIMOFFCD
   		/* This is the indication to GratingOnline that a successful trial was completed */
   		/* There should be no way to error after this */
   		time 50
   		random 50
		to fpoff
	fpoff:
		code FPOFFCD
		do hideFP()
		to rewon
	rewon:
		code REWCD
		do dio_on(REW)
		time 50
		to rewoff
	rewoff:
		do dio_off(REW)
		time 200							/* YES: added pause before extra reward state, and commented out to getnframes */
		to rewonextra
		/* to getnframes */				
	rewonextra:			
		do dio_on(REW)
		time 50
		to rewoffextra
	rewoffextra:
		do dio_off(REW)
		to getnframes	
	getnframes:
		do getNumFrames()
		time 10
		to markgoodtrial
	markgoodtrial:
		do markgoodtrial()
		to droptrialparams
	error:
		code ABORTCD
		do errortrial()
		to droptrialparams
	droptrialparams:
		do droptrialparams()
		to marktrialend
	marktrialend:
		code EOT
		time 10
		to plxstop
	requestPrefColor:
		do prefColorRequest()
		time 100
		to plxstop
	plxstop:
        do dio_off(PLXREC)
		time 140
	 	to timeout on 0 = gl_getnewparams
		to pause0
	timeout:
		time 1000
		to pause0
	finished:
	    time 1000
		to pause0
}

msg_set {
status ON
   begin
   msgFirst:
       to msgSecond
   msgSecond:
       do checkmsg();
       to msgSecond
abort list:
}

