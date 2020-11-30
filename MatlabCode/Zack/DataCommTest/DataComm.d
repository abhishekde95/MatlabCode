/*
 * DataComm: an all inclusive test of the MATLAB <-> REX UDP framework
 */

#include <stdlib.h>
#include "ldev.h"
#include "labcodes.h"
#include "notUDPtransfer.h"
#include "rigconsts.h"
#include "DataComm.h"

// "public" globals
int gl_fp_x = 0;
int gl_fp_y = 0;
int gl_rf_x = 0;
int gl_rf_y = 0;
int gl_eyewin_x  = 8;
int gl_eyewin_y = 8;
int gl_targwin_x = 13;
int gl_targwin_y = 13;
int gl_fp_size = 2;
int gl_targ_size = 4;
int gl_fp_rgb[] = {0, 0, 0};
int gl_targ_rgb[] = {0, 0, 0};

// "private" globals
int gl_done = 0;
int gl_trial_count = 0;
int gl_firsttrial = 1;
int gl_mac_ready = 0;

/*
 * Stuff for TrialParams array
 */

// Indices for gl_trialparams
#define NTRIALPARAMS    0

int gl_trialparamcodes[NTRIALPARAMS] = {};
int gl_trialparamtypes[NTRIALPARAMS] = {};
ecode_struct gl_trialparams[NTRIALPARAMS];

/*
 * Functions for opening/closing UDP channels
 */

int udpOpen()
{
    rigconstants(whichrig());
    udp_open_2con(REXIP, MACIP, PLEXIP, UDPPORT);
    return 0;
}

int udpClose()
{
    udp_close();
    return 0;
}

/*
 * Functions for sending requests via UDP
 */
 
#define DECL_S2R_STR( name, var, type ) const char name[] = "notsendToRex(udpCom, " #var ", '" #type "', message);";
// header-bound variables of 1 element
DECL_S2R_STR( RFXY_REQ, gl.bar.xy, double )
DECL_S2R_STR( HDI_REQ, gl.hd.i, int )
DECL_S2R_STR( HDD_REQ, gl.hd.d, double )
DECL_S2R_STR( HDF_REQ, gl.hd.f, float )
DECL_S2R_STR( HDL_REQ, gl.hd.l, long )
DECL_S2R_STR( HDC_REQ, gl.hd.c, char )
DECL_S2R_STR( HDS_REQ, gl.hd.s, short )
DECL_S2R_STR( HDUS_REQ, gl.hd.us, ushort )
DECL_S2R_STR( HDUL_REQ, gl.hd.ul, ulong )
DECL_S2R_STR( HDUI_REQ, gl.hd.ui, uint )
DECL_S2R_STR( HDLL_REQ, gl.hd.ll, llong )
DECL_S2R_STR( HDULL_REQ, gl.hd.ull, ullong )
// header-bound array variables
DECL_S2R_STR( HDAI_REQ, gl.hd.ai, int )
DECL_S2R_STR( HDAD_REQ, gl.hd.ad, double )
DECL_S2R_STR( HDAF_REQ, gl.hd.af, float )
DECL_S2R_STR( HDAL_REQ, gl.hd.al, long )
DECL_S2R_STR( HDAC_REQ, gl.hd.ac, char )
DECL_S2R_STR( HDAS_REQ, gl.hd.as, short )
DECL_S2R_STR( HDAUS_REQ, gl.hd.aus, ushort )
DECL_S2R_STR( HDAUL_REQ, gl.hd.aul, ulong )
DECL_S2R_STR( HDAUI_REQ, gl.hd.aui, uint )
DECL_S2R_STR( HDALL_REQ, gl.hd.all, llong )
DECL_S2R_STR( HDAULL_REQ, gl.hd.aull, ullong )
#undef DECL_S2R_STR

#define DECL_UDPSEND_FUNC( funcname, reqstr, dest ) int funcname##() { udp_send_to(reqstr, dest); return 0; }
DECL_UDPSEND_FUNC( startSlavePrgm, "DataCommSlave();", "MAC" )
DECL_UDPSEND_FUNC( stopSlavePrgm, "return;", "MAC" )
DECL_UDPSEND_FUNC( macAllOff, "AllOff();", "MAC" )
DECL_UDPSEND_FUNC( hideFP, "HideFP();", "MAC" )
// get the header-bound variables of 1 element
DECL_UDPSEND_FUNC( getHDi, HDI_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDd, HDD_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDf, HDF_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDl, HDL_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDc, HDC_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDs, HDS_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDus, HDUS_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDul, HDUL_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDui, HDUI_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDll, HDLL_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDull, HDULL_REQ, "MAC" )
// get the header-bound array variables
DECL_UDPSEND_FUNC( getHDai, HDAI_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDad, HDAD_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDaf, HDAF_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDal, HDAL_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDac, HDAC_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDas, HDAS_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDaus, HDAUS_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDaul, HDAUL_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDaui, HDAUI_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDall, HDALL_REQ, "MAC" )
DECL_UDPSEND_FUNC( getHDaull, HDAULL_REQ, "MAC" )
#undef DECL_UDPSEND_FUNC

/*
 * Mac commands via UDP
 */

int getRFxy()
{
    if ( gl_rf_x == 0 && gl_rf_y == 0 ) {
        udp_send_to(RFXY_REQ, "MAC");
    }
    return 0;
}

int setupMac()
{
    char buf[256];
    gl_mac_ready = 0;
    sprintf(buf, "InitDisplay(%d, %d, %s)", GL_MONDIST, GL_SCREENWIDTH, GL_CALFILENAME);
    udp_send_to(buf, "MAC");
    return 0;
}

int initHeaderVars()
{
    gl_mac_ready = 0;
    udp_send_to("InitHeaderVars()", "MAC");
    return 0;
}

int initTrialVars()
{
    gl_mac_ready = 0;
    udp_send_to("InitTrialVars()", "MAC");
    return 0;
}

/*
 * Getting UDP commands
 */

int checkMsg()
{
    char inMsg[REXBUFF];
    static message_struct message;
    int messageAvailable;
    double dvals[2]; // adjust this if you need more storage
    int inttmp[20];
    double doubletmp[768];
    float floattmp[722];
    long longtmp[72];
    char chartmp[32];
    short shorttmp[52];
    unsigned short ushorttmp[123];
    unsigned long ulongtmp[191];
    unsigned int uinttmp[222];
    long long llongtmp[435];
    unsigned long long ullongtmp[511];    

    inMsg[0] = '\0';
    messageAvailable = udp_check(0);

    if ( messageAvailable ) {
        udp_read(inMsg, REXBUFF);
        prepareMsg(&message, inMsg);

        if ( !strcmp(message.header, "DISPLAYINIT") ) {
            gl_mac_ready = 1;
        } else if ( !strcmp(message.header, RFXY_REQ) ) {
            hex2double(message.contents, (void *) &dvals);
            if ( message.size && (dvals[0] || dvals[1]) ) {
                gl_rf_x = my_round(10. * dvals[0]);
                gl_rf_y = my_round(10. * dvals[1]);
            }
        } else if ( !strcmp(message.header, HDI_REQ) ) {
            SEND_VAL_TO_STREAM(int, HDICD);
        } else if ( !strcmp(message.header, HDD_REQ) ) {
            SEND_VAL_TO_STREAM(double, HDDCD);
        } else if ( !strcmp(message.header, HDF_REQ) ) {
            SEND_VAL_TO_STREAM(float, HDFCD);
        } else if ( !strcmp(message.header, HDL_REQ) ) {
            SEND_VAL_TO_STREAM(long, HDLCD);
        } else if ( !strcmp(message.header, HDC_REQ) ) {
            SEND_VAL_TO_STREAM(char, HDCCD);
        } else if ( !strcmp(message.header, HDS_REQ) ) {
            SEND_VAL_TO_STREAM(short, HDSCD);
        } else if ( !strcmp(message.header, HDUS_REQ) ) {
            SEND_VAL_TO_STREAM(ushort, HDUSCD);
        } else if ( !strcmp(message.header, HDUL_REQ) ) {
            SEND_VAL_TO_STREAM(ulong, HDULCD);
        } else if ( !strcmp(message.header, HDUI_REQ) ) {
            SEND_VAL_TO_STREAM(uint, HDUICD);
        } else if ( !strcmp(message.header, HDLL_REQ) ) {
            SEND_VAL_TO_STREAM(llong, HDLLCD);
        } else if ( !strcmp(message.header, HDULL_REQ) ) {
            SEND_VAL_TO_STREAM(ullong, HDULLCD);
        } else if ( !strcmp(message.header, HDAI_REQ) ) {
            SEND_ARY_TO_STREAM(int, HDAICD);
        } else if ( !strcmp(message.header, HDAD_REQ) ) {
            SEND_ARY_TO_STREAM(double, HDADCD);
        } else if ( !strcmp(message.header, HDAF_REQ) ) {
            SEND_ARY_TO_STREAM(float, HDAFCD);
        } else if ( !strcmp(message.header, HDAL_REQ) ) {
            SEND_ARY_TO_STREAM(long, HDALCD);
        } else if ( !strcmp(message.header, HDAC_REQ) ) {
            SEND_ARY_TO_STREAM(char, HDACCD);
        } else if ( !strcmp(message.header, HDAS_REQ) ) {
            SEND_ARY_TO_STREAM(short, HDASCD);
        } else if ( !strcmp(message.header, HDAUS_REQ) ) {
            SEND_ARY_TO_STREAM(ushort, HDAUSCD);
        } else if ( !strcmp(message.header, HDAUL_REQ) ) {
            SEND_ARY_TO_STREAM(ulong, HDAULCD);
        } else if ( !strcmp(message.header, HDAUI_REQ) ) {
            SEND_ARY_TO_STREAM(uint, HDAUICD);
        } else if ( !strcmp(message.header, HDALL_REQ) ) {
            SEND_ARY_TO_STREAM(llong, HDALLCD);
        } else if ( !strcmp(message.header, HDAULL_REQ) ) {
            SEND_ARY_TO_STREAM(ullong, HDAULLCD);
        } 
    }
    return 0;
}

/*
 * Functions for dealing with the header and trial parameters
 */

int dropHeader()
{
    ec_send_code_tagged(PARADIGMIDENTCD , VALOFFSET + PARADIGMID);
    gl_firsttrial = 0;
    return 0;
}

// currently just sets the dummy trial parameter to 1
int initTrialParams()
{
    int i;

    for ( i = 0; i < NTRIALPARAMS; ++i ) {
        gl_trialparams[i].flagcodeval = gl_trialparamcodes[i];
        gl_trialparams[i].casttype = gl_trialparamtypes[i];
        gl_trialparams[i].nelements = 1;
        switch (gl_trialparamtypes[i]) {
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
    }
    return 0;
}

int dropTrialParams()
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

/*
 * REX setup stuff
 */

// this function is called whenever you press Reset States
void rinitf()
{
    PARADIGM_PRELUDE("DataComm");

    macAllOff();
    wd_src_pos(WIND0, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_src_check(WIND0, WD_SIGNAL, EYEH_SIG, WD_SIGNAL, EYEV_SIG);
    wd_cntrl(WIND0, WD_ON);

    freeMemory();

    gl_firsttrial = 1;
    gl_mac_ready = 0;
    i_b->t_wrate = 1;

    updateEyewin(WIND0);
}

void updateEWHelper(int window, int loc_x, int loc_y, int size_x, int size_y)
{
    wd_src_pos(window, WD_DIRPOS, 0, WD_DIRPOS, 0);
    wd_pos(window, loc_x, loc_y);
    wd_siz(window, size_x, size_y);
}

void updateEyewin(int window)
{
    if ( window == WIND0 )
        updateEWHelper(window, gl_fp_x, gl_fp_y,
                gl_eyewin_x, gl_eyewin_y);
}

// Incoming codes are queued at index lx (_l_oad) and are dequeued starting at index dx (_d_ump).
// So lx is strictly greater than dx when there's codes waiting to be dumped, and lx == dx
// when the queue is empty.
int lopri_clear()
{
    static int prev_load_idx = 0;
    if (prev_load_idx != i_b->pl_lolx && i_b->pl_lodx == i_b->pl_lolx) {
        prev_load_idx = i_b->pl_lolx;
        return 1;
    }
    return 0;
}

/*
 * Functions specific to this paradigm
 */

int showFP()
{
    char buf[256];
    sprintf(buf, "ShowFP(%d, %d, %d, %d, %d, %d)",
            gl_fp_x, gl_fp_y, gl_fp_size,
            gl_fp_rgb[0], gl_fp_rgb[1], gl_fp_rgb[2]);
    udp_send_to(buf, "MAC");
    updateEyewin(WIND0);
    return 0;
}

int showTarg()
{
    char buf[256];
    sprintf(buf, "ShowTarg(%d, %d, %d, %d, %d, %d)",
                gl_rf_x, gl_rf_y, gl_targ_size,
                gl_targ_rgb[0], gl_targ_rgb[1], gl_targ_rgb[2]);
    udp_send_to(buf, "MAC");
    updateEyewin(WIND1);
    return 0;
}

int errorTrial()
{
    macAllOff();
    return 0;
}

int correctTrial()
{
    // here's where you increment the trial counter and
    // check if you've done enough trials for that block
    // or if you're done with the experiment

    if ( ++gl_trial_count >= 20 )
        gl_done = 1;
    return 0;
}

int freeMemory()
{
    // here's where dynamically allocated memory is freed
    return 0;
}

/*
 * Functions called by menus
 */

// nothing here at the moment

/*
 * Menu declarations
 */

VLIST menu_fp[] =
{
    {"X",             &gl_fp_x,      NP, NP, 0, ME_DEC},
    {"Y",             &gl_fp_y,      NP, NP, 0, ME_DEC},
    {"Size",          &gl_fp_size,   NP, NP, 0, ME_DEC},
    {"Red",           &gl_fp_rgb[0], NP, NP, 0, ME_DEC},
    {"Green",         &gl_fp_rgb[1], NP, NP, 0, ME_DEC},
    {"Blue",          &gl_fp_rgb[2], NP, NP, 0, ME_DEC},
    {"Window width",  &gl_eyewin_x,  NP, NP, 0, ME_DEC},
    {"Window height", &gl_eyewin_y,  NP, NP, 0, ME_DEC},
    {NS},
};

VLIST menu_targ[] =
{
    {"X",             &gl_rf_x,        NP, NP, 0, ME_DEC},
    {"Y",             &gl_rf_y,        NP, NP, 0, ME_DEC},
    {NS},
};

MENU umenus[] = {
    {"Fixation",    &menu_fp,    NP, NP, 0, NP, NP},
    {"Target",      &menu_targ,  NP, NP, 0, NP, NP},
    {NS},
};

VLIST state_vl[] = {{NS}};

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
    time 100
    to startSlave
startSlave:
    do startSlavePrgm()
    time 100
    to initMac
initMac:
    do setupMac()
    time 5000
    to getRFxy on 1 = gl_mac_ready
getRFxy:
    do getRFxy()
    time 10
    to initHeaderVars
initHeaderVars:
    do initHeaderVars()
/*    to initTrialVars on 1 = gl_mac_ready
initTrialVars:
    do initTrialVars()*/
    to pause0 on 1 = gl_mac_ready
pause0:
    to finished on 1 = gl_done
    to pause1 on +PSTOP & softswitch
    to initTrial
pause1:
    to initTrial on -PSTOP & softswitch
initTrial:
    do initTrialParams()
    time 200
    to finished on 1 = gl_done
    to plxStart
plxStart:
    do dio_on(PLXREC)
    time 100
    to dropHeader on 1 = gl_firsttrial
    to fpOn
dropHeader:
    do dropHeader()
    time 10
    to getHDi /* start the header param requests */
getHDi:
    do getHDi()
    to getHDd on 1 % lopri_clear
getHDd:
    do getHDd()
    to getHDf on 1 % lopri_clear
getHDf:
    do getHDf()
    to getHDl on 1 % lopri_clear
getHDl:
    do getHDl()
    to getHDc on 1 % lopri_clear
getHDc:
    do getHDc()
    to getHDs on 1 % lopri_clear
getHDs:
    do getHDs()
    to getHDus on 1 % lopri_clear
getHDus:
    do getHDus()
    to getHDul on 1 % lopri_clear
getHDul:
    do getHDul()
    to getHDui on 1 % lopri_clear
getHDui:
    do getHDui()
    to getHDll on 1 % lopri_clear
getHDll:
    do getHDll()
    to getHDull on 1 % lopri_clear
getHDull:
    do getHDull()
    to getHDai on 1 % lopri_clear
getHDai:
    do getHDai()
    to getHDad on 1 % lopri_clear
getHDad:
    do getHDad()
    to getHDaf on 1 % lopri_clear
getHDaf:
    do getHDaf()
    to getHDal on 1 % lopri_clear
getHDal:
    do getHDal()
    to getHDac on 1 % lopri_clear
getHDac:
    do getHDac()
    to getHDas on 1 % lopri_clear
getHDas:
    do getHDas()
    to getHDaus on 1 % lopri_clear
getHDaus:
    do getHDaus()
    to getHDaul on 1 % lopri_clear
getHDaul:
    do getHDaul()
    to getHDaui on 1 % lopri_clear
getHDaui:
    do getHDaui()
    to getHDall on 1 % lopri_clear
getHDall:
    do getHDall()
    to getHDaull on 1 % lopri_clear
getHDaull:
    do getHDaull()
    to fpOn on 1 % lopri_clear
fpOn:
    code FPONCD
    do showFP()
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
    do hideFP()
    to rewOn
rewOn:
    code REWCD
    do dio_on(REW)
    time 70
    to rewOff
rewOff:
    do dio_off(REW)
    to correctTrial
error:
    code ABORTCD
    do errorTrial()
    time 200
    to dropTrialParams
correctTrial:
    do correctTrial()
    to dropTrialParams
dropTrialParams:
    do dropTrialParams()
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
