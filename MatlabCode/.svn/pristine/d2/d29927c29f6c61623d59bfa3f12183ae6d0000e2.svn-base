#ifndef DATACOMM_INCLUDED
#define DATACOMM_INCLUDED

#define INT 0
#define LONG 1
#define DOUBLE 2
#define FLOAT 3
#define PARADIGMID 888 // change this to your paradigm id

// Codes that get dropped in the header (on the first trial only; 8000s range)
#define HDICD    8000
#define HDDCD    8001
#define HDFCD    8002
#define HDLCD    8003
#define HDCCD    8004
#define HDSCD    8005
#define HDUSCD   8006
#define HDULCD   8007
#define HDUICD   8008
#define HDLLCD   8009
#define HDULLCD  8010
#define HDAICD   8011
#define HDADCD   8012
#define HDAFCD   8013
#define HDALCD   8014
#define HDACCD   8015
#define HDASCD   8016
#define HDAUSCD  8017
#define HDAULCD  8018
#define HDAUICD  8019
#define HDALLCD  8020
#define HDAULLCD 8021
#define PARADIGMIDENTCD 8999

// Codes that get dropped on every trial (in the 7000s range)
#define TPICD    7000
#define TPDCD    7001
#define TPFCD    7002
#define TPLCD    7003
#define TPCCD    7004
#define TPSCD    7005
#define TPUSCD   7006
#define TPULCD   7007
#define TPUICD   7008
#define TPLLCD   7009
#define TPULLCD  7010
#define TPAICD   7011
#define TPADCD   7012
#define TPAFCD   7013
#define TPALCD   7014
#define TPACCD   7015
#define TPASCD   7016
#define TPAUSCD  7017
#define TPAULCD  7018
#define TPAUICD  7019
#define TPALLCD  7020
#define TPAULLCD 7021

typedef struct {
    int flagcodeval;
    int casttype;
    int nelements;
    union {
        long lval;
        float fval;
        double dval;
        long* lpval;
    } val;
} ecode_struct;

#define PARADIGM_PRELUDE( name ) \
do { \
    static int _____do_once = 1; \
    if ( _____do_once ) { \
        dprintf(name " compiled on %s at %s\n", __DATE__, __TIME__); \
        _____do_once = 0; \
    } \
} \
while ( 0 )

// the variable names in the following macros must match the declarations in checkMsg()
#define SEND_VAL_TO_STREAM( type, code ) \
do { \
    type##tmp[0] = hex2##type(message.contents, NULL); \
    sendToPlexDataStream(#type, (void*) &type##tmp[0], 1, code, 0); \
} while (0)

#define SEND_ARY_TO_STREAM( type, code ) \
do { \
    hex2##type(message.contents, (void*) &type##tmp); \
    sendToPlexDataStream(#type, (void*) &type##tmp, message.size, code, 1); \
} while (0)

#endif /* DATACOMM_INCLUDED */
