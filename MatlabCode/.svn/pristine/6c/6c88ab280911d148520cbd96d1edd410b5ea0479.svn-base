/* WhiteNoise.h - WhiteNoise specific ecodes */

/* Codes that get dropped in the header (on the first trial only) */
#define FIXXCD	8000
#define FIXYCD		8001
#define FPSIZECD		8002
#define FPRCD		8003
#define FPGCD		8004
#define FPBCD		8005
#define EYEWINXCD	8006
#define EYEWINYCD	8007
#define RFXCD		8008
#define RFYCD		8009
#define BKGNDRGBCD		8010
#define GAMMATABLECD	8011
#define MONSPDCD	8012
#define FRAMERATECD		8013
#define FUNDAMENTALSCD		8014
#define INITPROTOCOLCD		8015
#define NTRLSPERCONDCD		8016
#define INITDIAMCD	8017
#define INITSFCD	8018
#define INITTFCD	8019
#define INITORIENTCD	8020
#define INITPHASECD	8021
#define INITLCONTCD	8022
#define INITMCONTCD	8023
#define INITSCONTCD	8024
#define INITNCYCLESCD	8025
#define NSIGMASCD			8026
#define NCYCLESPLATCD		8027
#define NCYCLESRAMPCD		8028
#define MAXCONTSCALEFACTCD		8029
#define MINCONTSCALEFACTCD		8030
#define NGABORREPSCD		8031
#define NCONTRASTSCD		8032
#define STIMDURCD		8033
#define OPTOFLAGCD		8034



#define PARADIGMIDENTCD	8999

/* Codes that get dropped on every trial */
#define DIAMCD				7000
#define SFCD					7001
#define TFCD					7002
#define LCONTCD				7003
#define MCONTCD			7004
#define SCONTCD				7005
#define ORIENTCD			7006
#define PHASECD				7007
#define NCYCLESCD			7008
#define NFRAMESCD			7009
#define PROTOCOLCD		7010
#define STIMTYPECD			7011

typedef struct {
	int flagcodeval;
	int casttype;
	int nelements;
	union  {
		long lval;
		double dval;
		long *lpval;
	} val;
} ecode_struct;

#define INT 0
#define LONG 1
#define DOUBLE 2
#define PARADIGMID 150