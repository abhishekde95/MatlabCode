/*
 *
 * FUNCTIONS FOR DEALING WITH UDP COMMUNICATIONS ON REX.
 *
 *	CAH 12/07
 *  GDLH 1/08 created this header file
 *  CAH & GDLH Error correction and changed the name of the header file from CharlieUDP.h
 *  ZALB 06/14 Adding support for more data types
 */


/***************************/
/*       DEFINITIONS       */
/***************************/

#define VALOFFSET	4000
#define REXBUFF 8500
#define MAXHEADERSIZE 500
#define UINTSPERDOUBLE 8
#define UINTSPERFLOAT 4
#define UINTSPERLONG 4
typedef struct {
	char *contents;
	char header[MAXHEADERSIZE];
	long size;
	}message_struct;

#define my_ceil(A)  ((A)==(float)(int)(A) ? (int)(A) : (int)((A)+1))
#define my_round(A)  ((A)>=0 ? (int)((A)+0.5) : (int)((A)-0.5))

/***************************/
/*       DECLARATIONS      */
/***************************/

#define DECL_HEX2_FUNC( funcname, type ) \
type funcname(char* hex, void* dump) { \
    char *start, *end; \
    type* p = dump; \
    union { unsigned long long x; type y; } val ; \
    for (start = hex; *start != '\0'; start = end) { \
        val.x = strtoull(start, &end, 16); \
        if (p) *p++ = val.y; \
    } \
    return val.y; \
}

// for all of these hex2* functions, make sure your buffer has enough allocated space!
DECL_HEX2_FUNC(hex2ullong, unsigned long long);
DECL_HEX2_FUNC(hex2llong, long long);
DECL_HEX2_FUNC(hex2double, double);
DECL_HEX2_FUNC(hex2ulong, unsigned long);
DECL_HEX2_FUNC(hex2long, long);
DECL_HEX2_FUNC(hex2uint, unsigned int);
DECL_HEX2_FUNC(hex2int, int);
DECL_HEX2_FUNC(hex2float, float);
DECL_HEX2_FUNC(hex2ushort, unsigned short);
DECL_HEX2_FUNC(hex2short, short);
#undef DECL_HEX2_FUNC

char hex2char(char* hex, void* dump)
{
    char *start, *end, *p = dump;
    char c;
    for (start = hex; *start != '\0'; start = end) {
        c = strtol(start, &end, 16);
        if (p) *p++ = c;
    }
    if (p) *p = 0;
    return c;
}

int sendToPlexDataStream(const char * , const void * , long, int, int);  // Drop multiple codes to Plexon
void prepareMsg(message_struct *, char *);		// Split incoming message up into three fields
int udpBulkTransfer(message_struct, char *, char *);


/***************************/
/*       FUNCTIONS      */
/***************************/


/*
 *
 *  PREPARE THE MESSAGE WITH prepareMsg
 *
 */
void prepareMsg(message_struct * message, char * udpMessage){

	int a=0;
	int index[2];
	int placeHolder=0;
	int chunkSize = 0;
	int newSize = 0;
	char sizeInfo[50];   //should be way more than enough space!
	char inMsg[REXBUFF];
	char *tmpPointer;


	/* message should be a static message_struct. previous calls to this function could have
	   called malloc() so free up that space now!! While your at it, initialize all the
	   appropriate variable, otherwise.... bad juju. */
	free(message->contents);
	message->contents = NULL;   // to avoid "double free" errors
	message->header[0] = '\0';
	message->size = 0;
	index[1] = 0;


	/* find the indicies to the size information and the message contents */
	while(!index[1] && (a<strlen(udpMessage))){
		if (!strncmp(&udpMessage[a], ">>", 2)){
			index[placeHolder++] = a+2;              /* index to the first character after a ">>" */
		}
		a++;
	}

	/* fill up the message_struct */
	if (index[0] > (MAXHEADERSIZE-4)){
		printf("Header length exceeds MAXHEADERSIZE\n");
	} else {
		strncat(message->header, udpMessage, index[0]-2);
		strncat(sizeInfo, &udpMessage[index[0]], (index[1] - index[0] - 2)); /* -2 due to >> */
		message->size = strtol(sizeInfo, NULL, 16);
	}

	/* Did all the data get sent accross in one transmission? If so, stop here. Otherwise, reallocate space and proceed */
	chunkSize = strlen(&udpMessage[index[1]]);
	if (chunkSize){
		newSize = chunkSize;
	} else{
		newSize = 17 * message->size;
	}

	tmpPointer = realloc(message->contents, (size_t)(newSize+1));
	if (tmpPointer == NULL){
		printf("problem with realloc\n");
	} else{
		message->contents = (char *)tmpPointer;
		message->contents[0] = '\0';
		strncat(message->contents, &udpMessage[index[1]], chunkSize);
	}


	if (!chunkSize && message->size){
		while (1){
			inMsg[0] = '\0';
			udp_read(inMsg, REXBUFF);

			/* if the transmission is complete, break, otherwise concatinate udp messages into a character array */
			if (!strcmp(inMsg, "transferComplete")){
				break;
			}
			strncat(message->contents, inMsg, newSize-strlen(message->contents));
		}
	} else {
		return;
	}

	return; //after concatenating the entire message during the while loop
}

/*
 * This function extracts each byte from `data`, adds VALOFFSET, and sends the
 * result to the low priority buffer. These bytes end up in the PLEXON data stream
 * and our PLX/NEX files. Then nex2stro reconstructs the original array of values.
 */
int sendToPlexDataStream(const char* type, const void* data, long nelements, int eventcode, int bookend)
{
    int i, encode_directly = 0, width = 0;
    short msg = 0;

    int* intp = (int*) data;
    unsigned char* bytes = (unsigned char*) data;

    if (!strcmp(type, "double")
        || !strcmp(type, "ullong")
        || !strcmp(type, "llong")) {
        width = sizeof(double);
    } else if (!strcmp(type, "float")
        || !strcmp(type, "ulong")
        || !strcmp(type, "long")
        || !strcmp(type, "uint")) {
        width = sizeof(float); /* long is 4 bytes wide on the REX machines */
    } else if (!strcmp(type, "ushort")
        || !strcmp(type, "short")) {
        width = sizeof(short);
    } else if (!strcmp(type, "int")) {
        width = 1;
        encode_directly = 1;
    } else if (!strcmp(type, "char")) {
        width = sizeof(char);
    } else {
        printf("sendToPlexDataStream: unrecognized data type \"%s\"\n", type);
        return 0;
    }

    /* send the bytes and mark the beginning with the eventcode */
    ec_send_lo(eventcode, 1);
    for (i = 0; i < nelements * width; ++i) {
        if (encode_directly)
            msg = intp[i] + VALOFFSET;
        else
            msg = bytes[i] + VALOFFSET;

        ec_send_lo(msg, 1);

        if (msg < 2000 || msg > 6999) {
            printf("Event codes must be between 2000 and 6999."
                " Code #%d is equal to: %d (overall type is %s)\n",
                i, msg, type);
        }
    }

    /* mark the end with the eventcode to delimit the array of data in between */
    if (bookend) ec_send_lo(eventcode, 1);

    return 0;
}

/*
 *
 *  SEND A BULK UDP CHARACTER ARRAY
 *
 */
int udpBulkTransfer(message_struct message, char * hostName, char * functionCall){
	char outMsg[REXBUFF];
	int startPos = 0;
	int numChunks = 0;
	int a;


	/*initialize the strings */
	outMsg[0] = '\0';

	/* alert the appropriate slave program that a big message is about to come accross */
	udp_send_to(functionCall, hostName);

	/*determine the number of transmissions necessary to send the message in it's entirety */
	numChunks = my_ceil((double)strlen(message.contents) / (REXBUFF-1));


	/* send away! */
	for (a=1;a <= numChunks; a++){
		outMsg[0] = '\0';
		strncat(outMsg, &message.contents[startPos], (REXBUFF-1));
		udp_send_to(outMsg, hostName);
		startPos += (REXBUFF-1);
	}

	return(0);
}
