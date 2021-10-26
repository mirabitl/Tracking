#ifndef TDC_MAPPING_HH
#define TDC_MAPPING_HH
#include<map>
#include<stdint.h>

//const int LEMO2PR[32]={14,13,12,11,10,9,8,7,6,4,2,0,1,3,5,24,
//		 26,28,30,31,29,27,25,23,22,21,20,19,18,17,16,15};
//const int PR2LEMO[32]={11,12,10,13,9,14,8,7,6,5,4,3,2,1,0,31,
//		 30,29,28,27,26,25,24,23,15,22,16,21,17,20,18,19};
//const int PR2LEMO[32]={11,12,10,13,9,14,8,7,6,5,4,3,2,1,0,31,
//		 30,29,28,27,26,25,24,23,15,22,16,21,17,20,18,19};
#define FW_RETURN
#ifdef FW_COAX
// From injection tests
const int TDC2LEMO[32]={19,11,18,12,20,10,17,13,21,9,16,14,22,8,15,7,23,6,24,5,25,3,26,4,-1,-1,-1,-1,-1,-1,-1,-1};

const int LEMO2TDC[32]={-1,-1,-1,21,23,19,17,15,13,9,5,1,3,7,11,
			14,10,6,2,0,4,8,12,16,18,20,22,-1,-1,-1,-1,-1};
#define FW_COAX_INITIAL
#ifdef FW_COAX_INITIAL
// From Xushian firmware

const int PR2TDC[32]={1,3,5,7,9,11,13,15,
		17,19,21,23,-1,-1,-1,-1,
		      -1,-1,-1,-1,22,20,18,16,
		14,12,10,8,6,4,2,0};

const int TDC2PR[32]={31,0,30,1,29,2,28,3,27,4,26,5,25,6,24,7,
		 23,8,22,9,21,10,20,11,-1,-1,-1,-1,-1,-1,-1,-1,};

const int PR2LEMO[32]={11,12,10,13,9,14,8,7,6,5,3,4,-1,-1,-1,-1,-1,-1,-1,-1,26,25,24,23,15,22,16,21,17,20,18,19
};

const int LEMO2PR[32]={-1,-1,-1,10,11,9,8,7,6,4,2,0,1,3,5,24,26,28,30,31,29,27,25,23,22,21,20,-1,-1,-1,-1,-1};
#endif
#ifdef FW_COAX_CORRECTED
// From Sirley firmware
const int PR2TDC[32]={1,3,5,7,9,11,13,15,
		17,19,21,-1,-1,-1,-1,-1,
		      -1,-1,-1,23,22,20,18,16,
		14,12,10,8,6,4,2,0};

const int TDC2PR[32]={31,0,30,1,29,2,28,3,27,4,26,5,25,6,24,7,
		      23,8,22,9,21,10,20,19,-1,-1,-1,-1,-1,-1,-1,-1,};

const int PR2LEMO[32]={11,12,10,13,9,14,8,7,6,5,3,4,-1,-1,-1,-1,-1,-1,-1,-1,26,25,24,23,15,22,16,21,17,20,18,19
};

const int LEMO2PR[32]={-1,-1,-1,10,11,9,8,7,6,4,2,0,1,3,5,24,26,28,30,31,29,27,25,23,22,21,20,-1,-1,-1,-1,-1};

#endif
// So we can build


// Example connection to be put in JSON geometry 
const int LEMO2STRIP[32]={86,85,84,83,82,81,80,79,78,77,76,75,74,73,72,71,
		    71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86};

const int FEB2STRIP[24]={255,255,255,255,255,0,24,36,12,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255};

// 1	10	14
// 2	6	11
// 3	2	7
// 4	0	3
// 5	4	1
// 6	8	5
// 7	12	9
// 8	16	13
// 9	18	15
// 10	20	17
// 11	22	19
// 12	missing (23)	21
const int SIDE[48]={1,0,1,0,
		    1,0,1,0,
		    1,0,1,0,
		    1,0,0,0,
		    1,0,1,0,
		    1,0,1,1,1,0,1,0,
		    1,0,1,0,
		    1,0,1,0,
		    1,0,0,0,
		    1,0,1,0,
		    1,0,1,1};
const int STRIP[24]={4,
		     5,
		     3,
		     4,
		     5,
		     6,
		     2,
		     3,6,7,1,2,7,8,1,9,8,10,9,11,10,12,11,12};


#endif
#ifdef FW_RETURN

// From injection tests
const int TDC2LEMO[32]={19,11,18,12,20,10,17,13,21,9,16,14,22,8,15,7,23,6,24,5,25,3,26,4,-1,-1,-1,-1,-1,-1,-1,-1};

const int LEMO2TDC[32]={-1,-1,-1,21,23,19,17,15,13,9,5,1,3,7,11,
			14,10,6,2,0,4,8,12,16,18,20,22,-1,-1,-1,-1,-1};

// From Sirley firmware


const int PR2TDC[32]={31,29,27,25,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,-1,-1,-1,-1,-1,-1};

const int TDC2PR[32]={32*-1};

// So we can build

const int PR2LEMO[32]={11,12,10,13,9,14,8,7,6,5,3,4,-1,-1,-1,-1,-1,-1,-1,-1,26,25,24,23,15,22,16,21,17,20,18,19
};

const int LEMO2PR[32]={-1,-1,-1,10,11,9,8,7,6,4,2,0,1,3,5,24,26,28,30,31,29,27,25,23,22,21,20,-1,-1,-1,-1,-1};

// Example connection to be put in JSON geometry 
const int LEMO2STRIP[32]={86,85,84,83,82,81,80,79,78,77,76,75,74,73,72,71,
		    71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86};

const int FEB2STRIP[24]={255,255,255,255,0,255,255,255,255,12,24,36,255,255,255,255,255,255,255,255,255,255,255,255};

// 0 B20 31 1 L
// 1 B19 29 1 H
// 2 C18 27 2 L
// 3 D16 25 2 H
// 4 D15 23 3 L
// 5 C14 22 3 H
// 6 E11 21 4 L
// 7 D11 20 4 H
// 8 E9 19 5 L
// 9 D9 18 5 H
// 10 G7 17 6 L
// 11 H7 16 6 H
// 12 B6 15 7 L
// 13 A6 14 7 H
// 14 B5 13 8 L
// 15 A5 12 8 H
// 16 B4 11 9 L
// 17 A4 10 9 H
// 18 B3 9 10 L
// 19 A3 8 10 H
// 20 B18 7 11 L
// 21 A18 6 11 H

const int SIDE[48]={1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,
		    1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0};
const int STRIP[24]={1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12};

#endif
// #define DEBUG_PRINT_ENABLED 1  // uncomment to enable DEBUG statements
#define INFO_PRINT_ENABLED 1
#if DEBUG_PRINT_ENABLED
#define INFO_PRINT_ENABLED 1
#define DEBUG_PRINTF printf
#else
#define DEBUG_PRINTF(format, args...) ((void)0)
#endif
#if INFO_PRINT_ENABLED
#define INFO_PRINTF printf
#else
#define INFO_PRINTF(format, args...) ((void)0)
#endif


#endif
