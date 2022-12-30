//
// Pibero Djawotho <pibero@iucf.indiana.edu>
// Indiana University
// March 15, 2006
//
// L2Result must not exceed 6 words
//

#ifndef L2_GAMMARESULT_H
#define L2_GAMMARESULT_H

struct L2gammaResult {
   unsigned char result_version;           /* 1-15 */

  unsigned char threshold;                /* 0x1=ht1 0x2=tp1 0x4=ht2 0x8=tp2 */
  unsigned char elapsed;                  /* cpu time [kTicks], 255=overflow for this event */

  unsigned char trigger;                  /* 0x0=abort, 0x1=ht 0x2=tp 0x4=preAccept 0x8=accept */
  unsigned char phibin;                   /* 0=none, 0-119 bemc, 0-59 eemc, 0x8=eemc bit  */
  unsigned char etabin;                   /* 0=none, 1-40 bemc, 1-12 eemc             */
  unsigned char pttowerx2;                /* E_T of the high tower times 2              */
  unsigned char ptclusterx2;              /* E_T of the cluster times 2               */

};
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------

inline void print_L2gammaResult( const L2gammaResult &result )
{
  printf("----------------------------------------------------------------\n");
  printf("L2gammaResult version %d\n",          /* print header version */
	 result.result_version);
  printf("bemc trig=%d / eemc trig=%d\n\n",     /* is trigger b/eemc based? */
	 !(result.phibin&0x8),
	 result.phibin&0x8);
  printf("ht1 tested=%d / ht2 tested=%d\n",     /* what thresholds tested? */
	 result.threshold&0x1,
	 result.threshold&0x4);
  printf("cl1 tested=%d / cl2 tested=%d\n", 
	 result.threshold&0x2,
	 result.threshold&0x8);
  printf("ht met=%d cl met=%d trig met=%d\n",   /* what threshold met? */
	 result.trigger&0x1,
	 result.trigger&0x2,
	 result.trigger&0x8);
  printf("phibin=%d etabin=%d\n",               /* and where? */
	 result.phibin&0xef,
	 result.etabin);
  printf("cpu time [kTicks]=%d\n",              /* and how long? */
	 result.elapsed);

};


#endif
