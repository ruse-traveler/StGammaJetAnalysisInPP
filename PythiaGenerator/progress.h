// 'progress.h'
// Derek Anderson
// 06.03.2015
//
// This header file contains a couple of functions for outputting a progress indicator
// and some time stamps. 'cerr' guarantees that output is streamed to the command line.
//
// Last updated: 03.08.2016

#pragma once

#include <ctime>

using namespace std;



// progress bar ---------------------------------------------------------------

static inline void progress(int m, int n) {

   int w = 50;                // total width of the bar

   float p = m / (float) n;   // ratio of complete to incomplete
   int   l = p * w;           // length of bar to display

   cerr << setw(3) << (int) (p * 100) << "% [";
   for (int i = 1; i < l; ++i)   cerr << "=";
   for (int i = l+1; i < w; ++i) cerr << " ";
   cerr << "] \r" << flush;

   if (m == n - 1) cerr << endl;

}


// start / end time -----------------------------------------------------------

void timeBegin() {

   time_t start;
   struct tm *local;

   time(&start);
   local = localtime(&start);
   cerr << "Start: " << asctime(local);

}


void timeEnd() {

   time_t end;
   struct tm *local;

   time(&end);
   local = localtime(&end);
   cerr << "End:   " << asctime(local);

}

// End ------------------------------------------------------------------------
