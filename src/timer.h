#ifndef TIMER_SEEN
#define TIMER_SEEN

#include <time.h>


static clock_t timer_store[100];
static clock_t old_clock,zero_clock;

void init_timer();

void clear_timer(int no);

void add_timer(int i);

void print_timer(int i);

/* prints all the timers to stdout */

void print_timer_all(int n);


#endif
