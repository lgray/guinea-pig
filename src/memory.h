#ifndef MEMORY_SEEN
#define  MEMORY_SEEN

struct memory_pointer
{
  void *memory;
 struct memory_pointer  *next;
};
typedef struct memory_pointer MEMORY_POINTER;


typedef struct
{
  long amount,calls,size;
  MEMORY_POINTER *first;
}MEMORY_ACCOUNT;


/* extern int store_memory_monitor[10]; */

void init_memory_account(MEMORY_ACCOUNT*,long);
void* get_memory(MEMORY_ACCOUNT*,int);
void clear_memory_account(MEMORY_ACCOUNT*);
long memory_account_amount(MEMORY_ACCOUNT*);
long memory_account_size(MEMORY_ACCOUNT*);
long memory_account_calls(MEMORY_ACCOUNT*);


#endif
