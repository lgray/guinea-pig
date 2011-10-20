#include "memory.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

/**********************************************************************/
/* memory_account routines                                            */
/**********************************************************************/


/* Storage and routines for the memory requirements */


void init_memory_account(MEMORY_ACCOUNT *account,long size)
{
  account->size=size;
  account->calls=0;
  account->amount=0;
  account->first=NULL;
}

void* get_memory(MEMORY_ACCOUNT *account,int size)
{
  MEMORY_POINTER *temp;
  //  printf(" get_memory : size= %d\n", size);
  //  printf(" get_memory : account->size= %d\n", account->size);
  if (account->size>0)
    {
      printf(" get_memory : account->amount = %ld\n", account->amount);
      if (account->size < account->amount+size)
	{
	  fprintf(stderr,"error: accout exceeded\n");
	  fflush(stderr);
	  exit(10);
	}
    }
  account->calls++;
  account->amount += size;
  temp=(MEMORY_POINTER*)malloc(sizeof(MEMORY_POINTER));
  if (temp==NULL){
      fprintf(stderr,"not eneough memory\n");
      fprintf(stderr,"%d bytes required\n",size);
      fflush(stderr);
      exit(10);
  }
  temp->next=account->first;
  temp->memory=malloc(size);
  if(temp->memory==NULL){
      fprintf(stderr,"not eneough memory\n");
      fprintf(stderr,"%d bytes required\n",size);
      fflush(stderr);
      exit(10);
  }
  account->first=temp;
  return temp->memory;
}

void clear_memory_account(MEMORY_ACCOUNT *account)
{
  MEMORY_POINTER *temp;
  account->amount=0;
  account->calls=0;
  temp=account->first;
  while(temp!=NULL)
    {
      account->first=temp->next;
      free(temp->memory);
      free(temp);
      temp=account->first;
    }
}

long memory_account_amount(MEMORY_ACCOUNT *account)
{
  return account->amount;
}

long memory_account_size(MEMORY_ACCOUNT *account)
{
  return account->size;
}

long memory_account_calls(MEMORY_ACCOUNT *account)
{
  return account->calls;
}

