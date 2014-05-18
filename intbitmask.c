#include <stdio.h>

void highword(long *longword, int *high)
{
  (*high) =((*longword) >> 32);
}
void lowword(long *longword, int *low)
{
  (*low) = ((*longword) & (0xFFFF));
}
