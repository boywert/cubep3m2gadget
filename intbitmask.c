#include <stdio.h>
void highword_(long *longword, int *high);
void lowword_(long *longword, int *low);
void highword_(long *longword, int *high)
{
  (*high) =((*longword) >> 32);
}
void lowword_(long *longword, int *low)
{
  (*low) = ((*longword) & (0xFFFF));
}
