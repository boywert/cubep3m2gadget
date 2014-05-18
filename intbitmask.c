#include <stdio.h>

int highword_(long *longword)
{
  long w = *longword;
  int k;
  k = (w >> 32);
  printf("high:%ld->%d\n",w,k);
  return k;
}
int lowword_(long *longword)
{
  long w = *longword;
  int k;
  printf("low:%ld->%d\n",w,k);
  k =  (w & (0xFFFF));
  printf("low:%ld->%d\n",w,k);
  return k;
}
