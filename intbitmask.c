#include <stdio.h>

int highword_(long *longword)
{
  long w = *longword;
  printf("high:%ld\n",w);
  return (w >> 32);
}
int lowword_(long *longword)
{
  long w = *longword;
  printf("low:%ld\n",w);
  return (w & (0xFFFF));
}
