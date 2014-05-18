#include <stdio.h>

int highword(long *longword)
{
  long w = *longword;
  return =(w >> 32);
}
int lowword(long *longword)
{
  long w = *longword;
  return (w & (0xFFFF));
}
