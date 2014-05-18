#include <stdio.h>

int highword_(long *longword)
{
  long w = *longword;
  return =(w >> 32);
}
int lowword_(long *longword)
{
  long w = *longword;
  return (w & (0xFFFF));
}
