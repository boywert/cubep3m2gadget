#include <stdio.h>

int highword(long longword)
{
  return (longword >> 32);
}
int lowword(long longword)
{
  return  (longword & (0xFFFF));
}
