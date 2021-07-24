#include <stdio.h>
#include <stdint.h>

int main()
{
  uint32_t m1=((uint32_t)(1<<17))+((uint32_t) (1<<31));
  uint32_t m2=((uint32_t) (1<<(58-32)));

  printf("%u %u \n",~m1,~m2);
  
}
