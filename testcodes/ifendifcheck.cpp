#include <stdio.h>

#define C_LANG    'D'
#define NO_ERROR  0

int main(void)
{
  int x = 4;
  if (x == 4)
  {
    #define B_LANG    'B'
  }
   #if C_LANG == 'C' && B_LANG == 'B'
     #define C_LANG_VALUE "I know the C language.\n"
     #define B_LANG_VALUE "I know BASIC.\n"
     printf("%s%s", C_LANG_VALUE, B_LANG_VALUE);
   #elif C_LANG == 'C'
     #define C_LANG_VALUE "I only know C language.\n"
     printf("%s", C_LANG_VALUE);
   #elif B_LANG == 'B'
     #define B_LANG_VALUE "I only know BASIC.\n"
     printf("%s", B_LANG_VALUE);
   #else
     printf("I don't know C or BASIC.\n");
   #endif

   return NO_ERROR;
}