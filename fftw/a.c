#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>

int main(void)
{
int i;
for (i = 0; i < 100000; ++i) {
void *p = memalign(32, rand() % 100 + 1);
if (((size_t) p) % 32 != 0) {
fprintf(stderr, "misaligned alloc %p\n", p);
return 1;
}
}
printf("success.\n");
return 0;
}
