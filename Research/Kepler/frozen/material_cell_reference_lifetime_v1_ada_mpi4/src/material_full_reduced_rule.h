#ifndef MATERIAL_FULL_REDUCED_RULE_H
#define MATERIAL_FULL_REDUCED_RULE_H
#include <stdint.h>
#include <stddef.h>
extern "C" int mf_evolved_rule_upload(void*,const uint64_t*,size_t,const uint32_t*,const double*,size_t);
extern "C" int mf_evolved_rule_use(void*,int);
#endif
