enum xyzslicetype { xyslice = 0, xzslice = 1, yzslice = 2, volume = 3};
const static char *outdenname[4]={"xyslice", "xzslice", "yzslice", "volume"};

void volumeout(SimParameters *, GridInfo *, enum xyzslicetype);
