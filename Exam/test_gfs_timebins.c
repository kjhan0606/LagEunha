#include "gfs_timebins.h"
#include <stdio.h>

#define CHECK(c) do { if (!(c)) { \
    fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #c); \
    return 1; \
} } while (0)

int main(void)
{
    /* Symmetric six-cell chain; closure must propagate from the fine end. */
    const size_t offset[] = {0, 1, 3, 5, 7, 9, 10};
    const uint32_t neighbor[] = {1, 0, 2, 1, 3, 2, 4, 3, 5, 4};
    const double dt[] = {1.0, 0.25, 0.125, 0.5, 0.25, 0.0625};
    const unsigned expected[] = {1, 2, 3, 2, 3, 4};
    unsigned level[6] = {99, 99, 99, 99, 99, 99};

    CHECK(gfs_timebins_assign(6, dt, 1.0, 4, offset, neighbor, 10, level)
          == GFS_TIMEBINS_OK);
    for (size_t i = 0; i < 6; ++i) {
        CHECK(level[i] == expected[i]);
        CHECK(gfs_timebins_interval(1.0, level[i]) <= dt[i]);
        for (size_t q = offset[i]; q < offset[i + 1]; ++q) {
            unsigned a = level[i], b = level[neighbor[q]];
            CHECK(a <= b + 1u && b <= a + 1u);
        }
    }

    /* Exercise the integer event schedule across one macro interval. */
    unsigned active_updates = 0;
    for (uint64_t tick = 0; tick < 16; ++tick)
        for (size_t i = 0; i < 6; ++i)
            active_updates += (unsigned)gfs_timebins_active(level[i], 4, tick);
    CHECK(active_updates == 42);
    CHECK(gfs_timebins_active(0, 4, 16));
    CHECK(!gfs_timebins_active(0, 4, 17));

    /* Invalid input and insufficient depth must not partially publish bins. */
    unsigned sentinel[6] = {91, 92, 93, 94, 95, 96};
    const double bad_dt[] = {1.0, 0.0, 0.125, 0.5, 0.25, 0.0625};
    CHECK(gfs_timebins_assign(6, bad_dt, 1.0, 4, offset, neighbor, 10,
                              sentinel) == GFS_TIMEBINS_INVALID);
    for (size_t i = 0; i < 6; ++i) CHECK(sentinel[i] == 91u + i);
    const double too_small[] = {1.0, 0.25, 0.125, 0.5, 0.25, 0.001};
    CHECK(gfs_timebins_assign(6, too_small, 1.0, 4, offset, neighbor, 10,
                              sentinel) == GFS_TIMEBINS_TOO_FINE);
    for (size_t i = 0; i < 6; ++i) CHECK(sentinel[i] == 91u + i);
    const uint32_t bad_neighbor[] = {99, 0, 2, 1, 3, 2, 4, 3, 5, 4};
    CHECK(gfs_timebins_assign(6, dt, 1.0, 4, offset, bad_neighbor, 10,
                              sentinel) == GFS_TIMEBINS_INVALID);
    for (size_t i = 0; i < 6; ++i) CHECK(sentinel[i] == 91u + i);

    puts("PASS gfs_timebins: 2:1 closure, exact power-of-two schedule, transactional errors");
    return 0;
}
