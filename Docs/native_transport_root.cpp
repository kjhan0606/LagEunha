// Ordered kinetic mass-root primitives. Compile without fast-math or contraction.
#include <cmath>
#include <cstdint>
#include <vector>

extern "C" int lag_transport_root(
    int64_t n, int64_t cells, const int64_t *i, const int64_t *j,
    const double *scale, const double *base, const double *distance,
    const double *vl, const double *vr, const double *rl, const double *rr,
    const double *z, double *wl, double *wr, double *wn, double *fm,
    double *primitive, double *rho, double *residual, double *weight) {
  if (n < 1 || cells < 1) return 1;
  // Independent node arithmetic may parallelize; owner reductions remain ordered.
  #pragma omp parallel for schedule(static)
  for (int64_t k = 0; k < n; ++k) {
    double speed = base[k];
    if (j[k] >= 0) speed += (z[i[k]] - z[j[k]]) / (2.0 * distance[k]);
    wn[k] = speed;
    double sum_l = 0.0, sum_r = 0.0, abs_l = 0.0, abs_r = 0.0;
    double active_l = 0.0, active_r = 0.0;
    for (int a = 0; a < 4; ++a) {
      const int64_t q = 4*k + a;
      const double sl = vl[q] - speed, sr = vr[q] - speed;
      const double left = (0.25 * rl[k]) * (sl > 0.0 ? sl : 0.0);
      const double right = (0.25 * rr[k]) * (sr < 0.0 ? sr : 0.0);
      wl[q] = left; wr[q] = right;
      sum_l += left; sum_r += right;
      abs_l += std::abs(left); abs_r += std::abs(right);
      active_l += sl > 0.0 ? 1.0 : (sl < 0.0 ? 0.0 : 0.5);
      active_r += sr < 0.0 ? 1.0 : (sr > 0.0 ? 0.0 : 0.5);
    }
    fm[k] = scale[k] * (sum_l + sum_r);
    primitive[k] = scale[k] * (abs_l + abs_r);
    rho[k] = 0.25 * (rl[k]*active_l + rr[k]*active_r);
    weight[k] = j[k] >= 0 ? scale[k]*rho[k]/(2.0*distance[k]) : 0.0;
  }
  // Match the TWO bincount accumulations and their final subtraction exactly.
  std::vector<double> left(cells, 0.0), right(cells, 0.0);
  for (int64_t k = 0; k < n; ++k) left[i[k]] += fm[k];
  for (int64_t k = 0; k < n; ++k) if (j[k] >= 0) right[j[k]] += fm[k];
  for (int64_t k = 0; k < cells; ++k) residual[k] = left[k] - right[k];
  return 0;
}
