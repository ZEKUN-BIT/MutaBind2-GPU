#ifndef CHARMMPBEQ_CUDA_H
#define CHARMMPBEQ_CUDA_H
#ifdef __cplusplus
extern "C" {
#endif
/* Bounds are the adjusted Fortran one-based nonperiodic PBEQ1 bounds.
 * One iteration means one checkerboard color, as in PBEQ1.
 * Return 0 after a successful GPU calculation, even if MAXITS was reached.
 * iterations is MAXITS+1 on nonconvergence, otherwise the converged iteration.
 * On any error return nonzero, describe it in error, and leave host PHI intact.
 * Caller checks the returned convergence information before accepting a solve.
 */
int charmm_pbeq1_cuda(
    int nx, int ny, int nz,
    int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax,
    int maxits, double omega, double tol, double epsw, double rjac,
    float *phi, const float *epsx, const float *epsy, const float *epsz,
    const float *cden, const float *fkapa,
    int *iterations, double *anorm, char *error, int error_len);
/* Same solve, additionally record the original per-color PRNLEV>10 trace.
 * Arrays must all have trace_capacity>=MAXITS. Omega is recorded before its
 * recurrence update, and center_phi is recorded after the color update.
 */
int charmm_pbeq1_cuda_trace(
    int nx, int ny, int nz,
    int ixmin, int ixmax, int iymin, int iymax, int izmin, int izmax,
    int maxits, double omega, double tol, double epsw, double rjac,
    float *phi, const float *epsx, const float *epsy, const float *epsz,
    const float *cden, const float *fkapa,
    int *iterations, double *anorm, char *error, int error_len,
    double *norm_trace, double *omega_trace, float *center_trace,
    int trace_capacity);
#ifdef __cplusplus
}
#endif
#endif
