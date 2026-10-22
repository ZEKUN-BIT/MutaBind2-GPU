#ifndef CHARMM_PB_BOUNDARY_CUDA_H
#define CHARMM_PB_BOUNDARY_CUDA_H
#ifdef __cplusplus
extern "C" {
#endif
/* Replacement for only MAYER_BP_INT's atom loop, restricted to:
 * odd grid dimensions >=5, zero salt, no membrane/images/basis selection.
 * Sources are packed in original LSTPRP order, with original double charge
 * and already transformed xi=x(i)+tranx-xbcen (likewise yi,zi).
 * Zero charges may be retained. PHI is float, z-contiguous, and is not cleared.
 * CPU caller performs the original bilinear interpolation afterward.
 * 0 = success; nonzero = error with original host PHI unchanged.
 */
int charmm_pb_boundary_cuda(
    int nx, int ny, int nz, int ncharges, double dcel, double epsw,
    const double *x, const double *y, const double *z, const double *q,
    float *phi, char *error, int error_len);
#ifdef __cplusplus
}
#endif
#endif
