#ifndef CHARMM_PB_REEN_CUDA_H
#define CHARMM_PB_REEN_CUDA_H
#ifdef __cplusplus
extern "C" {
#endif
/* Accelerate only MAYER_REEN's accessible-probe loop and negative-grid restore.
 * CPU retains original GENPOINT/RNG, neighbor construction/QCKSRT and cleanup.
 * selected and neighbors contain the original one-based atom identifiers.
 * neighbors is Fortran(neighbor_stride,nselected), column major.
 * x/y/z/pbrad have natom entries. pox/poy/poz have mapt entries.
 * ml/nl have nselected entries; buried atoms have ml=0 and nl=-1.
 * qfkap must be 0 or 1. No final negative-grid cleanup is done here.
 * On return 0 all four output arrays are committed; on any error none changes.
 */
int charmm_pb_reen_cuda(
    int nx, int ny, int nz, int natom, int nselected, int mapt,
    int neighbor_stride, int qfkap,
    double dcel, double radw, double epsp, double rsmall,
    double tranx, double trany, double tranz,
    double xbcen, double ybcen, double zbcen,
    const int *selected, const int *ml, const int *nl, const int *neighbors,
    const double *x, const double *y, const double *z, const double *pbrad,
    const double *pox, const double *poy, const double *poz,
    float *rmay, float *rmayx, float *rmayy, float *rmayz,
    char *error, int error_len);
#ifdef __cplusplus
}
#endif
#endif
