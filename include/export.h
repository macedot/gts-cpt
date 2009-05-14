#ifndef __EXPORT_H__
#define __EXPORT_H__

size_t cpt_eulerian_mesh_write(SignedDistance *pSignedDistance, FILE *fp);
size_t cpt_eulerian_mesh_write_silo(SignedDistance *pSignedDistance, const char* szFileOut);
size_t cpt_eulerian_mesh_write_ze(SignedDistance *pSignedDistance, GtsSurface *surface, const char* szFileOut);

#endif // __CPT_H__
