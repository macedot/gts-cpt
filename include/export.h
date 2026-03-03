#ifndef GTS_CPT_EXPORT_H
#define GTS_CPT_EXPORT_H

size_t cpt_eulerian_mesh_write(SignedDistance *pSignedDistance, FILE *fp);
size_t cpt_eulerian_mesh_write_silo(SignedDistance *pSignedDistance, const char* szFileOut);
size_t cpt_eulerian_mesh_write_ze(SignedDistance *pSignedDistance, GtsSurface *surface, const char* szFileOut);

#endif // GTS_CPT_EXPORT_H
