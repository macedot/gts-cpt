#ifndef __CPT_H__
#define __CPT_H__

typedef gdouble (*CptFuncDistance)(GtsPoint* mesh_point, gpointer data);

typedef struct param_distance_function {
	gpointer    pNormal;
	gpointer    pObj;
	GtsPoint*   pSurfacePoint;
} ParamDistFunc;


#define BOX_BEGIN(box_begin, point, normal, sign) \
	(box_begin)[0] = cpt_min((box_begin)[0], cpt_min((point)->x, (point)->x + (sign) * (normal)[0])); \
	(box_begin)[1] = cpt_min((box_begin)[1], cpt_min((point)->y, (point)->y + (sign) * (normal)[1])); \
	(box_begin)[2] = cpt_min((box_begin)[2], cpt_min((point)->z, (point)->z + (sign) * (normal)[2])); \

#define BOX_END(box_begin, point, normal, sign) \
	(box_end)[0]   = cpt_max((box_end)[0]  , cpt_max((point)->x, (point)->x + (sign) * (normal)[0])); \
	(box_end)[1]   = cpt_max((box_end)[1]  , cpt_max((point)->y, (point)->y + (sign) * (normal)[1])); \
	(box_end)[2]   = cpt_max((box_end)[2]  , cpt_max((point)->z, (point)->z + (sign) * (normal)[2]));

gdouble cpt_min(gdouble a, gdouble b);
gdouble cpt_max(gdouble a, gdouble b);

gdouble cpt_get_dist_cut(gdouble size_sup, gdouble delta_max);
gdouble cpt_get_dist_max(gdouble dist_cut, gdouble dist_extra);

gdouble cpt_vector_angle (GtsVector vector[]);

gdouble cpt_round_value(gdouble valueParam, gdouble valueMin);

int cpt_init_distance_function(SignedDistance* pSignedDistance);

void cpt_vector_normalize(GtsVector *pVector, gdouble factor);

void cpt_foreach_face(GtsTriangle *pTriangle, gpointer pEulerian);
void cpt_foreach_edge(GtsEdge *pEdge, gpointer pEulerian);
void cpt_foreach_vertex(GtsVertex *pVertex, gpointer pEulerian);

void surface_findbox_by_list(GtsSurface *surface, gpointer pEulerian);
void surface_findbox(GtsVertex *pVertex, gpointer pEulerian);

size_t cpt_eulerian_mesh_write(SignedDistance *pSignedDistance, FILE *fp);
size_t cpt_eulerian_mesh_write_silo(SignedDistance *pSignedDistance, const char* szFileOut);
size_t cpt_eulerian_mesh_write_ze(SignedDistance *pSignedDistance, GtsSurface *surface, const char* szFileOut);

#endif // __CPT_H__
