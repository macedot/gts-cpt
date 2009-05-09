// version 1.0 : Jose Eduardo

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gts.h"
#include "gtstools.h"

#ifdef __cplusplus
}
#endif

#ifdef PI
#define PI 3.14159265359
#endif

/*
 *  Description: bluid a GSList.
 *  Status:      fully working.
 */
#undef __FUNC__
#define __FUNC__ "build_list"
void build_list(gpointer data, GSList ** list)
{
	*list = g_slist_prepend (*list, data);
}

/*
 *  Description: calculate curvature at each vertex of a closed surface.
 *  Status:      working, print the curvature
 */
#undef __FUNC__
#define __FUNC__ "Curvature_Test"
void Curvature_Test(GtsSurface *surface)
{
	GSList		*vertex = NULL;
	GSList		*item	= NULL;
	GtsVertex	*v		= NULL;
	GtsVector	Kh;
	gboolean	success;
	double		khm;
	//int			count	= 0;


	// build a list of vertes
	gts_surface_foreach_vertex(surface, (GtsFunc) build_list, &vertex);

	// walk by the list
	for(item = vertex ; item ; item = item->next)
	{
		v = GTS_VERTEX(item->data);

		// Mean curvature in each triangle vertex
		success = gts_vertex_mean_curvature_normal(v, surface, Kh);
		khm = 0.5*gts_vector_norm(Kh);

		printf ("Curv: %f \n", khm);

		//count++;
	}
	g_slist_free (vertex);
}

/*
 *  Description: scale a surface, (GtsFunc)ScaleSurface // gts_surface_foreach_vertex
 *  Status:      fully working.
 */
#undef __FUNC__
#define __FUNC__ "ScaleSurface"
void ScaleSurface(GtsVertex * v, double data[4])
{
	// coord times the radius;
	v->p.x *= data[3];
	v->p.y *= data[3];
	v->p.z *= data[3];
}

/*
 *  Description: move a surface, (GtsFunc)MoveSurface // gts_surface_foreach_vertex
 *  Status:      fully working.
 */
#undef __FUNC__
#define __FUNC__ "MoveSurface"
void MoveSurface(GtsVertex * v, double data[4])
{
	v->p.x += data[0];
	v->p.y += data[1];
	v->p.z += data[2];
}


