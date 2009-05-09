// cpt.cpp

#include "gtscpt.h"
#include "cpt.h"
#include <silo.h>
#include <VisItControlInterface_V1.h>

//#define CPT_VERBOSE_DEBUG

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
typedef gdouble (*CptFuncDistance)(GtsPoint* mesh_point, gpointer data);

typedef struct param_distance_function {
	gpointer    pNormal;
	gpointer    pObj;
	GtsPoint*   pSurfacePoint;
} ParamDistFunc;

#ifndef SQRT_3
static const gdouble SQRT_3 = sqrt(3.0);
#endif

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/******************************************************************************
 *
 ******************************************************************************/
gdouble cpt_min(gdouble a, gdouble b)
{
	return (a < b ? a : b);
}
/******************************************************************************
 *
 ******************************************************************************/
gdouble cpt_max(gdouble a, gdouble b)
{
	return (a > b ? a : b);
}
/******************************************************************************
 * epson = sqrt(3) * (suporte da delta de dirac / 2) * h;
 * h = delta maximo
 * suporte da delta de dirac = 2
 ******************************************************************************/
gdouble cpt_get_dist_cut(gdouble size_sup, gdouble delta_max)
{
	return SQRT_3 * 0.5 * size_sup * delta_max;
}
/******************************************************************************
 *
 ******************************************************************************/
gdouble cpt_get_dist_max(gdouble dist_cut, gdouble dist_extra)
{
	static const gdouble c = 1.0;
	return dist_cut + c * dist_extra;
}

/**
 * cpt_vector_angle:
 * @vector: #GtsVector of vectors.
 *
 * Returns: the value (in radians) of the angle between @t1 and @t2.
 * NOTE: This function is similar to the function function presented as 
 *       gts_triangle_angle at GTS library, except this version do not check 
 *       if theta <  0.0;
 */
gdouble cpt_vector_angle (GtsVector v1, GtsVector v2)
{
	gdouble pvx, pvy, pvz;
	gdouble theta;

	pvx = v1[1]*v2[2] - v1[2]*v2[1];
	pvy = v1[2]*v2[0] - v1[0]*v2[2];
	pvz = v1[0]*v2[1] - v1[1]*v2[0];

	theta = atan2 (sqrt (pvx*pvx + pvy*pvy + pvz*pvz) ,
	               v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);

	return theta;
}

/**
 * cpt_point_angle:
 */
gdouble cpt_point_angle(GtsPoint *p1, GtsPoint *p2, GtsPoint *p3)
{
	gdouble  n1x, n1y, n1z, n2x, n2y, n2z;
	gdouble  dist;
	gdouble  angle;

	dist = gts_point_distance(p2, p1);
	n1x = (p2->x - p1->x) / dist;
	n1y = (p2->y - p1->y) / dist;
	n1z = (p2->z - p1->z) / dist;

	dist = gts_point_distance(p3, p1);
	n2x = (p3->x - p1->x) / dist;
	n2y = (p3->y - p1->y) / dist;
	n2z = (p3->z - p1->z) / dist;

	angle = acos(n1x*n2x + n1y*n2y + n1z*n2z);
	return angle;
}

/**
 * cpt_vertex_pseudonormal:
 * @v: #GtsVector
 */
void cpt_vertex_pseudonormal(GtsVertex *v, gdouble *x, gdouble *y, gdouble *z)
{
	GtsTriangle  *t;
	GSList       *item;
	GSList       *triangles;
	GtsVertex    *v1, *v2, *v3;
	GtsVector    normalPseudo = {0.0 , 0.0 , 0.0};
	GtsVector    normalFace;
	gdouble      angle;

	g_return_if_fail(v != NULL);

	triangles = gts_vertex_triangles(v, NULL);

	// Angle between edges shared vertex;
	for ( item = triangles ; item ; item = item->next )
	{
		t = GTS_TRIANGLE(item->data);
		gts_triangle_vertices(t, &v1, &v2, &v3);

		if ( v == v1 )
		{
			angle = cpt_point_angle (GTS_POINT(v), GTS_POINT(v2), GTS_POINT(v3));
		}
		else if ( v == v2 )
		{
			angle = cpt_point_angle (GTS_POINT(v), GTS_POINT(v3), GTS_POINT(v1));
		}
		else if ( v == v3 )
		{
			angle = cpt_point_angle (GTS_POINT(v), GTS_POINT(v1), GTS_POINT(v2));
		}
		else
		{
			fprintf(stderr, "#! @cpt_vertex_pseudonormal: Something is *REALLY* wrong!\n");
			return;
		}
	
		// Triangle normal
		gts_triangle_normal (t, &normalFace[0], &normalFace[1], &normalFace[2]);
		gts_vector_normalize(normalFace);

		normalPseudo[0] += normalFace[0] * angle;
		normalPseudo[1] += normalFace[1] * angle;
		normalPseudo[2] += normalFace[2] * angle;
	}
	g_slist_free(triangles);

	gts_vector_normalize(normalPseudo);

	*x = normalPseudo[0];
	*y = normalPseudo[1];
	*z = normalPseudo[2];
}

/******************************************************************************
 * ;
 ******************************************************************************/
void cpt_distance_eulerian_mesh (SignedDistance* pSignedDistance,
	SizeVector idx_begin, SizeVector idx_end, CptFuncDistance funcDistance,
	gpointer paramFuncDist, gdouble inside_outside
)
{
	GtsPoint  mesh_point;
	gdouble   absDist;
	gdouble   signDist;
	gdouble   distance;
	size_t    i, j, k;
	size_t    idxZ, idxY, mesh_idx;

	// loop through the points inside the box;
	for(k = idx_begin[2]; k <= idx_end[2]; k++)
	{
		idxZ = k * (pSignedDistance->size[0] * pSignedDistance->size[1]);
		mesh_point.z = pSignedDistance->coordMin[2] + (k + 0.5)*pSignedDistance->delta[2];
		
		for(j = idx_begin[1]; j <= idx_end[1]; j++)
		{
			idxY = j * pSignedDistance->size[0];
			mesh_point.y = pSignedDistance->coordMin[1] + (j + 0.5)*pSignedDistance->delta[1];
			
			for(i = idx_begin[0]; i <= idx_end[0]; i++)
			{
				mesh_idx = idxZ + idxY + i;
				mesh_point.x = pSignedDistance->coordMin[0] + (i + 0.5)*pSignedDistance->delta[0];

				distance = funcDistance(&mesh_point , paramFuncDist);
				absDist = fabs(distance);
				
				if( absDist < pSignedDistance->distCut )
				{
					if( absDist < fabs(pSignedDistance->value[ mesh_idx ]) )
					{
						signDist = (inside_outside > 0.0 ? distance / absDist : inside_outside);
						if(absDist > pSignedDistance->sigma) 
							absDist = pSignedDistance->sigma;
							
//						if( signDist * pSignedDistance->value[ mesh_idx ] < 0.0 )
//							printf("wrong\n");
							
						pSignedDistance->value[ mesh_idx ] = signDist * absDist;
					}
				}
			}
		}
	}
}


/******************************************************************************
 * ;
 ******************************************************************************/
gdouble cpt_distance_sign(GtsPoint* pPointMesh, ParamDistFunc* pParam)
{
	GtsVector*  pNormal       = (GtsVector*)pParam->pNormal;
	GtsPoint*   pSurfacePoint = (GtsPoint*)pParam->pSurfacePoint;
	gdouble     distance;
	
	distance = (*pNormal)[0] * (pPointMesh->x - pSurfacePoint->x)
		     + (*pNormal)[1] * (pPointMesh->y - pSurfacePoint->y)
		     + (*pNormal)[2] * (pPointMesh->z - pSurfacePoint->z);
		     
	return distance / fabs(distance);
}

/******************************************************************************
 * ;
 ******************************************************************************/
gdouble cpt_signed_distance_face(GtsPoint* pPointMesh, ParamDistFunc* pParam)
{
	gdouble    distance;
	gdouble    signDist;
	
	distance = gts_point_triangle_distance(pPointMesh, (GtsTriangle*)pParam->pObj);
	signDist = cpt_distance_sign(pPointMesh, pParam);
		     
	return signDist * distance;
}

/******************************************************************************
 * ;
 ******************************************************************************/
gdouble cpt_signed_distance_edge(GtsPoint* pPointMesh, ParamDistFunc* pParam)
{
	gdouble    distance;
	gdouble    signDist;
	
	distance = gts_point_triangle_distance(pPointMesh, (GtsTriangle*)pParam->pObj);
	signDist = cpt_distance_sign(pPointMesh, pParam);
		     
	return signDist * distance;
}


/******************************************************************************
 * ;
 ******************************************************************************/
gdouble cpt_signed_distance_vertex(GtsPoint* pPointMesh, ParamDistFunc* pParam)
{
	gdouble    distance;
	gdouble    signDist;
	
	distance = gts_point_distance(pPointMesh, (GtsPoint*)pParam->pObj);
	signDist = cpt_distance_sign(pPointMesh, pParam);
		     
	return signDist * distance;
}

/******************************************************************************
 * ;
 ******************************************************************************/
// loop at dimensions to build the "index" box;
int calc_idx_box(SizeVector idx_begin, SizeVector idx_end, GtsVector box_begin, GtsVector box_end, SignedDistance* pSignedDistance)
{
	for ( int i = 0 ; i < 3 ; i++ )
	{
		// distance in delta unit;
		idx_begin[i] = (size_t)floor( fabs(box_begin[i] - pSignedDistance->coordMin[i]) / pSignedDistance->delta[i] );
		idx_end[i] =  idx_begin[i] + (size_t)ceil( fabs(box_end[i] - box_begin[i]) / pSignedDistance->delta[i] );

		// Expand by one the box;
		if(idx_begin[i] > 1) idx_begin[i]--;
		if(idx_end[i] < pSignedDistance->size[i]-1) idx_end[i]++;
	}

	return 0;
}

/******************************************************************************
 * ;
 ******************************************************************************/
int cpt_init_distance_function(SignedDistance* pSignedDistance)
{
	GtsPoint  mesh_point;
	size_t  i, j, k;
	size_t  posZ, posY, mesh_pos;
	
	GNode *tree = gts_bb_tree_surface(pSignedDistance->pSurface);
	gboolean inside, open = FALSE;
		
	// loop through the points inside the box;
	for(k = 0; k < pSignedDistance->size[2]; k++)
	{
		posZ = k * (pSignedDistance->size[0] * pSignedDistance->size[1]);
		mesh_point.z = pSignedDistance->coordMin[2] + (k + 0.5)*pSignedDistance->delta[2];

		for(j = 0; j < pSignedDistance->size[1]; j++)
		{
			posY = j * pSignedDistance->size[0];
			mesh_point.y = pSignedDistance->coordMin[1] + (j + 0.5)*pSignedDistance->delta[1];
	
			for(i = 0; i < pSignedDistance->size[0]; i++)
			{
				mesh_pos = posZ + posY + i;
				mesh_point.x = pSignedDistance->coordMin[0] + (i + 0.5)*pSignedDistance->delta[0];
		
#if 0
				if(mesh_point.x * mesh_point.x + mesh_point.y * mesh_point.y + mesh_point.z * mesh_point.z < 1.0)
					pSignedDistance->value[ mesh_pos ] = -pSignedDistance->distMax;
				else
					pSignedDistance->value[ mesh_pos ] = pSignedDistance->distMax;
#endif					
					
				inside = gts_point_is_inside_surface(&mesh_point, tree, open);
				pSignedDistance->value[ mesh_pos ] = pSignedDistance->distMax - inside * 2.0 * pSignedDistance->distMax;
			}
		}
	}
	
	gts_bb_tree_destroy (tree, TRUE);
					
	return 0;
}

/******************************************************************************
 * ;
 ******************************************************************************/
gdouble cpt_round_value(gdouble valueParam, gdouble valueMin)
{
	gdouble  valueRet = valueParam;
	int      i;
	
	while(valueRet > 1.0)
		valueRet /= 10.0;
		
	// Assert delta not too small, otherwise an overflow may occurs at mesh size;
	for(i = 0; valueRet < 1.0; i++)
		valueRet *= 10.0;
		
	valueRet = round(valueRet);
	
	while(i-- > 0 && (valueRet / 10.0) >= valueMin)
		valueRet /= 10.0;
		
	return valueRet;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/******************************************************************************
 ******************************************************************************/
void cpt_face(GtsTriangle *pTriangle, GtsVector normal, SignedDistance* pSignedDistance, gdouble inside_outside)
{
	GtsVertex    *triangle_vertex[3];
	GtsPoint     *triangle_point;
	GtsVector    box_begin = { pSignedDistance->size[0] , pSignedDistance->size[1] , pSignedDistance->size[2] }; 
	GtsVector    box_end   = { 0 , 0 , 0};
	SizeVector   idx_begin, idx_end;
	size_t      i;

	gts_triangle_vertices (pTriangle, &triangle_vertex[0], &triangle_vertex[1], &triangle_vertex[2]);

	// loop at points to build the prism;
	for	( i = 0 ; i < 3 ; i++)
	{
		triangle_point = GTS_POINT (triangle_vertex[i]);
		BOX_BEGIN(box_begin, triangle_point, normal, inside_outside);
		BOX_END  (box_end  , triangle_point, normal, inside_outside);
	}

	calc_idx_box(idx_begin, idx_end, box_begin, box_end, pSignedDistance);

	ParamDistFunc sParam;
	sParam.pNormal       = normal;
	sParam.pObj          = pTriangle;
	sParam.pSurfacePoint = triangle_point;

	cpt_distance_eulerian_mesh(pSignedDistance, idx_begin, idx_end,
		(CptFuncDistance)cpt_signed_distance_face, &sParam, 1.0);
//		(CptFuncDistance)gts_point_triangle_distance, t, 1.0);
}

/******************************************************************************
 * @EULERIAN
 ******************************************************************************/
void cpt_foreach_face(GtsTriangle *t, gpointer pEulerian)
{
	SignedDistance  *pSignedDistance = (SignedDistance*)pEulerian;
	GtsVector     normal;
	int           i;
	static int   idx = 0;

#ifdef CPT_VERBOSE_DEBUG
	printf("# \n");
	printf("# Face #%d:\n", idx);
#endif

	// Computes the coordinates of the oriented normal of @t as the
	// cross-product of two edges, using the left-hand rule. The normal is
	// not normalized.  If this triangle is part of a closed and oriented
	// surface, the normal points to the outside of the surface.
	// void gts_triangle_normal (GtsTriangle * t, gdouble * x, gdouble * y, gdouble * z)
	gts_triangle_normal (t, &normal[0], &normal[1], &normal[2]);
	gts_vector_normalize(normal);

	for ( i = 0 ; i < 3 ; i++ )
		normal[i] *= sqrt(pSignedDistance->distMax);

	cpt_face(t, normal, pSignedDistance,  1.0);
	cpt_face(t, normal, pSignedDistance, -1.0);

	idx++;
}

/******************************************************************************
 * Use the inside_outside parameter to choose between the outside (1.0) and inside (-1.0) normal/prism;
 ******************************************************************************/
void cpt_edge(GtsEdge *e, GtsVector normal[], SignedDistance* pSignedDistance, gdouble inside_outside)
{
	GtsSegment   *edge_segment;
	GtsPoint     *edge_point;
	GtsVector    box_begin = { pSignedDistance->size[0] , pSignedDistance->size[1] , pSignedDistance->size[2] }; 
	GtsVector    box_end   = { 0 , 0 , 0};
	SizeVector   idx_begin, idx_end;


	edge_segment = GTS_SEGMENT(e);

	edge_point = GTS_POINT(edge_segment->v1);	
	BOX_BEGIN(box_begin, edge_point, normal[0], inside_outside);
	BOX_END  (box_end  , edge_point, normal[0], inside_outside);

	edge_point = GTS_POINT(edge_segment->v2);	
	BOX_BEGIN(box_begin, edge_point, normal[1], inside_outside);
	BOX_END  (box_end  , edge_point, normal[1], inside_outside);

#ifdef CPT_VERBOSE_DEBUG
	for ( int i = 0 ; i < 3 ; i++ )
		printf("#\tbox { %g , %g };\n", box_begin[i], box_end[i]);
#endif

	calc_idx_box(idx_begin, idx_end, box_begin, box_end, pSignedDistance);

	cpt_distance_eulerian_mesh(pSignedDistance, idx_begin, idx_end,
		(CptFuncDistance)gts_point_segment_distance, edge_segment, inside_outside);

}

/******************************************************************************
 * @EULERIAN
 ******************************************************************************/
void cpt_foreach_edge(GtsEdge *e, gpointer pEulerian)
{
	SignedDistance  *pSignedDistance = (SignedDistance*)pEulerian;
	GtsVector     normal[2];
	GtsTriangle   *t;
	GSList        *item          = NULL;
	gdouble       inside_outside = 1.0;
	gdouble       theta;
	int           i , j;
	static int   idx = 0;


#ifdef CPT_VERBOSE_DEBUG
	printf("# \n");
	printf("# Edge #%d:\n", idx);
#endif

	// walk by the triangle list
	for(item = e->triangles , j = 0 ; item; item = item->next , j++)
	{
		if (j > 1)
		{
			// TODO : Mensagem de erro: Superficie nao eh valida para aplicacao do CPT; (AINDA!)
			// Exemplo de caso: Estrela
			return;
		}
		t = GTS_TRIANGLE(item->data);
		gts_triangle_normal (t, &normal[j][0], &normal[j][1], &normal[j][2]);
		gts_vector_normalize(normal[j]);
	}

	theta = cpt_vector_angle(normal[0], normal[1]);

	// Check if we have it as inside or outside;
	if ( theta < 0.0)
	{
		// theta is updated, but its useless for the algorithm;
		theta += 2.*M_PI;
		inside_outside = -1.0;
	}

	// Expand normals to create the "canaleta" with the correct maximum size;
	for ( j = 0 ; j < 2 ; j++ )
		for ( i = 0 ; i < 3 ; i++ )
			normal[j][i] *= sqrt(pSignedDistance->distMax);

	cpt_edge(e, normal, pSignedDistance, inside_outside);

	idx++;
}

/******************************************************************************
 *
 ******************************************************************************/
void cpt_vertex(GtsVertex *v, GtsVector normal, SignedDistance* pSignedDistance)
{
	ParamDistFunc sParam;
	GtsPoint     *point = GTS_POINT(v);
	GtsVector    box_begin = { pSignedDistance->size[0] , pSignedDistance->size[1] , pSignedDistance->size[2] }; 
	GtsVector    box_end   = { 0 , 0 , 0};
	SizeVector   idx_begin, idx_end;

	BOX_BEGIN(box_begin, point, normal,  1.0);
	BOX_END  (box_end  , point, normal,  1.0);
	BOX_BEGIN(box_begin, point, normal, -1.0);
	BOX_END  (box_end  , point, normal, -1.0);

#ifdef CPT_VERBOSE_DEBUG
	printf("#\n");
	for ( i = 0 ; i < 3 ; i++ )
	{
		printf("#\tbox { %g , %g };\n", box_begin[i], box_end[i]);
	}
#endif

	// loop at dimensions to build the box;
	calc_idx_box(idx_begin, idx_end, box_begin, box_end, pSignedDistance);

	sParam.pNormal = normal;
	sParam.pObj    = point;

	cpt_distance_eulerian_mesh(pSignedDistance, idx_begin, idx_end,
		(CptFuncDistance)cpt_signed_distance_vertex, &sParam, 1.0);
}

/******************************************************************************
 * @EULERIAN
 ******************************************************************************/
void cpt_foreach_vertex(GtsVertex *v, gpointer pEulerian)
{
	SignedDistance  *pSignedDistance = (SignedDistance*)pEulerian;
	GtsVector     normal;
	static int    idx = 0;


#ifdef CPT_VERBOSE_DEBUG
	printf("# \n");
	printf("# Vertex #%d:\n", idx);
#endif

	cpt_vertex_pseudonormal(v, &normal[0], &normal[1], &normal[2]);

	for ( int i = 0 ; i < 3 ; i++ )
		normal[i] *= sqrt(pSignedDistance->distMax);
		
	cpt_vertex(v, normal, pSignedDistance);

	idx++;
}

/******************************************************************************
 *
 ******************************************************************************/
void surface_findbox(GtsVertex *v, gpointer pEulerian)
{
	g_return_if_fail (v != NULL);
	g_return_if_fail (pEulerian != NULL);

	SignedDistance  *pSignedDistance = (SignedDistance*)pEulerian;

	try {

		if(pSignedDistance->coordMin[0] > v->p.x)
			pSignedDistance->coordMin[0] = v->p.x;

		if(pSignedDistance->coordMin[1] > v->p.y)
			pSignedDistance->coordMin[1] = v->p.y;

		if(pSignedDistance->coordMin[2] > v->p.z)
			pSignedDistance->coordMin[2] = v->p.z;

		if(pSignedDistance->coordMax[0] < v->p.x)
			pSignedDistance->coordMax[0] = v->p.x;

		if(pSignedDistance->coordMax[1] < v->p.y)
			pSignedDistance->coordMax[1] = v->p.y;

		if(pSignedDistance->coordMax[2] < v->p.z)
			pSignedDistance->coordMax[2] = v->p.z;

	}
	catch(...) {

		fprintf(stderr, ">>> Exception at surface_findbox;");
		return;

	}
}

/******************************************************************************
 *
 ******************************************************************************/
size_t cpt_eulerian_mesh_write_silo(SignedDistance *pSignedDistance,  const char* szFileOut)
{
	static const char *coordnames[] = { "X" , "Y" , "Z" };    // Name the coordinate axes 'X' and 'Y'
	DBfile  *dbfile= NULL;            // The Silo file pointer
	float   *nodex, *nodey, *nodez;   // The coordinate arrays
	float   *var;                     // Field var
	int     dims[3], size;            // The number of nodes in each dimension
	int     i, j, k;
	int     idxZ, idxY, mesh_idx;
	size_t  valid = 0;


	// How many nodes in each direction?
	dims[0] = pSignedDistance->size[0];
	dims[1] = pSignedDistance->size[1];
	dims[2] = pSignedDistance->size[2];

	nodex = (float*) malloc(dims[0]*sizeof(float));
	nodey = (float*) malloc(dims[1]*sizeof(float));
	nodez = (float*) malloc(dims[2]*sizeof(float));

	for(i = 0; i < dims[0]; i++)	nodex[i] = (float) (pSignedDistance->coordMin[0] + (i + 0.5)*pSignedDistance->delta[0]);
	for(j = 0; j < dims[1]; j++)	nodey[j] = (float) (pSignedDistance->coordMin[1] + (j + 0.5)*pSignedDistance->delta[1]);
	for(k = 0; k < dims[2]; k++)	nodez[k] = (float) (pSignedDistance->coordMin[2] + (k + 0.5)*pSignedDistance->delta[2]);

	// Define the array of coordinate and assign coordinates to coordinates array
	float *coordinates[] = {(float*)nodex, (float*)nodey, (float*)nodez};

	dbfile = DBCreate(szFileOut, DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
	if(dbfile == NULL)
	{
		fprintf(stderr, "Cant create SILO file: '%s'!", szFileOut);
		return 0;
	}

	DBPutQuadmesh(dbfile, "CartesianMesh", (char**)coordnames, coordinates, dims, 3, DB_FLOAT, DB_COLLINEAR, NULL);

	size = dims[0]*dims[1]*dims[2];
	var = (float*) malloc(size*sizeof(float));

	// Field function: colorfunction
	for(k = 0; k < dims[2]; k++)
	{
		idxZ = k * dims[0] * dims[1];
		for(j = 0; j < dims[1]; j++)
		{
			idxY = j * dims[0];
			for(i = 0; i < dims[0]; i++)
			{
				mesh_idx = idxZ + idxY + i;

				//var[mesh_idx] = (float)(pSignedDistance->value[ mesh_idx ] > 0.0 ? 1.0 : -1.0 );
				var[mesh_idx] = (float)(pSignedDistance->value[ mesh_idx ]);
				if( fabs(pSignedDistance->value[ mesh_idx ]) < pSignedDistance->distMax)
				{
					valid++;
				}
			}
		}
	}

	DBPutQuadvar1(dbfile   , "colorfunction", "CartesianMesh", var   , dims, 3, NULL, 0, DB_FLOAT, DB_NODECENT, NULL);

	DBClose(dbfile);

	free(var);

	free(nodex);
	free(nodey);
	free(nodez);

	return valid;
}
/******************************************************************************
 *
 ******************************************************************************/
#if 0
size_t cpt_eulerian_mesh_write_ze(SignedDistance *pSignedDistance, GtsSurface *surface, const char* szFileOut)
{
	static const char *coordnames[] = { "X" , "Y" , "Z" };    // Name the coordinate axes 'X' and 'Y'
	DBfile  *dbfile_ze= NULL;            // The Silo file pointer
	float   *nodex, *nodey, *nodez;   // The coordinate arrays
	float   *var_ze;                     // Field var
	int     dims[3], size;            // The number of nodes in each dimension
	int     i, j, k;
	int     idxZ, idxY, mesh_idx;
	size_t  valid = 0;


	// How many nodes in each direction?
	dims[0] = pSignedDistance->size[0];
	dims[1] = pSignedDistance->size[1];
	dims[2] = pSignedDistance->size[2];

	nodex = (float*) malloc(dims[0]*sizeof(float));
	nodey = (float*) malloc(dims[1]*sizeof(float));
	nodez = (float*) malloc(dims[2]*sizeof(float));

	for(i = 0; i < dims[0]; i++)	nodex[i] = (float) (pSignedDistance->coordMin[0] + (i + 0.5)*pSignedDistance->delta[0]);
	for(j = 0; j < dims[1]; j++)	nodey[j] = (float) (pSignedDistance->coordMin[1] + (j + 0.5)*pSignedDistance->delta[1]);
	for(k = 0; k < dims[2]; k++)	nodez[k] = (float) (pSignedDistance->coordMin[2] + (k + 0.5)*pSignedDistance->delta[2]);

	// Define the array of coordinate and assign coordinates to coordinates array
	float *coordinates[] = {(float*)nodex, (float*)nodey, (float*)nodez};

	dbfile_ze = DBCreate(szFileOut, DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
	if(dbfile_ze == NULL)
	{
		fprintf(stderr, "Cant create SILO file: '%s'!", szFileOut);
		return 0;
	}
	DBPutQuadmesh(dbfile_ze, "CartesianMesh", (char**)coordnames, coordinates, dims, 3, DB_FLOAT, DB_COLLINEAR, NULL);

	// The data must be (NX-1) * (NY-1) * (NZ-1), since colorfunction is in cell center.
	dims[0]--; dims[1]--; dims[2]--;

	size = dims[0]*dims[1]*dims[2];
	var_ze = (float*) malloc(size*sizeof(float));

	/* ---------------- jeso.Test  ---------------- */
	GNode *tree = gts_bb_tree_surface(surface);
	gboolean inside, open = FALSE;
	GtsPoint point;
	/* ---------------- jeso.Test  ---------------- */

	// Field function: colorfunction
	for(k = 0; k < dims[2]; k++)
	{
		idxZ = k * dims[0] * dims[1];
		for(j = 0; j < dims[1]; j++)
		{
			idxY = j * dims[0];
			for(i = 0; i < dims[0]; i++)
			{
				mesh_idx = idxZ + idxY + i;

				/* ---------------- jeso.Test  ---------------- */
				gts_point_set(&point, nodex[i], nodey[j], nodez[k]);
				inside = gts_point_is_inside_surface(&point, tree, open);
				if (inside) {
					var_ze[mesh_idx] = 0.0;
					valid++;
				} else {
					var_ze[mesh_idx] = 1.0;
				}
				/* ---------------- jeso.Test  ---------------- */
			}
		}
	}

	DBPutQuadvar1(dbfile_ze, "colorfunction", "CartesianMesh", var_ze, dims, 3, NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);

	DBClose(dbfile_ze);

	free(var_ze);

	free(nodex);
	free(nodey);
	free(nodez);

	return valid;
}
#endif

// EoF

