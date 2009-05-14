// export.cpp

#include "gtscpt.h"
#include "cpt.h"
#include <silo.h>
#include <VisItControlInterface_V1.h>

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

