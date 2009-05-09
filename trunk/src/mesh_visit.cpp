// mesh_visit.cpp
// version 1.0 : Thiago Macedo
// Interface to send eulerian mesh to VisIt;

/* Simulation mesh */
float mesh_x[] = {0., 1., 2.5, 5.};
float mesh_y[] = {0., 2., 2.25, 2.55,   5.};
int    mesh_dims[] = {4, 5, 1};
int    mesh_ndims = 2;

VisIt_MeshData *VisItGetMesh(int domain, const char *name)
{
	VisIt_MeshData *mesh = NULL;
	size_t sz = sizeof(VisIt_MeshData);
	
	if(strcmp(name, "mesh2d") == 0)
	{
		/* Allocate VisIt_MeshData. */
		mesh = (VisIt_MeshData *)malloc(sz);
		memset(mesh, 0, sz);
		
		/* Make VisIt_MeshData contain a VisIt_RectilinearMesh. */
		sz = sizeof(VisIt_RectilinearMesh);
		mesh->rmesh = (VisIt_RectilinearMesh *)malloc(sz);
		memset(mesh->rmesh, 0, sz);
		
		/* Tell VisIt which mesh object to use. */
		mesh->meshType = VISIT_MESHTYPE_RECTILINEAR;
		
		/* Set the mesh’s number of dimensions. */
		mesh->rmesh->ndims = mesh_ndims;
		
		/* Set the mesh dimensions. */
		mesh->rmesh->dims[0] = mesh_dims[0];
		mesh->rmesh->dims[1] = mesh_dims[1];
		mesh->rmesh->dims[2] = mesh_dims[2];
		mesh->rmesh->baseIndex[0] = 0;
		mesh->rmesh->baseIndex[1] = 0;
		mesh->rmesh->baseIndex[2] = 0;
		mesh->rmesh->minRealIndex[0] = 0;
		mesh->rmesh->minRealIndex[1] = 0;
		mesh->rmesh->minRealIndex[2] = 0;
		mesh->rmesh->maxRealIndex[0] = mesh_dims[0]-1;
		mesh->rmesh->maxRealIndex[1] = mesh_dims[1]-1;
		mesh->rmesh->maxRealIndex[2] = mesh_dims[2]-1;
		
		/* Let VisIt use simulation’s copy of the mesh coordinates. */
		mesh->rmesh->xcoords = VisIt_CreateDataArrayFromFloat(VISIT_OWNER_SIM, mesh_x);
		mesh->rmesh->ycoords = VisIt_CreateDataArrayFromFloat(VISIT_OWNER_SIM, mesh_y);
	}
	
	return mesh;
}

#define NPTS 100
float  angle = 0.;
int    pmesh_ndims = 3;
float  pmesh_x[NPTS], pmesh_y[NPTS], pmesh_z[NPTS];

void simulate_one_timestep(void)
{
	int i;
	
	for(i = 0; i < NPTS; ++i)
	{
		float t = ((float)i) / ((float)(NPTS-1));
		float a = 3.14159 * 10. * t;
		pmesh_x[i] = t * cos(a + (0.5 + 0.5 * t) *angle);
		pmesh_y[i] = t * sin(a + (0.5 + 0.5 * t) * angle);
		pmesh_z[i] = t;
	}
	angle = angle + 0.05;
}

VisIt_MeshData *VisItGetMesh(int domain, const char *name)
{
	VisIt_MeshData *mesh = NULL;
	size_t sz = sizeof(VisIt_MeshData);

	if(strcmp(name, "point3d") == 0)
	{
		/* Allocate VisIt_MeshData. */
		mesh = (VisIt_MeshData *)malloc(sz);
		memset(mesh, 0, sz);
		
		/* Make VisIt_MeshData contain a VisIt_PointMesh. */
		sz = sizeof(VisIt_PointMesh);
		mesh->pmesh = (VisIt_PointMesh *)malloc(sz);
		memset(mesh->pmesh, 0, sz);
		
		/* Tell VisIt which mesh object to use. */
		mesh->meshType = VISIT_MESHTYPE_POINT;
		
		/* Set the mesh’s number of dimensions. */
		mesh->pmesh->ndims = pmesh_ndims;
		
		/* Set the number of points in the mesh. */
		mesh->pmesh->nnodes = NPTS;
		
		/* Let VisIt use simulation’s copy of the mesh coordinates. */
		mesh->pmesh->xcoords = VisIt_CreateDataArrayFromFloat(VISIT_OWNER_SIM, (float *)pmesh_x);
		mesh->pmesh->ycoords = VisIt_CreateDataArrayFromFloat(VISIT_OWNER_SIM, (float *)pmesh_y);
		mesh->pmesh->zcoords = VisIt_CreateDataArrayFromFloat(VISIT_OWNER_SIM, (float *)pmesh_z);
	}
	return mesh;
}


// EoF
