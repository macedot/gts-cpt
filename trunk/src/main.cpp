// main.cpp
// version 1.0 : Thiago Macedo
// for use with c++ code see : http://geneura.ugr.es/~jmerelo/c++-faq/mixing-c-and-cpp.html
// base code structure from Jose Eduardo;

/*
TODO:
- Dominio (omega) [A1, B1] x [A2, B2] x [A3, B3] fornecido externamente;
- encaixar box da superficie na malha euleriana que JAH EXISTE!!!
*/

#include "gtscpt.h"
#include "cpt.h"

#define CHECK_DELTA

// "The quality of a triangle is defined as the ratio of its surface
// to its perimeter relative to this same ratio for an equilateral
// triangle with the same area. The quality is then one for equilateral
// triangle and tends to zero for a very stretched triangle."
static const gdouble  DBL_CLOCKS_PER_SEC  = (double)CLOCKS_PER_SEC;
       const gdouble  MIN_TRANGLE_QUALITY = 0.5;   // TODO: Quanto eh bom?
       const gdouble  MIN_DELTA_MESH      = 0.001;
       const gdouble  size_sup            = 2.0;   // suporte da delta de dirac;
static      gboolean  verbose             = TRUE;

int main(int argc, char *argv[])
{
	SignedDistance  *pSignedDistance = NULL;
//	GtsFile     *fp = NULL;
	time_t      theClock;
	time_t      timeCpt[3];
	size_t     i;

	try {

		if(!setlocale(LC_ALL, "POSIX"))
			g_warning("cannot set locale to POSIX");

		////////////////////////////////////////////////////////////////////////
			
		// @TODO
		pSignedDistance = new SignedDistance;
		memset(pSignedDistance, 0, sizeof(SignedDistance));

		for(i = 0; i < 3; i++)
		{
			pSignedDistance->coordMin[i] =  DBL_MAX;
			pSignedDistance->coordMax[i] = -DBL_MAX;
		}

		////////////////////////////////////////////////////////////////////////

		pSignedDistance->pSurface = gts_surface_new(gts_surface_class(), gts_face_class(), gts_edge_class(), gts_vertex_class());
#if 0
		fp = gts_file_new(stdin);
		if(gts_surface_read(pSignedDistance->pSurface, fp))
		{
			fputs("gtscpt: the file on standard input is not a valid GTS file\n", stderr);
			fprintf(stderr, "stdin:%d:%d: %pSignedDistance->pSurface\n", fp->line, fp->pos, fp->error);
			return 1; // failure
		}
#endif

		guint geodesation_order = 7;

		if(verbose)
		{
			printf("#\n");
			printf("# Generation unit sphere (geodesation_order = %d)...", geodesation_order);

			fflush(stdout);
			fflush(stderr);
		}

		theClock = clock();
		gts_surface_generate_sphere(pSignedDistance->pSurface, geodesation_order);
		theClock = (clock() - theClock);

		if(verbose)
		{
			printf(" done [%+1.8f]\n", theClock / DBL_CLOCKS_PER_SEC);
			printf("#\n");
		}

		fflush(stdout);
		fflush(stderr);

		////////////////////////////////////////////////////////////////////////

		if(verbose)
		{
			gts_surface_print_stats(pSignedDistance->pSurface, stdout);
			printf("#   Total area: %g\n", gts_surface_area(pSignedDistance->pSurface));
		}

		////////////////////////////////////////////////////////////////////////
		
		// test if surface is a closed and orientable manifold.
		// we don't need to test if pSignedDistance->pSurface is a manifold since both tests below implies that.
		if(!gts_surface_is_closed(pSignedDistance->pSurface))
		{
			fprintf(stderr, ">>> WARNING: surface is not closed -> not a manifold;\n");
			return 1;
		}
		if(!gts_surface_is_orientable(pSignedDistance->pSurface))
		{
			fprintf(stderr, ">>> WARNING: surface is not orientable -> not a manifold;\n");
			return 1;
		}

		GtsSurfaceQualityStats   qstats;
		gts_surface_quality_stats(pSignedDistance->pSurface, &qstats);

		// "The quality of a triangle is defined as the ratio of its surface
		// to its perimeter relative to this same ratio for an equilateral
		// triangle with the same area. The quality is then one for equilateral
		// triangle and tends to zero for a very stretched triangle."
		if(qstats.face_quality.max < MIN_TRANGLE_QUALITY)
		{
			fprintf(stderr, ">>> FATAL: global trangle quality ratio below minimum triangle quality parameter;\n");
			return 1;
		}
		else if(qstats.face_quality.mean < MIN_TRANGLE_QUALITY)
		{
			fprintf(stderr, ">>> WARNING: mean of trangle quality ratio below minimum triangle quality parameter;\n");
		}
		else if(qstats.face_quality.min < MIN_TRANGLE_QUALITY)
		{
			fprintf(stderr, ">>> WARNING: some trangle have quality ratio below minimum triangle quality parameter;\n");
		}

		// @TODO:Definir funcoes usadas como parametro;
		// @TODO: Melhorar o objeto discretizado tal que:
		//	aresta media menor que deltas da malha euleriana (faixa para tamanhos de arestas);
		// gts_surface_refine(pSignedDistance->pSurface, NULL, NULL, NULL, NULL, refine_stop, refine_stop_param);

		////////////////////////////////////////////////////////////////////////

		// Use the minimum valid delta;
		pSignedDistance->delta[0] = 1.1 * (qstats.edge_length.mean + 2.0*qstats.edge_length.stddev);
		pSignedDistance->delta[0] = cpt_round_value(pSignedDistance->delta[0], MIN_DELTA_MESH);

		pSignedDistance->delta[1] = (1.0 + 0.05) * pSignedDistance->delta[0];
		pSignedDistance->delta[2] = (1.0 - 0.05) * pSignedDistance->delta[0];
		
		//pSignedDistance->sigma = 3.0 * cpt_min(pSignedDistance->delta[0], cpt_min(pSignedDistance->delta[1], pSignedDistance->delta[2]));
		pSignedDistance->sigma = 2.0 * cpt_min(pSignedDistance->delta[0], cpt_min(pSignedDistance->delta[1], pSignedDistance->delta[2]));

		gdouble  delta_max = cpt_max(pSignedDistance->delta[0], cpt_max(pSignedDistance->delta[1], pSignedDistance->delta[2]));
		pSignedDistance->distCut = cpt_get_dist_cut(size_sup, delta_max);
//		pSignedDistance->distCut = cpt_round_value(pSignedDistance->distCut, FLT_EPSILON);
		pSignedDistance->distMax = cpt_get_dist_max(pSignedDistance->distCut, delta_max);

		if(verbose)
		{
			printf("#\tdistCut = { %+1.8f }; \n", pSignedDistance->distCut);
			printf("#\tdistMax = { %+1.8f }; \n", pSignedDistance->distMax);
			printf("#\tsigma   = { %+1.8f }; \n", pSignedDistance->sigma);
		}
		
		////////////////////////////////////////////////////////////////////////

		fflush(stdout);
		fflush(stderr);

		// find the box than contains the surface;
		theClock = clock();
		gts_surface_foreach_vertex(pSignedDistance->pSurface, (GtsFunc)surface_findbox, pSignedDistance);
		theClock = (clock() - theClock);

		if(verbose)
		{
			printf("#\n");
			printf("# Finding box...\n");
			printf("#\tCallback process time: [%+1.8f]\n", theClock / DBL_CLOCKS_PER_SEC);
			printf("#\tcoordMin = { %+1.8f , %+1.8f , %+1.8f}; \n", pSignedDistance->coordMin[0], pSignedDistance->coordMin[1], pSignedDistance->coordMin[2]);
			printf("#\tcoordMax = { %+1.8f , %+1.8f , %+1.8f}; \n", pSignedDistance->coordMax[0], pSignedDistance->coordMax[1], pSignedDistance->coordMax[2]);
			printf("#\tdelta    = { %+1.8f , %+1.8f , %+1.8f}; \n", pSignedDistance->delta[0]   , pSignedDistance->delta[1]   , pSignedDistance->delta[2]   );
		}

		////////////////////////////////////////////////////////////////////////

		// @TODO : encaixar box na malha euleriana;
		
		// POR HORA FAREMOS COMO ABAIXO;

		gdouble size_box;

		// round coords to create a "big" box containg the surface;
		for(i = 0; i < 3; i++)
		{
			// make sure we have some space between surface and the box boundary;
			size_box = fabs(pSignedDistance->coordMax[i] - pSignedDistance->coordMin[i]);
			pSignedDistance->coordMin[i]-= size_box;
			pSignedDistance->coordMax[i]+= size_box;
			// ... but assert the maximun coord build a box with
			// a dimention multiple of the delta in this axis;
			pSignedDistance->size[i]     = (size_t)ceil(fabs(pSignedDistance->coordMax[i] - pSignedDistance->coordMin[i]) / pSignedDistance->delta[i]);
			pSignedDistance->coordMax[i] = pSignedDistance->coordMin[i] + pSignedDistance->size[i] * pSignedDistance->delta[i];
		}

		if(verbose)
		{
			printf("#\n");
			printf("# Expanding box...\n");
			printf("#\tcoordMin = { %+1.8f , %+1.8f , %+1.8f}; \n", pSignedDistance->coordMin[0], pSignedDistance->coordMin[1], pSignedDistance->coordMin[2]);
			printf("#\tcoordMax = { %+1.8f , %+1.8f , %+1.8f}; \n", pSignedDistance->coordMax[0], pSignedDistance->coordMax[1], pSignedDistance->coordMax[2]);
			printf("#\tsize     = {   %6lu ,   %6lu ,   %6lu}; \n", pSignedDistance->size[0]    , pSignedDistance->size[1]    , pSignedDistance->size[2]    );
		}

		////////////////////////////////////////////////////////////////////////

		if(pSignedDistance->size[0] < 1 ||
		   pSignedDistance->size[1] < 1 ||
		   pSignedDistance->size[2] < 1)
		{
			g_warning("INVALID MESH SIZE!!");
			return -1;
		}

		size_t  mesh_size = (pSignedDistance->size[0] * pSignedDistance->size[1] * pSignedDistance->size[2]);

		if(verbose)
		{
			printf("#\n");
			printf("# Creating the eulerian mesh [size = %lu]...", mesh_size);
		}

		fflush(stdout);
		fflush(stderr);

		pSignedDistance->value = new gdouble[mesh_size];

		theClock = clock();
		cpt_init_distance_function(pSignedDistance);
		theClock = (clock() - theClock);

		if(verbose)
		{
			printf("done [%+1.8f];\n", theClock / DBL_CLOCKS_PER_SEC);
			printf("#\n");
		}

		fflush(stdout);
		fflush(stderr);

		////////////////////////////////////////////////////////////////////////
		static const int maxLoop = 1;
		size_t qt_pts;

		timeCpt[0] = timeCpt[1] = timeCpt[2] = 0;

		if(verbose)
			printf("#\n");

		for ( int t = 0 ; t < maxLoop; t++)
		{
			theClock = clock();
			gts_surface_foreach_face(pSignedDistance->pSurface,(GtsFunc) cpt_foreach_face, pSignedDistance);
			theClock = (clock() - theClock);
			timeCpt[0] += theClock;

			printf("# gts_surface_foreach_face: [%+1.8f sec];\n", timeCpt[0] / (DBL_CLOCKS_PER_SEC * maxLoop));

			theClock = clock();
			gts_surface_foreach_edge(pSignedDistance->pSurface,(GtsFunc) cpt_foreach_edge, pSignedDistance);
			theClock = (clock() - theClock);
			timeCpt[1] += theClock;

			printf("# gts_surface_foreach_edge: [%+1.8f sec];\n", timeCpt[1] / (DBL_CLOCKS_PER_SEC * maxLoop));

			theClock = clock();
			gts_surface_foreach_vertex(pSignedDistance->pSurface,(GtsFunc) cpt_foreach_vertex, pSignedDistance);
			theClock = (clock() - theClock);
			timeCpt[2] += theClock;

			printf("# gts_surface_foreach_vertex: [%+1.8f sec];\n", timeCpt[2] / (DBL_CLOCKS_PER_SEC * maxLoop));
		}

		if(verbose)
		{
			printf("#\n");
			printf("# total time: [%+1.8f sec];\n", (timeCpt[0] + timeCpt[1] + timeCpt[2]) / DBL_CLOCKS_PER_SEC);
			if(maxLoop > 1)
				printf("# avg time: [%+1.8f sec];\n", (timeCpt[0] + timeCpt[1] + timeCpt[2]) / (DBL_CLOCKS_PER_SEC * maxLoop));
		}

		fflush(stdout);
		fflush(stderr);

		////////////////////////////////////////////////////////////////////////

		theClock = clock();
		qt_pts = cpt_eulerian_mesh_write_silo(pSignedDistance, "out/eulerianmesh.silo");
		theClock = (clock() - theClock);
		if(verbose)
		{
			printf("#\n");
			printf("# Total points in the interface: [%lu ; %+1.8f sec];\n", qt_pts, theClock / DBL_CLOCKS_PER_SEC);
		}

#if 0
		theClock = clock();
		qt_pts = cpt_eulerian_mesh_write_ze(pSignedDistance, pSignedDistance->pSurface, "out/eulerianmesh_ze.silo");
		theClock = (clock() - theClock);
		if(verbose)
		{
			printf("#\n");
			printf("# Total points inside surface: [%lu ; %+1.8f sec];\n", qt_pts, theClock / DBL_CLOCKS_PER_SEC);
		}
#endif
		////////////////////////////////////////////////////////////////////////
#if 0
		printf("#\n");
		for(mesh_pos = 0; mesh_pos < mesh_size; mesh_pos++)
		{
			printf("#\t%lu = %+1.8f\n", mesh_pos, pSignedDistance->value[mesh_pos]);
		}
#endif
		////////////////////////////////////////////////////////////////////////

		FILE* fpVtk = fopen("out/surface.vtk", "w+b");
		gts_surface_write_vtk(pSignedDistance->pSurface, fpVtk);
		fclose(fpVtk);

		////////////////////////////////////////////////////////////////////////

		gts_object_destroy(GTS_OBJECT(pSignedDistance->pSurface));
		gts_finalize();

		delete [] pSignedDistance->value;  pSignedDistance->value = NULL;
		delete pSignedDistance;           pSignedDistance       = NULL;

		fflush(stdout);
		fflush(stderr);

	}
	catch(...) {

		g_warning("unhandle exception");
		return -1;

	}

	return 0;
}
// EoF
