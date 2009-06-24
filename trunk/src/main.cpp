// main.cpp
// version 1.0 : Thiago Macedo
// for use with c++ code see : http://geneura.ugr.es/~jmerelo/c++-faq/mixing-c-and-cpp.html
// base code structure from Jose Eduardo;

#include "gtscpt.h"
#include "cpt.h"
#include "export.h"
#include <getopt.h>

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif /* HAVE_UNISTD_H */

#define CHECK_DELTA

// If defined, use an runtime generated sphere as input example wit the given geodesical order;
// if > 8 -> visit cant plot pseudocolor (mem > 4GiB (?) !!)
//#define INPUT_RUNTIME_SPHERE    8

// "The quality of a triangle is defined as the ratio of its surface
// to its perimeter relative to this same ratio for an equilateral
// triangle with the same area. The quality is then one for equilateral
// triangle and tends to zero for a very stretched triangle."
static const gdouble  DBL_CLOCKS_PER_SEC  = (double)CLOCKS_PER_SEC;
       const gdouble  MIN_TRANGLE_QUALITY = 0.5;   // TODO: Quanto eh bom?
       const gdouble  MIN_DELTA_MESH      = 0.001;
       const gdouble  size_sup            = 2.0;   // suporte da delta de dirac;
static      gboolean  verbose             = TRUE;
static      gboolean  normalize           = FALSE;
static         guint  geodesation_order   = 0;

int main(int argc, char *argv[])
{
	SignedDistance  *pSignedDistance = NULL;
	time_t      theClock;
	time_t      timeCpt[3];
	GtsFile*    fp;
	GtsVector   beginMesh = { 0.0 , 0.0 , 0.0 };     // begin of eulerian mesh (x,y,z)
	GtsVector   endMesh   = { 0.0 , 0.0 , 0.0 };     // end of eulerian mesh (x,y,z)
	SizeVector  sizeMesh  = {   0 ,   0 ,   0 };     // site in delta unit;
	size_t      i;
	int         c = 0;

	try {

		if(!setlocale(LC_ALL, "POSIX"))
			g_warning("cannot set locale to POSIX");

		////////////////////////////////////////////////////////////////////////

		static struct option long_options[] = {
			{"begin-x"  , required_argument, NULL, 'A'},
			{"begin-y"  , required_argument, NULL, 'B'},
			{"begin-z"  , required_argument, NULL, 'C'},
			{"end-x"    , required_argument, NULL, 'D'},
			{"end-y"    , required_argument, NULL, 'E'},
			{"end-z"    , required_argument, NULL, 'F'},
			{"size-x"   , required_argument, NULL, 'G'},
			{"size-y"   , required_argument, NULL, 'H'},
			{"size-z"   , required_argument, NULL, 'I'},
			{"normalize", no_argument      , NULL, 'o'},
			{"sphere"   , required_argument, NULL, 's'},
			{"verbose"  , no_argument      , NULL, 'v'},
			{"help"     , no_argument      , NULL, 'h'},
			{ NULL }
		};
		int option_index = 0;

		/* parse options using getopt */
		while (c != EOF) {
			switch ((c = getopt_long (argc, argv, "A:B:C:D:E:F:G:H:I:os:vh", long_options, &option_index))) {
				case 'A': /* x */
					beginMesh[0] = (gdouble) atof (optarg);
					break;
				case 'B': /* sy */
					beginMesh[1] = (gdouble) atof (optarg);
					break;
				case 'C': /* sz */
					beginMesh[2] = (gdouble) atof (optarg);
					break;
				case 'D': /* x */
					endMesh[0] = (gdouble) atof (optarg);
					break;
				case 'E': /* sy */
					endMesh[1] = (gdouble) atof (optarg);
					break;
				case 'F': /* sz */
					endMesh[2] = (gdouble) atof (optarg);
					break;
				case 'G': /* sx */
					sizeMesh[0] = (size_t) atoi (optarg);
					break;
				case 'H': /* sy */
					sizeMesh[1] = (size_t) atoi (optarg);
					break;
				case 'I': /* sz */
					sizeMesh[2] = (size_t) atoi (optarg);
					break;
				case 'o': /* normalize */
					normalize = TRUE;
					break;
				case 's': /* verbose */
					geodesation_order = (guint) atoi (optarg);
					if( geodesation_order <= 0 )
					{
						fputs("gtscpt: you must supply a valid integer great than zero value as geodesation order parameter!\n", stderr);
						return 1; // failure
					}
					break;
				case 'v': /* verbose */
					verbose = TRUE;
					break;
				case 'h': /* help */
				case '?': /* wrong options */
					fprintf (stderr,
					"Usage: gts-cpt [OPTION] < file.gts\n"
					"CPT using the GTS library.\n"
					"\n"

					"  --begin-x VALUE              begin of eulerian mesh (x)\n"
					"  --begin-y VALUE              begin of eulerian mesh (x)\n"
					"  --begin-z VALUE              begin of eulerian mesh (x)\n"

					"  --end-x VALUE                end of eulerian mesh (x)\n"
					"  --end-y VALUE                end of eulerian mesh (x)\n"
					"  --end-z VALUE                end of eulerian mesh (x)\n"

					"  --size-x VALUE               number of cells for X axis\n"
					"  --size-y VALUE               number of cells for Y axis\n"
					"  --size-z VALUE               number of cells for Z axis\n"

					"  --normalize                  fit the resulting surface in a cube of\n"
					"                               size 1 centered at the origin\n"

					"  --verbose                    print statistics about the surface\n"
					"  --help                       display this help and exit\n"
					"\n"
					"Reports bugs to tmacedo@usp.br\n");
					return (c != 'h'); /* success or failure */
				}
		}

		////////////////////////////////////////////////////////////////////////		

		// do not allow 2D eulerian mesh ( no reason in special... just i dont care about it, yet =] );
		if ( beginMesh[0] >= endMesh[0] || beginMesh[1] >= endMesh[1] || beginMesh[2] >= endMesh[2])
		{
			fputs("gtscpt: you must correctly especify the eulerian mesh!\n", stderr);
			return 1; // failure
		}

		// valid number of cells in each eulerian mesh axis is a must!		
		if ( sizeMesh[0] <= 0 || sizeMesh[1] <= 0 || sizeMesh[2] <= 0)
		{
			fputs("gtscpt: you must correctly especify the size of eulerian mesh!\n", stderr);
			return 1; // failure
		}

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

		if(geodesation_order > 0)
		{
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
		}
		else
		{
			fp = gts_file_new(stdin);
			if(gts_surface_read(pSignedDistance->pSurface, fp))
			{
				fputs("gtscpt: the file on standard input is not a valid GTS file\n", stderr);
				fprintf(stderr, "stdin:%d:%d: %p\n", fp->line, fp->pos, fp->error);
				return 1; // failure
			}

		}

		////////////////////////////////////////////////////////////////////////

		if(verbose)
		{
			gts_surface_print_stats(pSignedDistance->pSurface, stdout);
			printf("#   Total area: %g\n", gts_surface_area(pSignedDistance->pSurface));
		}

		if (normalize) 
		{
			GtsBBox * bb = gts_bbox_surface (gts_bbox_class (), pSignedDistance->pSurface);
			gdouble scale = bb->x2 - bb->x1;
			GtsMatrix * sc;

			if (bb->y2 - bb->y1 > scale) scale = bb->y2 - bb->y1;
			if (bb->z2 - bb->z1 > scale) scale = bb->z2 - bb->z1;
			if (scale > 0.) scale = 1./scale;
			else scale = 1.;
			sc = gts_matrix_identity (NULL);
			sc[0][3] = - (bb->x1 + bb->x2)/2.;
			sc[1][3] = - (bb->y1 + bb->y2)/2.;
			sc[2][3] = - (bb->z1 + bb->z2)/2.;
			gts_surface_foreach_vertex (pSignedDistance->pSurface, (GtsFunc) gts_point_transform, sc);
			sc[0][0] = sc[1][1] = sc[2][2] = scale;    
			sc[0][3] = sc[1][3] = sc[2][3] = 0.;
			gts_surface_foreach_vertex (pSignedDistance->pSurface, (GtsFunc) gts_point_transform, sc);
			gts_matrix_destroy (sc);


			gts_surface_print_stats(pSignedDistance->pSurface, stdout);
			printf("#   Total area (normalized): %g\n", gts_surface_area(pSignedDistance->pSurface));
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

		// @TODO: Definir funcoes usadas como parametro;
		// @TODO: Melhorar o objeto discretizado tal que:
		//	aresta media menor que deltas da malha euleriana (faixa para tamanhos de arestas);
		// gts_surface_refine(pSignedDistance->pSurface, NULL, NULL, NULL, NULL, refine_stop, refine_stop_param);

		////////////////////////////////////////////////////////////////////////

		for(i = 0; i < 3; i++)
		{
			pSignedDistance->coordMin[i] = beginMesh[i];
			pSignedDistance->coordMax[i] = endMesh[i];
			pSignedDistance->size[i]     = sizeMesh[i];
			pSignedDistance->delta[i]    = fabs( (pSignedDistance->coordMax[i] - pSignedDistance->coordMin[i]) / (gdouble)pSignedDistance->size[i] );
		}

#if 0
		// Use the minimum valid delta;
		pSignedDistance->delta[0] = 1.1 * (qstats.edge_length.mean + 2.0*qstats.edge_length.stddev);
		pSignedDistance->delta[0] = cpt_round_value(pSignedDistance->delta[0], MIN_DELTA_MESH);
		pSignedDistance->delta[1] = (1.0 + 0.05) * pSignedDistance->delta[0];
		pSignedDistance->delta[2] = (1.0 - 0.05) * pSignedDistance->delta[0];
#endif

		pSignedDistance->sigma = 3.0 * cpt_min(pSignedDistance->delta[0], cpt_min(pSignedDistance->delta[1], pSignedDistance->delta[2]));

		gdouble  delta_max = cpt_max(pSignedDistance->delta[0], cpt_max(pSignedDistance->delta[1], pSignedDistance->delta[2]));
		pSignedDistance->distCut = cpt_get_dist_cut(size_sup, delta_max);
		pSignedDistance->distMax = cpt_get_dist_max(pSignedDistance->distCut, delta_max);

		if(verbose)
		{
			printf("#\tdistCut = { %+1.8f }; \n", pSignedDistance->distCut);
			printf("#\tdistMax = { %+1.8f }; \n", pSignedDistance->distMax);
			printf("#\tsigma   = { %+1.8f }; \n", pSignedDistance->sigma);

			fflush(stdout);
			fflush(stderr);
		}

		////////////////////////////////////////////////////////////////////////

#if 0		

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
			pSignedDistance->coordMin[i]-= 0.5*size_box;
			pSignedDistance->coordMax[i]+= 0.5*size_box;
			// ... but assert the maximun coord build a box with
			// a dimention multiple of the delta in this axis;
			pSignedDistance->size[i]     = (size_t)ceil(fabs(pSignedDistance->coordMax[i] - pSignedDistance->coordMin[i]) / pSignedDistance->delta[i]);
			pSignedDistance->coordMax[i] = pSignedDistance->coordMin[i] + pSignedDistance->size[i] * pSignedDistance->delta[i];
		}

#endif
		if(verbose)
		{
			printf("#\n");
			printf("# Expanding box...\n");
			printf("#\tcoordMin = { %+1.8f , %+1.8f , %+1.8f}; \n", pSignedDistance->coordMin[0], pSignedDistance->coordMin[1], pSignedDistance->coordMin[2]);
			printf("#\tcoordMax = { %+1.8f , %+1.8f , %+1.8f}; \n", pSignedDistance->coordMax[0], pSignedDistance->coordMax[1], pSignedDistance->coordMax[2]);
			printf("#\tsize     = {   %6u ,   %6u ,   %6u}; \n", pSignedDistance->size[0]    , pSignedDistance->size[1]    , pSignedDistance->size[2]    );
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
			printf("# Creating the eulerian mesh [size = %u]...", mesh_size);
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

		////////////////////////////////////////////////////////////////////////

#if 1
		if(verbose)
		{
			printf("#\n");
			printf("# Normalizing distances...");
		}
		for(i = 0; i < mesh_size; i++)
		{
			pSignedDistance->value[ i ] /= pSignedDistance->distMax;
		}
		if(verbose)
		{
			printf("done!\n");
		}
#endif

		////////////////////////////////////////////////////////////////////////

		fflush(stdout);
		fflush(stderr);

		////////////////////////////////////////////////////////////////////////

		theClock = clock();
		qt_pts = cpt_eulerian_mesh_write_silo(pSignedDistance, "out/eulerianmesh.silo");
		theClock = (clock() - theClock);
		if(verbose)
		{
			printf("#\n");
			printf("# Total points in the interface: [%u ; %+1.8f sec];\n", qt_pts, theClock / DBL_CLOCKS_PER_SEC);
		}

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
