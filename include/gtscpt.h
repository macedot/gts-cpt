#ifndef __GTSCPT_H__
#define __GTSCPT_H__

// Standart C++ library
#include <iostream>
#include <cstdlib>
#include <unistd.h>
// Blitz related
//#include <blitz/numinquire.h>
//#include <blitz/array.h>

#ifdef __cplusplus
extern "C" {
#endif

// Standart C library
#include <stdlib.h>
#include <locale.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <memory.h>

// GTS related
#include "gts.h"
#include "gtstools.h"

#ifdef __cplusplus
}
#endif

//#define SIGN(x)                 ( (x) < 0.0 ? (-1.0) : (1.0) )
#define SIGN(x)                 ( (x) / fabs(x) )


typedef size_t  SizeVector[3];

// Defines one mesh; 
typedef struct eulerianMesh {
	gdouble      *value;         // raw mesh [z, y, x] : mesh_pos = k*(sizeZ * sizeY) + j*sizeY + i;
	gdouble      distCut;       // (epson_1) half size of band;
	gdouble      distMax;       // (epson_2) maximum distance between an eulerian point and the surface;
	gdouble      sigma;         // = 3 * MIN (HX , HY , HZ);
	GtsSurface*  pSurface;      // Lagrangian Mesh

// @TODO
//	struct {
	GtsVector   coordMin;      // minimum coord in the mesh
	GtsVector   coordMax;      // maximum coord in the mesh
	GtsVector   delta;         // step in the mesh;
	SizeVector  size;          // mesh size in delta unit;
//	} eulerianMesh;
	
} SignedDistance;


#endif // __GTSCPT_H__
