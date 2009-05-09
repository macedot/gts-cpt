#ifndef __GTSTOOLS_H__
#define __GTSTOOLS_H__

void build_list(gpointer data, GSList ** list);
void Curvature_Test(GtsSurface *surface);
void ScaleSurface(GtsVertex * v, double data[4]);
void MoveSurface(GtsVertex * v, double data[4]);
void Surface_FindBox(GtsSurface *surface, double coordMin[3], double coordMax[3]);

#endif //__GTSTOOLS_H__
