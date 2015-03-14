The signed distance function between an arbitrary point in 3D space and a given closed surface returns the minimum distance from that point to the collection of triangles representing the surface. By convention, the sign is positive if the point is outside and negative if the point is inside the region determined by the surface.

In the context of Computational Fluid Dynamics, not only the signed distance function is useful to locate the separation interface between two different fluids with constant (but possibly different) material properties as its zero level-set surface but also it determines for each Eulerian grid point its instantaneous mass density and viscosity (one fluid is represented by the positive values and the other one by the negative ones).

Although there are other ways to compute the signed distance function, when the number of Lagrangian points (vertices of the triangles) is less than N\*log(N), where N is the number of points in the Eulerian grid, a very robust, most efficient way of computing the signed distance function is obtained through the application of techniques from Computational Geometry (Mauch´s Algorithm). The computational cost for the **entire Eulerian domain** is linear in the number of the geometrical elements composing the triangulated surface (that is, vertices, edges, and faces) with a O(1) proportion constant.

The current GTS-CPT code creates an Eulerian grid in run-time based on the average size of the triangle edges and computes the (cutoff) signed distance function. For details on its 2D computation and use for a two-phase flow simulations refer to Ceniceros and Roma 2005, and to Ceniceros et al. 2009.

The current stage is being developed by the undergraduate student Thiago P. Macedo, under the supervision of Prof. Alexandre M. Roma (roma (at) ime usp br) in the Instituto de Matematica e Estatistica at the Universidade de Sao Paulo (IME-USP).


REFERENCES: (If you used this implementation and want to thank us, please site us :- )

[1](1.md) Ceniceros and Roma 2005: (preprint) http://www.math.ucsb.edu/~hdc/public/left.pdf
> Journal of Computational Physics, Volume 205, [Issue 2](https://code.google.com/p/gts-cpt/issues/detail?id=2), 20 May 2005, Pages 391-400.
> doi:  10.1016/j.jcp.2004.11.013


[2](2.md) Ceniceros et al. 2009: http://www.math.ucsb.edu/~hdc/public/Adaptive_LEFT.pdf
> Commun. Comput. Phys., 8 (2010), pp. 51-94.
> doi:  10.4208/cicp.050509.141009a


[3](3.md) Mauch´s Algorithm: http://www.cacr.caltech.edu/~sean/projects/cpt/html3/index.html


---


The current version can be tested by downloading some GTS sample files examples at http://gts.sourceforge.net/samples.html. As of now, the current program version runs only with closed manifolds. One can remove this constraint in the source code and test it with other cases at his own risk. No warranties on what´s gonna happen are given =]

The computed **cut** distance function (indicator function) is output and may be visualized with VisIt (https://wci.llnl.gov/codes/visit/). To compile the example program, you'll need the SILO library installed (see downloads, or grab it at https://wci.llnl.gov/codes/visit/3rd_party/silo060605.sh). We added also at download section some dependences to run VisIt (1.11.1, 1.11.2, i386) with Silo support.


---


Please let us know if you decide to download and to use our implementation for this is important to us for future reference. We appreciate any comments about the above topics, the life, the universe and everything =]

Please write about problems/bugs or improvement suggestions to Thiago P. Macedo (tmacedo at usp br).


---


Some results:

https://www.dropbox.com/s/ag0pp12fcdlj3f8/horse01.png?dl=1

https://www.dropbox.com/s/jri4y9he93c01yc/horse05.png?dl=1

https://www.dropbox.com/s/eoac3en8fb4oodo/horse03.png?dl=1

https://www.dropbox.com/s/9raaflj7pogwvfx/viewer%202013-01-24%2017-46-35-09.png?dl=1

https://www.dropbox.com/s/jm099z5l70wx8cm/viewer%202013-01-24%2017-55-02-47.png?dl=1

https://www.dropbox.com/s/6nk3v1whg2g7fsp/viewer%202013-01-24%2018-22-26-28.png?dl=1