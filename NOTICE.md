# Third-party notices

## ANN (Approximate Nearest Neighbor library)

`engine/src/ann_kdtree.cpp` and `engine/include/isvg/ann_kdtree.h` are a port of parts of ANN 1.1.1
(kd_tree, kd_split, kd_util, kd_search, pr_queue_k).

Copyright (c) 1997-2005 University of Maryland and Sunil Arya and David Mount.
Distributed under the GNU Lesser General Public License, version 2.1:
https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html
ANN: https://www.cs.umd.edu/~mount/ANN/

## stb_image

`engine/third_party/stb/stb_image.h` - Sean Barrett, public domain (or MIT, see the file's end).
https://github.com/nothings/stb

## Reference BLAS / LAPACK routines

`engine/include/isvg/blas_ref.h` reimplements a few reference BLAS/LAPACK routines (dgemm, dgemv,
dgesv and helpers) following the public netlib algorithms: https://www.netlib.org/lapack/
