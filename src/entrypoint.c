// Unified DLL initializer for extendr (Rust) and Rcpp (C++) routines.
#include <R_ext/Rdynload.h>

void R_init_Seurat_extendr(void *dll);
void register_extendr_panic_hook(void);

extern void R_init_Seurat_rcpp(DllInfo *dll);

void R_init_Seurat(DllInfo *dll) {
  register_extendr_panic_hook();
  R_init_Seurat_extendr(dll);
  R_init_Seurat_rcpp(dll);
}
