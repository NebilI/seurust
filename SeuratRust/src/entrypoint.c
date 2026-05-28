// Forward routine registration from C to Rust so the linker keeps the staticlib.
void R_init_SeuratRust_extendr(void *dll);
void register_extendr_panic_hook(void);

void R_init_SeuratRust(void *dll) {
  register_extendr_panic_hook();
  R_init_SeuratRust_extendr(dll);
}
