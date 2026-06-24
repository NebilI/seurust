// Forward routine registration from C to Rust so the linker keeps the staticlib.
void R_init_seurust_extendr(void *dll);
void register_extendr_panic_hook(void);

void R_init_seurust(void *dll) {
  register_extendr_panic_hook();
  R_init_seurust_extendr(dll);
}
