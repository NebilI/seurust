use std::path::Path;

fn main() {
    let seurat_src = Path::new("..");
    cc::Build::new()
        .cpp(true)
        .flag_if_supported("-std=c++17")
        .file(seurat_src.join("ModularityOptimizer.cpp"))
        .file("cpp/modularity_bridge.cpp")
        .include(seurat_src)
        .compile("seurat_modularity");

    if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=c++");
    } else {
        println!("cargo:rustc-link-lib=stdc++");
    }
}
