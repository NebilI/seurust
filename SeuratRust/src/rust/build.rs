use std::path::Path;

fn main() {
    let cpp_dir = Path::new("../cpp");
    cc::Build::new()
        .cpp(true)
        .flag_if_supported("-std=c++17")
        .file(cpp_dir.join("ModularityOptimizer.cpp"))
        .file("cpp/modularity_bridge.cpp")
        .include(cpp_dir)
        .compile("seurat_modularity");

    if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=c++");
    } else {
        println!("cargo:rustc-link-lib=stdc++");
    }
}
