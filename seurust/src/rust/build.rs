use std::path::PathBuf;
use std::process::Command;

fn rcpp_eigen_include() -> Option<PathBuf> {
    if let Ok(path) = std::env::var("RCPP_EIGEN_INCLUDE") {
        let p = PathBuf::from(path);
        if p.join("Eigen").exists() {
            return Some(p);
        }
    }

    let r_home = std::env::var("R_HOME").ok().or_else(|| {
        Command::new("R").arg("RHOME").output().ok().and_then(|o| {
            if o.status.success() {
                Some(String::from_utf8_lossy(&o.stdout).trim().to_string())
            } else {
                None
            }
        })
    })?;

    let candidates = [
        PathBuf::from(&r_home).join("library/RcppEigen/include"),
        PathBuf::from("/usr/lib/R/site-library/RcppEigen/include"),
        PathBuf::from("/usr/local/lib/R/site-library/RcppEigen/include"),
    ];

    for path in candidates {
        if path.join("Eigen").exists() {
            return Some(path);
        }
    }

    if let Ok(output) = Command::new(format!("{}/bin/Rscript", r_home))
        .args(["-e", "cat(system.file('include', package='RcppEigen'))"])
        .output()
    {
        if output.status.success() {
            let path = PathBuf::from(String::from_utf8_lossy(&output.stdout).trim());
            if path.join("Eigen").exists() {
                return Some(path);
            }
        }
    }

    None
}

fn r_home() -> Option<String> {
    std::env::var("R_HOME").ok().or_else(|| {
        Command::new("R").arg("RHOME").output().ok().and_then(|o| {
            if o.status.success() {
                Some(String::from_utf8_lossy(&o.stdout).trim().to_string())
            } else {
                None
            }
        })
    })
}

fn r_include() -> Option<PathBuf> {
    if let Some(r_home) = r_home() {
        let path = PathBuf::from(&r_home).join("include");
        if path.join("R.h").exists() {
            return Some(path);
        }
    }

    for path in [
        PathBuf::from("/usr/share/R/include"),
        PathBuf::from("/usr/local/lib/R/include"),
    ] {
        if path.join("R.h").exists() {
            return Some(path);
        }
    }

    None
}

fn rcpp_include() -> Option<PathBuf> {
    let r_home = r_home()?;
    let path = PathBuf::from(&r_home).join("library/Rcpp/include");
    if path.join("Rcpp").exists() {
        return Some(path);
    }

    for path in [
        PathBuf::from("/usr/lib/R/site-library/Rcpp/include"),
        PathBuf::from("/usr/local/lib/R/site-library/Rcpp/include"),
    ] {
        if path.join("Rcpp").exists() {
            return Some(path);
        }
    }

    if let Ok(output) = Command::new(format!("{}/bin/Rscript", r_home))
        .args(["-e", "cat(system.file('include', package='Rcpp'))"])
        .output()
    {
        if output.status.success() {
            let path = PathBuf::from(String::from_utf8_lossy(&output.stdout).trim());
            if path.join("Rcpp").exists() {
                return Some(path);
            }
        }
    }

    None
}

fn main() {
    println!("cargo:rustc-check-cfg=cfg(snn_eigen)");

    if let Some(eigen_inc) = rcpp_eigen_include() {
        println!("cargo:rerun-if-env-changed=RCPP_EIGEN_INCLUDE");
        println!("cargo:rerun-if-env-changed=R_HOME");
        println!("cargo:rustc-cfg=snn_eigen");
        let mut build = cc::Build::new();
        build
            .cpp(true)
            .flag_if_supported("-std=c++17")
            .flag_if_supported("-O3")
            .flag_if_supported("-Wno-ignored-attributes")
            .flag_if_supported("-Wno-cast-function-type")
            .define("SNN_BRIDGE_RCPP", None)
            .file("cpp/snn_bridge.cpp")
            .include(&eigen_inc);
        if let Some(r_inc) = r_include() {
            build.include(r_inc);
        }
        if let Some(rcpp_inc) = rcpp_include() {
            build.include(rcpp_inc);
        }
        build.compile("snn_eigen_bridge");
    } else {
        // Pure Rust fallback when RcppEigen headers are unavailable.
    }

    if cfg!(target_os = "macos") {
        println!("cargo:rustc-link-lib=c++");
    } else {
        println!("cargo:rustc-link-lib=stdc++");
    }
}
