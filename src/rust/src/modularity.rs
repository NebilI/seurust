//! FFI wrapper around the proven C++ ModularityOptimizer implementation.
extern "C" {
    fn modularity_cluster_from_snn(
        snn_x: *const f64,
        snn_x_len: i32,
        snn_i: *const i32,
        snn_i_len: i32,
        snn_p: *const i32,
        snn_p_len: i32,
        snn_nrows: i32,
        snn_ncols: i32,
        modularity_function: i32,
        resolution: f64,
        algorithm: i32,
        n_random_starts: i32,
        n_iterations: i32,
        random_seed: i64,
        print_output: i32,
        edge_filename: *const std::ffi::c_char,
        out_len: *mut i32,
        error_msg: *mut std::ffi::c_char,
        error_msg_len: i32,
    ) -> *mut i32;

    fn modularity_free(ptr: *mut i32);
}

pub fn run_modularity_clustering_impl(
    x: &[f64],
    i: &[i32],
    p: &[i32],
    nrows: i32,
    ncols: i32,
    modularity_function: i32,
    resolution: f64,
    algorithm: i32,
    n_random_starts: i32,
    n_iterations: i32,
    random_seed: i32,
    print_output: bool,
    edge_filename: &str,
) -> Result<Vec<i32>, String> {
    let mut out_len = 0i32;
    let mut err_buf = vec![0u8; 512];
    let edge_c = std::ffi::CString::new(edge_filename).map_err(|e| e.to_string())?;

    let ptr = unsafe {
        modularity_cluster_from_snn(
            x.as_ptr(),
            x.len() as i32,
            i.as_ptr(),
            i.len() as i32,
            p.as_ptr(),
            p.len() as i32,
            nrows,
            ncols,
            modularity_function,
            resolution,
            algorithm,
            n_random_starts,
            n_iterations,
            random_seed as i64,
            i32::from(print_output),
            edge_c.as_ptr(),
            &mut out_len,
            err_buf.as_mut_ptr() as *mut std::ffi::c_char,
            err_buf.len() as i32,
        )
    };

    if ptr.is_null() {
        let msg = err_buf.split(|&b| b == 0).next().unwrap_or(b"unknown error");
        return Err(String::from_utf8_lossy(msg).into_owned());
    }

    let clusters = unsafe { std::slice::from_raw_parts(ptr, out_len as usize).to_vec() };
    unsafe { modularity_free(ptr) };
    Ok(clusters)
}
