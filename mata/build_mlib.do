/* ==============================================================================
   build_mlib.do --- Rebuild trop_stata/ltrop.mlib from the current .mata sources.

   Produces trop_stata/ltrop.mlib containing every function defined across the
   dependency-ordered Mata source files.  Intended to be run from the project
   root (or any directory where trop_stata/ is reachable).
   ============================================================================== */

clear all
set more off

local mata_dir "trop_stata/mata"
local out_dir  "trop_stata"

capture confirm file "`mata_dir'/trop_constants.mata"
if _rc {
    local mata_dir "mata"
    local out_dir  "."
    capture confirm file "`mata_dir'/trop_constants.mata"
    if _rc {
        // Fall back to the absolute workspace path (stata-mp -b run from
        // an unrelated CWD) so this helper works regardless of invocation
        // directory.
        local mata_dir "/Users/cxy/Desktop/2026project/trop/trop_stata/mata"
        local out_dir  "/Users/cxy/Desktop/2026project/trop/trop_stata"
        capture confirm file "`mata_dir'/trop_constants.mata"
        if _rc {
            di as error "Cannot locate trop_stata Mata sources."
            exit 601
        }
    }
}

// Compile sources in dependency order.  This matches load_mata_once.do;
// any reordering must be mirrored there.
local mata_files ///
    trop_constants ///
    trop_rust_interface ///
    trop_data_transfer ///
    trop_lambda_grid ///
    trop_backend_select ///
    trop_ereturn_store ///
    trop_validation ///
    trop_loocv_validation ///
    trop_bootstrap_diagnostics ///
    trop_estat_helpers ///
    trop_main

foreach f of local mata_files {
    qui run "`mata_dir'/`f'.mata"
}

// Bundle all compiled functions into ltrop.mlib under trop_stata/.
mata: mata mlib create ltrop, dir("`out_dir'") replace
mata: mata mlib add ltrop *(), dir("`out_dir'")
mata: mata mlib index

di as result "Wrote `out_dir'/ltrop.mlib"
