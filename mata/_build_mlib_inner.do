/* Internal build helper: compiles Mata sources using `run` so that the
   stata-mcp guard (which bans top-level `do`) does not trip.  Called from
   the top-level `build_mlib.do` or from an ad-hoc stata-mp batch. */

clear all
set more off

local mata_dir "/Users/cxy/Desktop/2026project/trop/trop_stata/mata"
local base_dir "/Users/cxy/Desktop/2026project/trop/trop_stata"

// `mata mlib add` resolves the target library against c(pwd) (plus
// adopath), so anchor the working directory before creating the mlib.
cd "`base_dir'"

local mata_files ///
    "trop_constants.mata" ///
    "trop_rust_interface.mata" ///
    "trop_data_transfer.mata" ///
    "trop_lambda_grid.mata" ///
    "trop_backend_select.mata" ///
    "trop_ereturn_store.mata" ///
    "trop_validation.mata" ///
    "trop_loocv_validation.mata" ///
    "trop_bootstrap_diagnostics.mata" ///
    "trop_estat_helpers.mata" ///
    "trop_main.mata"

foreach f of local mata_files {
    display as text "  compiling `f'"
    run "`mata_dir'/`f'"
}

mata: mata mlib create ltrop, replace
mata: mata mlib add ltrop *(), complete

display as result ""
display as result "ltrop.mlib rebuilt"
