package_list <- c(
    "ggplot2",
    "grid",
    "gridBase",
    "gridExtra",
    "stringr",
    "reshape2"
)

success <- TRUE
for (pkg in package_list) {
   #suppressWarnings(suppressPackageStartupMessages({
    if (!require(pkg, character.only=T)) {
        print(paste0("Package '",pkg,"' not installed. Attempting install..."));
        tryCatch(install.packages(pkg), error=function(e) print(paste0("Error installing package ",pkg)))
        if (!require(pkg, character.only=T)) {
            print(paste0("Package '",pkg,"' install failed..."));
            success <- FALSE
        } else {
            print(paste0("Package '",pkg,"' successfully installed..."));
        }
    } else {
        print(paste0("Package '",pkg,"' loaded."));
    }
  #}))
}
if (success) {
    print("All R packages loaded successfully.")
} else {
    warning(paste0("At least one required R package could not be loaded or installed (see error messages)\n",
            "These packages are required to generate KaryoScan diagnostic plots.\n",
            "Fix these install errors or plots will not be generated.\n",
            "Required packages: ",paste(package_list,sep="",collapse=", ")))
}
