args <- commandArgs(trailingOnly = TRUE)

argument <- function(name, default = "") {
  prefix <- paste0("--", name, "=")
  match <- args[startsWith(args, prefix)]
  if (!length(match)) return(default)
  sub(prefix, "", match[[1]], fixed = TRUE)
}

script_arg <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", script_arg[startsWith(script_arg, "--file=")][[1]])
app_dir <- normalizePath(dirname(script_path), mustWork = TRUE)

config <- argument(
  "config",
  Sys.getenv("GETPHENO_SHINY_CONFIG", unset = file.path(app_dir, "projects.conf"))
)
host <- argument("host", "127.0.0.1")
port <- suppressWarnings(as.integer(argument("port", "3838")))
if (is.na(port) || port < 1L || port > 65535L) stop("Invalid --port value")
if (!file.exists(config)) stop("Project registry not found: ", config)
if (!requireNamespace("shiny", quietly = TRUE)) {
  stop("The R package 'shiny' is required. Run install.packages('shiny').")
}

Sys.setenv(GETPHENO_SHINY_CONFIG = normalizePath(config, mustWork = TRUE))
shiny::runApp(app_dir, host = host, port = port, launch.browser = FALSE)
