systemc <- function(command, env = CONDA_ENV, ignore_stdout = FALSE, ignore_stderr = FALSE, conda_sh_path = CONDA_PATH) {
  com <- glue::glue("bash -c '. {conda_sh_path}; conda activate {env}; {command}'")
  system(com, ignore.stdout = ignore_stdout, ignore.stderr = ignore_stderr)
}


cmh_palette <- c('#50514f', '#559389', '#f25f5c', '#69af6d', '#247ba0', '#e2dd3f', '#7f6699', '#e08831')