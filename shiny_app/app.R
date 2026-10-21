library(shiny)

app_dir <- normalizePath(getwd(), mustWork = TRUE)
source(file.path(app_dir, "R", "job_utils.R"), local = TRUE)

options(shiny.maxRequestSize = 2 * 1024^2)

registry_path <- Sys.getenv(
  "GETPHENO_SHINY_CONFIG",
  unset = file.path(app_dir, "projects.conf")
)
project_registry <- read_project_registry(registry_path)
app_user <- current_unix_user()
user_registry <- accessible_registry(project_registry)

help_label <- function(label, help) {
  tagList(
    span(label),
    tags$span(
      class = "field-help",
      `data-tooltip` = help,
      tabindex = "0",
      role = "img",
      `aria-label` = help,
      icon("circle-info")
    )
  )
}

definition_input <- function(prefix, label, help, required = FALSE, placeholder = NULL) {
  modes <- if (required) c("Paste" = "paste", "Upload" = "upload") else
    c("None" = "none", "Paste" = "paste", "Upload" = "upload")
  default <- if (required) "paste" else "none"
  if (is.null(placeholder)) {
    placeholder <- if (required) "MDD\tICD10:F32,ICD10:F33" else "Name\tICD10:code1,ICD10:code2"
  }
  tagList(
    radioButtons(
      paste0(prefix, "_mode"), help_label(label, help), choices = modes, selected = default,
      inline = TRUE
    ),
    conditionalPanel(
      condition = sprintf("input.%s_mode == 'upload'", prefix),
      fileInput(
        paste0(prefix, "_file"), NULL,
        accept = c(".txt", ".tsv", ".csv", "text/plain", "text/tab-separated-values")
      )
    ),
    conditionalPanel(
      condition = sprintf("input.%s_mode == 'paste'", prefix),
      textAreaInput(
        paste0(prefix, "_text"), NULL, rows = if (required) 9 else 5,
        placeholder = placeholder
      )
    )
  )
}

project_choices <- function(registry) {
  projects <- !duplicated(registry$Project)
  stats::setNames(registry$Project[projects], registry$ProjectLabel[projects])
}

optional_section <- function(title, ...) {
  tags$details(
    class = "form-section optional-section",
    tags$summary(
      span(title),
      span(class = "optional-badge", "Optional")
    ),
    div(class = "optional-content", ...)
  )
}

ui <- navbarPage(
  title = "GetPheno",
  id = "main_navigation",
  header = tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "app.css")
  ),
  tabPanel(
    "New run",
    fluidPage(
      div(
        class = "identity-bar",
        span(icon("user"), "Unix user"),
        tags$strong(app_user),
        span(class = "identity-separator", "|"),
        span(sprintf("%d accessible protocol%s", nrow(user_registry), if (nrow(user_registry) == 1L) "" else "s"))
      ),
      conditionalPanel(
        condition = sprintf("%s", if (nrow(user_registry)) "false" else "true"),
        div(class = "alert alert-danger", "No readable project protocol with a writable output directory is configured for this Unix user.")
      ),
      fluidRow(
        column(
          width = 4,
          div(
            class = "form-section",
            h3("Run context"),
            selectInput(
              "project",
              help_label("Project", "Only projects readable by your current Unix account are listed."),
              choices = project_choices(user_registry)
            ),
            selectInput(
              "protocol",
              help_label("Protocol", "Selects the administrator-configured input files, INI settings, and Slurm resources."),
              choices = character()
            ),
            uiOutput("freeze_selection"),
            uiOutput("additional_selection"),
            uiOutput("context_access")
          ),
          div(
            class = "form-section",
            h3("Phenotypes"),
            definition_input(
              "phenotype", "Definition source",
              paste(
                "Paste or upload tab-separated definitions, one phenotype per line.",
                "Example: MDD [tab] ICD10:F32,ICD10:F33.",
                "Advanced example: MDD [tab] main=ICD10:F32,ICD10:F33;rule_out=ICD10:F25.",
                "Use sub= for supporting codes. A trailing * always requests prefix matching."
              ),
              required = TRUE
            ),
            uiOutput("phenotype_name"),
            textOutput("phenotype_summary")
          ),
          optional_section(
            "IID subset",
            definition_input(
              "iid_subset", "Subset source",
              "Optionally paste or upload one IID per line with no header. Only these IIDs can appear as cases or controls.",
              placeholder = "1000001\n1000002\n1000003"
            )
          )
        ),
        column(
          width = 4,
          optional_section(
            "Cohort and matching",
            checkboxInput(
              "ipsych", help_label(
                "iPSYCH comparison",
                "Builds the requested phenotypes from raw diagnosis data and applies the configured iPSYCH cohort, birth-date, and diagnosis-date restrictions."
              ), FALSE
            ),
            checkboxInput(
              "exact_match", help_label(
                "Exact matching (--eM)",
                "When selected, codes without * match exactly and codes ending in * remain prefix matches. Without it, all codes are prefix matches."
              ), FALSE
            ),
            radioButtons(
              "icd_mode", help_label(
                "Diagnostic classification",
                "Optionally restrict processing to ICD8 or ICD10. All configured keeps every supported classification."
              ),
              choices = c("All configured" = "all", "ICD8 only" = "ICD8", "ICD10 only" = "ICD10"),
              selected = "all"
            ),
            selectInput(
              "gender", help_label("Gender", "Optionally keep only individuals recorded as female or male."),
              choices = c("All" = "", "Female" = "F", "Male" = "M")
            ),
            textInput(
              "fyob", help_label("Earliest birth date", "Excludes individuals born before this date. Use YYYY-MM-DD."),
              placeholder = "YYYY-MM-DD"
            ),
            textInput(
              "min_max_age", help_label(
                "Age at first diagnosis",
                "Inclusive minimum and maximum age at first qualifying diagnosis, written min,max; 0,0 disables this filter."
              ), value = "0,0"
            ),
            uiOutput("death_filters")
          ),
          optional_section(
            "Diagnosis types",
            textInput(
              "diag_exclusions", help_label("Exclude", "Comma-separated diagnosis types to remove, for example H,M."),
              placeholder = "H,M"
            ),
            textInput(
              "diag_inclusions", help_label("Include", "Comma-separated diagnosis types to retain, for example A,B."),
              placeholder = "A,B"
            )
          )
        ),
        column(
          width = 4,
          optional_section(
            "Exclusions",
            definition_input(
              "lifetime", "Lifetime exclusion",
              "Removes a case if any listed code occurs at any time. Enter comma-separated codes, or Name [tab] code1,code2; for example BPD [tab] ICD10:F30*,ICD10:F31*."
            ),
            definition_input(
              "post", "Post-onset exclusion",
              "Excludes main diagnoses occurring after the first listed exclusion code. Enter comma-separated codes, or Name [tab] code1,code2."
            ),
            definition_input(
              "prior", "One-year-prior exclusion",
              "Excludes main-diagnosis entries occurring within one year before a listed exclusion code. Enter comma-separated codes, or Name [tab] code1,code2."
            )
          ),
          optional_section(
            "Outputs",
            checkboxGroupInput(
              "formats", help_label(
                "Additional formats",
                "Optionally write fastGWA and/or PLINK2 phenotype files in addition to the standard TSV output."
              ),
              choices = c("fastGWA" = "fastgwa", "PLINK2" = "plink2")
            )
          )
        )
      ),
      div(
        class = "submit-band",
        actionButton(
          "submit", "Submit Slurm job", icon = icon("play"), class = "btn-primary",
          title = "Validate the selections, create a unique run code, and submit GetPheno through Slurm."
        ),
        uiOutput("submitted_run")
      )
    )
  ),
  tabPanel(
    "Retrieve run",
    fluidPage(
      div(
        class = "retrieve-toolbar",
        textInput(
          "lookup_code", help_label("Run code", "Enter the unique code shown when the Slurm job was submitted."),
          placeholder = "GP-YYYYMMDDTHHMMSSZ-xxxxxxxxxxxx"
        ),
        actionButton("lookup", "Find run", icon = icon("search"), title = "Find this run in your accessible project outputs."),
        actionButton("refresh", "Refresh", icon = icon("refresh"), title = "Refresh the run state and available files.")
      ),
      uiOutput("run_details"),
      tableOutput("run_metadata"),
      fluidRow(
        column(8, selectInput(
          "download_file", help_label("Run file", "Select one generated result, log, or run metadata file."),
          choices = character()
        )),
        column(4, br(), downloadButton("download", "Download", icon = icon("download")))
      )
    )
  )
)

server <- function(input, output, session) {
  submitted <- reactiveVal(NULL)
  located_run <- reactiveVal(NULL)
  refresh_counter <- reactiveVal(0L)

  observeEvent(input$project, {
    rows <- user_registry[user_registry$Project == input$project, , drop = FALSE]
    choices <- stats::setNames(rows$Protocol, rows$ProtocolLabel)
    selected <- if (length(choices)) choices[[1]] else ""
    updateSelectInput(session, "protocol", choices = choices, selected = selected)
  }, ignoreInit = FALSE)

  selected_record <- reactive({
    req(input$project, input$protocol)
    rows <- user_registry[
      user_registry$Project == input$project & user_registry$Protocol == input$protocol,
      , drop = FALSE
    ]
    validate(need(nrow(rows) == 1L, "Select an accessible project and protocol"))
    rows[1, , drop = FALSE]
  })

  available_freezes <- reactive(discover_freezes(selected_record()))

  output$freeze_selection <- renderUI({
    record <- selected_record()
    freezes <- available_freezes()
    if (!nzchar(trim_value(record$FreezePattern))) return(NULL)
    if (!nrow(freezes)) {
      return(div(class = "context-warning", icon("triangle-exclamation"), "No readable cohort freeze is available"))
    }
    selectInput(
      "freeze_file", help_label(
        "Cohort freeze",
        "Selects the cohort file supplied as -i. Available freezes are sorted newest first, with the newest selected by default."
      ),
      choices = stats::setNames(freezes$path, freezes$label),
      selected = freezes$path[[1]]
    )
  })

  output$additional_selection <- renderUI({
    record <- selected_record()
    standard <- split_config_paths(record$StandardAdditionalFiles)
    optional <- parse_optional_additional_files(record$OptionalAdditionalFiles)
    optional <- optional[file.exists(optional$path) & file.access(optional$path, 4L) == 0L, , drop = FALSE]
    if (!length(standard) && !nrow(optional)) return(NULL)
    tagList(
      if (length(standard)) {
        div(class = "standard-files", tags$strong("Standard information"), basename(standard))
      },
      if (nrow(optional)) {
        checkboxGroupInput(
          "optional_additional", help_label(
            "Additional information (optional)",
            "Adds selected protocol-approved files through -j; GetPheno merges them by IID with the standard information files."
          ),
          choices = stats::setNames(optional$key, optional$label)
        )
      }
    )
  })

  output$death_filters <- renderUI({
    record <- selected_record()
    death_column <- trim_value(record$DeathDateColumn)
    if (!nzchar(death_column)) return(NULL)
    tagList(
      textInput(
        "death_age_range",
        help_label(
          "Age at death",
          "Inclusive completed age-at-death range written min,max. 0,0 disables this filter. By default, an active range keeps only deceased people in range."
        ),
        value = "0,0"
      ),
      checkboxInput(
        "keep_alive",
        help_label(
          "Also retain living people",
          "With an active age-at-death range, also retain IIDs whose death date is blank."
        ),
        FALSE
      ),
      checkboxInput(
        "only_dead",
        help_label(
          "Only deceased people",
          "Keep only IIDs with a recorded death. This can be combined with the age-at-death range, but not with Also retain living people."
        ),
        FALSE
      ),
      div(class = "standard-files", paste("Death-date column:", death_column))
    )
  })

  output$context_access <- renderUI({
    record <- selected_record()
    div(
      class = "access-ok",
      icon("lock"),
      span("Access checked with the current Unix account")
    )
  })

  active_definition_text <- reactive({
    if (identical(input$phenotype_mode, "upload")) {
      if (is.null(input$phenotype_file)) return("")
      tryCatch(read_uploaded_text(input$phenotype_file$datapath), error = function(error) "")
    } else {
      normalize_definition_text(input$phenotype_text)
    }
  })

  active_phenotype_count <- reactive(phenotype_count(active_definition_text()))

  output$phenotype_name <- renderUI({
    if (active_phenotype_count() == 1L) {
      textInput(
        "phenotype_name_value",
        help_label("Phenotype name", "Output column name used when the definition contains exactly one phenotype."),
        value = "MainPheno"
      )
    }
  })

  output$phenotype_summary <- renderText({
    count <- active_phenotype_count()
    if (!count) return("No phenotype definition supplied")
    sprintf("%d phenotype%s detected", count, if (count == 1L) "" else "s")
  })

  uploaded_path <- function(file_input) {
    if (is.null(file_input)) "" else file_input$datapath %||% ""
  }

  materialize_optional <- function(prefix, run_dir, suffix) {
    materialize_definition(
      input[[paste0(prefix, "_mode")]],
      uploaded_path(input[[paste0(prefix, "_file")]]),
      input[[paste0(prefix, "_text")]],
      file.path(run_dir, suffix),
      paste(prefix, "exclusion"),
      required = FALSE
    )
  }

  observeEvent(input$submit, {
    tryCatch({
      record <- selected_record()
      if (!record_is_accessible(record)) stop("Project access is no longer available")

      count <- active_phenotype_count()
      if (count < 1L) stop("Supply a phenotype definition")
      run_id <- new_run_id()
      user_root <- file.path(record$OutputRoot, app_user)
      if (!dir.exists(user_root) && !dir.create(user_root, recursive = TRUE, mode = "0700")) {
        stop("Could not create the user output directory")
      }
      Sys.chmod(user_root, mode = "0700")
      run_dir <- file.path(user_root, run_id)
      while (dir.exists(run_dir)) {
        run_id <- new_run_id()
        run_dir <- file.path(user_root, run_id)
      }
      if (!dir.create(run_dir, mode = "0700")) stop("Could not create the run directory")

      phenotype_file <- materialize_definition(
        input$phenotype_mode,
        uploaded_path(input$phenotype_file),
        input$phenotype_text,
        file.path(run_dir, paste0(run_id, ".phenotypes.tsv")),
        "Phenotype definition",
        required = TRUE
      )
      iid_subset_file <- materialize_definition(
        input$iid_subset_mode,
        uploaded_path(input$iid_subset_file),
        input$iid_subset_text,
        file.path(run_dir, paste0(run_id, ".iids.txt")),
        "IID subset",
        required = FALSE
      )
      exclusions <- list(
        lifetime = materialize_optional("lifetime", run_dir, paste0(run_id, ".lifetime.tsv")),
        post = materialize_optional("post", run_dir, paste0(run_id, ".post.tsv")),
        prior = materialize_optional("prior", run_dir, paste0(run_id, ".one-year-prior.tsv"))
      )
      protocol_files <- resolve_protocol_files(
        record,
        input$freeze_file %||% "",
        input$optional_additional %||% character()
      )
      output_file <- file.path(run_dir, paste0(run_id, ".tsv"))
      has_death_data <- nzchar(trim_value(protocol_files$death_date_column))
      command <- build_getpheno_command(
        record, run_id, phenotype_file, output_file,
        options = list(
          ipsych = isTRUE(input$ipsych),
          fastgwa = "fastgwa" %in% input$formats,
          plink2 = "plink2" %in% input$formats,
          exact_match = isTRUE(input$exact_match),
          gender = input$gender,
          fyob = input$fyob,
          min_max_age = input$min_max_age,
          death_age_range = if (has_death_data) input$death_age_range %||% "0,0" else "0,0",
          keep_alive = has_death_data && isTRUE(input$keep_alive),
          only_dead = has_death_data && isTRUE(input$only_dead),
          diag_exclusions = input$diag_exclusions,
          diag_inclusions = input$diag_inclusions,
          icd_mode = input$icd_mode,
          single_phenotype = count == 1L,
          name = input$phenotype_name_value %||% ""
        ),
        exclusions = exclusions,
        protocol_files = protocol_files,
        iid_subset = iid_subset_file %||% ""
      )
      script <- write_slurm_script(record, run_id, run_dir, command)
      metadata_path <- file.path(run_dir, "metadata.tsv")
      metadata <- list(
        run_id = run_id,
        user = app_user,
        project = record$Project,
        project_label = record$ProjectLabel,
        protocol = record$Protocol,
        protocol_label = record$ProtocolLabel,
        submitted_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
        job_id = "",
        output_file = basename(output_file),
        cohort_freeze = basename(protocol_files$stam %||% ""),
        additional_files = paste(basename(protocol_files$additional %||% character()), collapse = ","),
        command = command$display
      )
      write_run_metadata(metadata_path, metadata)
      write_private_text(file.path(run_dir, "status.txt"), "SUBMITTING")

      job_id <- tryCatch(
        submit_slurm_job(script, record$Sbatch %||% "sbatch"),
        error = function(error) {
          write_private_text(file.path(run_dir, "status.txt"), "SUBMISSION_FAILED")
          stop(error)
        }
      )
      metadata$job_id <- job_id
      write_run_metadata(metadata_path, metadata)
      write_private_text(file.path(run_dir, "status.txt"), "SUBMITTED")
      submitted(list(run_id = run_id, job_id = job_id, run_dir = run_dir))
      updateTextInput(session, "lookup_code", value = run_id)
      showNotification(paste("Submitted Slurm job", job_id), type = "message")
    }, error = function(error) {
      showNotification(conditionMessage(error), type = "error", duration = NULL)
    })
  })

  output$submitted_run <- renderUI({
    run <- submitted()
    if (is.null(run)) return(NULL)
    div(
      class = "run-code",
      span("Run code"),
      tags$strong(run$run_id),
      span(paste("Slurm job", run$job_id))
    )
  })

  locate_requested_run <- function() {
    code <- trim_value(input$lookup_code)
    if (!is_valid_run_id(code)) stop("Enter a valid GetPheno run code")
    run_dir <- find_run_directory(user_registry, app_user, code)
    if (is.null(run_dir)) stop("Run not found in your accessible project outputs")
    located_run(run_dir)
    refresh_counter(refresh_counter() + 1L)
    files <- list_run_files(run_dir)
    updateSelectInput(session, "download_file", choices = basename(files))
  }

  observeEvent(input$lookup, {
    tryCatch(locate_requested_run(), error = function(error) {
      located_run(NULL)
      showNotification(conditionMessage(error), type = "error")
    })
  })

  observeEvent(input$refresh, {
    if (is.null(located_run())) {
      tryCatch(locate_requested_run(), error = function(error) {
        showNotification(conditionMessage(error), type = "error")
      })
    } else {
      refresh_counter(refresh_counter() + 1L)
      files <- list_run_files(located_run())
      updateSelectInput(session, "download_file", choices = basename(files))
    }
  })

  output$run_details <- renderUI({
    refresh_counter()
    run_dir <- located_run()
    if (is.null(run_dir)) return(div(class = "empty-state", "Enter a run code to inspect its status and files."))
    metadata <- read_run_metadata(file.path(run_dir, "metadata.tsv"))
    div(
      class = "run-status",
      div(span("Run"), tags$strong(metadata$run_id %||% basename(run_dir))),
      div(span("State"), tags$strong(query_job_state(run_dir))),
      div(span("Slurm job"), tags$strong(metadata$job_id %||% "pending"))
    )
  })

  output$run_metadata <- renderTable({
    run_dir <- located_run()
    req(run_dir)
    metadata <- read_run_metadata(file.path(run_dir, "metadata.tsv"))
    visible <- intersect(
      c("project_label", "protocol_label", "submitted_utc", "output_file"),
      names(metadata)
    )
    data.frame(Field = gsub("_", " ", visible), Value = unlist(metadata[visible]), row.names = NULL)
  }, striped = TRUE, bordered = FALSE, spacing = "s")

  output$download <- downloadHandler(
    filename = function() basename(input$download_file %||% "run-file"),
    content = function(destination) {
      run_dir <- located_run()
      req(run_dir, input$download_file)
      allowed <- list_run_files(run_dir)
      source <- allowed[basename(allowed) == basename(input$download_file)]
      validate(need(length(source) == 1L, "Selected file is unavailable"))
      if (!file.copy(source, destination, overwrite = TRUE)) stop("Could not copy run file")
    },
    contentType = "application/octet-stream"
  )
}

shinyApp(ui, server)
