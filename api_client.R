###########################################################################################################
## Thin R client for the Quickomics data API (see the QuickomicsAPI project).
##
## Only invoked for projects loaded via ?serverfile=... / ?testfile=... URLs
## (i.e. when ProjectInfo$source is set to "published"/"unpublished" in
## inputdata.R). Every other data source -- the local data/ folder, the
## unlisted/ folder, and user file uploads -- is untouched and continues to
## use base::load() exactly as before.
###########################################################################################################

api_base_url <- Sys.getenv("QUICKOMICS_API_URL",      "http://localhost:8499")
api_username <- Sys.getenv("QUICKOMICS_API_USERNAME", "admin")
api_password <- Sys.getenv("QUICKOMICS_API_PASSWORD", "")

# Directory containing the API's docker-compose.yml -- only used by
# ensure_api_running() below. Auto-start is opt-in: set
# QUICKOMICS_API_AUTOSTART=true to enable it. Off by default, since having
# application code spawn `docker` commands is a meaningful trade-off (see
# comment on ensure_api_running() below) that shouldn't happen silently.
api_compose_dir  <- Sys.getenv("QUICKOMICS_API_COMPOSE_DIR", "")
api_autostart_on <- tolower(Sys.getenv("QUICKOMICS_API_AUTOSTART", "false")) %in% c("true", "1", "yes")

# Simple in-process token cache -- shared across all Shiny sessions in this
# R process, since the API token represents the app's own service identity,
# not an individual user. Survives for the life of the R process; a fresh
# Shiny app start re-fetches a token on first use.
.api_token_cache <- new.env(parent = emptyenv())
.api_token_cache$token <- NULL

# Also track whether we've already tried to auto-start the API this R
# session, so a genuinely broken API doesn't get a `docker compose up`
# attempt on every single request.
.api_autostart_cache <- new.env(parent = emptyenv())
.api_autostart_cache$attempted <- FALSE

#' Check whether the API is currently reachable.
api_health_check <- function(timeout_sec = 2) {
  resp <- tryCatch(
    httr::GET(paste0(api_base_url, "/health"), httr::timeout(timeout_sec)),
    error = function(e) NULL
  )
  !is.null(resp) && httr::status_code(resp) == 200
}

#' If the API isn't reachable, try to start it via `docker compose up -d`
#' and poll /health until it responds (or we give up after max_wait_sec).
#'
#' This is a deliberate trade-off, not a default-on convenience: it means
#' the QuickOmics R process needs permission to run `docker` commands (the
#' OS user it runs as must be in the `docker` group, or have passwordless
#' sudo configured for it) -- effectively giving application code control
#' over infrastructure lifecycle, which is normally kept separate. The more
#' conventional alternative is simply leaving the API container always-on
#' (docker-compose.yml already has `restart: unless-stopped`), managed
#' independently of QuickOmics -- e.g. via systemd, or just... leaving it
#' running. This function exists for dev/testing convenience (e.g. running
#' QuickOmics from RStudio without remembering to start the API container
#' first), and is OFF by default -- set QUICKOMICS_API_AUTOSTART=true and
#' QUICKOMICS_API_COMPOSE_DIR=/path/to/QuickomicsAPI to enable it.
ensure_api_running <- function(max_wait_sec = 30, poll_interval_sec = 2) {
  if (api_health_check()) return(invisible(TRUE))

  if (!api_autostart_on) {
    stop(
      "Quickomics API is not reachable at ", api_base_url, " and auto-start ",
      "is disabled. Start it manually (`docker compose up -d` in the API ",
      "project directory), or set QUICKOMICS_API_AUTOSTART=true and ",
      "QUICKOMICS_API_COMPOSE_DIR to enable automatic startup."
    )
  }

  if (.api_autostart_cache$attempted) {
    stop(
      "Quickomics API is still not reachable at ", api_base_url,
      " after a previous auto-start attempt this session. Not retrying ",
      "automatically -- check `docker compose logs` in ", api_compose_dir,
      " for why the container isn't coming up."
    )
  }
  .api_autostart_cache$attempted <- TRUE

  if (!nzchar(api_compose_dir)) {
    stop(
      "QUICKOMICS_API_AUTOSTART is enabled but QUICKOMICS_API_COMPOSE_DIR ",
      "is not set -- don't know where to find docker-compose.yml."
    )
  }
  compose_file <- file.path(api_compose_dir, "docker-compose.yml")
  if (!file.exists(compose_file)) {
    stop("QUICKOMICS_API_COMPOSE_DIR is set to '", api_compose_dir,
         "' but no docker-compose.yml was found there.")
  }

  message("Quickomics API not reachable at ", api_base_url,
          " -- attempting `docker compose up -d` in ", api_compose_dir, " ...")

  out <- tryCatch(
    system2("docker", c("compose", "-f", shQuote(compose_file), "up", "-d"),
            stdout = TRUE, stderr = TRUE),
    error = function(e) {
      stop("Failed to run `docker compose up -d`: ", conditionMessage(e),
           ". Does this R process have permission to run docker commands?")
    }
  )
  message(paste(out, collapse = "\n"))

  waited <- 0
  while (waited < max_wait_sec) {
    Sys.sleep(poll_interval_sec)
    waited <- waited + poll_interval_sec
    if (api_health_check()) {
      message("Quickomics API is up (took ~", waited, "s).")
      return(invisible(TRUE))
    }
  }

  stop(
    "Started `docker compose up -d` but the API still isn't responding at ",
    api_base_url, " after ", max_wait_sec, "s. Check `docker compose logs` ",
    "in ", api_compose_dir, " for errors."
  )
}

# Fetch a fresh token from the API and cache it.
api_fetch_token <- function() {
  resp <- httr::POST(
    paste0(api_base_url, "/auth/token"),
    httr::content_type_json(),
    body = list(username = api_username, password = api_password),
    encode = "json"
  )
  if (httr::status_code(resp) != 200) {
    stop("Quickomics API: failed to obtain auth token (HTTP ", httr::status_code(resp), ")")
  }
  parsed <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
  # /auth/token currently box-wraps scalars as one-element arrays, e.g.
  # {"access_token":["..."]}, so unwrap defensively either way.
  token <- if (is.list(parsed$access_token) || length(parsed$access_token) > 1) {
    parsed$access_token[[1]]
  } else {
    parsed$access_token
  }
  .api_token_cache$token <- token
  token
}

# Return the cached token, fetching one if this is the first call since the
# R process started.
api_get_token <- function() {
  if (is.null(.api_token_cache$token)) api_fetch_token() else .api_token_cache$token
}

#' Load a project's full data bundle from the Quickomics API, assigning the
#' same object names base::load(RDataFile) would have created directly into
#' `envir` -- so existing code written against `objs <- load(RDataFile)` and
#' bare object names (MetaData, ProteinGeneName, data_wide, ...) keeps
#' working completely unchanged after switching the call site to
#' `objs <- api_load(project_id, source)`.
#'
#' @param project_id Project ID (matches ProjectInfo$ProjectID)
#' @param source "published" or "unpublished" (matches ProjectInfo$source)
#' @param envir Environment to assign loaded objects into; defaults to the
#'   caller's own environment, exactly like base::load()'s default behavior.
#' @return Character vector of object names assigned -- mirrors load()'s own
#'   return value, so `identical(objs, c("results_long", "ProteinGeneName"))`
#'   style checks downstream keep working unmodified.
api_load <- function(project_id, source, envir = parent.frame()) {
  ensure_api_running()

  fetch_once <- function(token) {
    httr::GET(
      paste0(api_base_url, "/datasets/", source, "/", project_id, "/bundle.rds"),
      httr::add_headers(Authorization = paste("Bearer", token))
    )
  }

  token <- api_get_token()
  resp <- fetch_once(token)

  # Cached token may have expired since it was minted -- refresh once and retry.
  if (httr::status_code(resp) == 401) {
    token <- api_fetch_token()
    resp <- fetch_once(token)
  }

  if (httr::status_code(resp) != 200) {
    stop(
      "Quickomics API: failed to load project '", project_id,
      "' (source='", source, "'), HTTP ", httr::status_code(resp),
      ". Response: ", httr::content(resp, "text", encoding = "UTF-8")
    )
  }

  bundle <- unserialize(httr::content(resp, "raw"))

  assigned <- character(0)
  assign_if_present <- function(var_name, bundle_key) {
    value <- bundle[[bundle_key]]
    if (!is.null(value)) {
      # data.table objects carry an internal self-reference pointer
      # (.internal.selfref) that isn't guaranteed to survive a
      # serialize()/unserialize() round trip cleanly -- this is a known
      # data.table quirk, not specific to our API, but this is the exact
      # path that round trip happens on. setDT() is data.table's own
      # official remedy; it's a no-op for anything that isn't actually a
      # data.table. This matters most for data_results specifically:
      # unlike results_long/data_long, DataReactive() never explicitly
      # rebuilds it downstream, so it's the one object with no other
      # safety net if this goes uncorrected.
      if (data.table::is.data.table(value)) {
        data.table::setDT(value)
      }
      assign(var_name, value, envir = envir)
      assigned <<- c(assigned, var_name)
    }
  }
  # Left side: object name base::load() would create.
  # Right side: corresponding key in the API's /bundle(.rds) response.
  assign_if_present("MetaData",        "samples")
  assign_if_present("ProteinGeneName", "genes")
  assign_if_present("data_wide",       "data_wide")
  assign_if_present("data_long",       "expression_long")
  assign_if_present("results_long",    "results")
  assign_if_present("data_results",    "data_results")
  assign_if_present("comp_info",       "comp_info")

  # inputdata.R checks `identical(objs, c("results_long", "ProteinGeneName"))`
  # to detect comparison-only projects -- identical() on a character vector
  # is order-sensitive, but assign_if_present() above always adds
  # ProteinGeneName before results_long (alphabetical-ish call order), which
  # would never match that literal ordering. The API already computes this
  # flag correctly and order-independently server-side (bundle$comp_only),
  # so use it directly here rather than trying to reconstruct load()'s
  # original object order, which isn't something the API can know anyway --
  # it depends on the order objects happened to be saved in the .RData file.
  if (isTRUE(bundle$comp_only)) {
    assigned <- c("results_long", "ProteinGeneName")
  }

  assigned
}

#' Fetch the complete, unfiltered co-expression network for a project from
#' the API. Unlike api_load(), this returns the network data.frame directly
#' rather than assigning into an environment -- network.R's own code
#' assigns the result to a local `network` variable itself, so there's no
#' load()-style multi-object semantics to replicate here.
#'
#' @param project_id Project ID
#' @param source "published" or "unpublished"
#' @return The network data.frame (columns like from/to/cor/p, depending on
#'   how the source .RData file was structured).
api_load_network <- function(project_id, source) {
  ensure_api_running()

  fetch_once <- function(token) {
    httr::GET(
      paste0(api_base_url, "/datasets/", source, "/", project_id, "/network/bundle.rds"),
      httr::add_headers(Authorization = paste("Bearer", token))
    )
  }

  token <- api_get_token()
  resp <- fetch_once(token)

  if (httr::status_code(resp) == 401) {
    token <- api_fetch_token()
    resp <- fetch_once(token)
  }

  # 404 means "no precomputed network for this project" -- the same
  # meaning as a missing local file in the non-API path, i.e. not an
  # error, just a signal to compute fresh instead. Any other non-200
  # (auth failure, 5xx, unreachable API) is a genuine problem and should
  # surface loudly rather than be silently treated the same way.
  if (httr::status_code(resp) == 404) {
    return(NULL)
  }

  if (httr::status_code(resp) != 200) {
    stop(
      "Quickomics API: failed to load network for project '", project_id,
      "' (source='", source, "'), HTTP ", httr::status_code(resp),
      ". Response: ", httr::content(resp, "text", encoding = "UTF-8")
    )
  }

  network <- unserialize(httr::content(resp, "raw"))
  # Same data.table serialization-safety fix as api_load() -- network.R
  # never explicitly rebuilds this object downstream (unlike results_long/
  # data_long in DataReactive()), so it has the same exposure to a stale
  # .internal.selfref if the source file's network object is a data.table.
  if (data.table::is.data.table(network)) {
    data.table::setDT(network)
  }
  network
}

#' Load a project's WGCNA bundle (netwk, dataExpr, picked_power, sft) from
#' the API, assigning the same object names that
#' `load(wgcnafile); load(load_wgcna_file)` would have created directly
#' into `envir` -- mirrors all three call sites in wgcna.R that use that
#' exact two-load pattern today.
#'
#' Unlike api_load(), there's no comp_only-style ordering concern here --
#' wgcna.R never checks `identical()` on an object-name vector for this
#' data, it just checks individual objects' existence (`exists("sft")`
#' etc.), so plain presence/absence is all that matters.
#'
#' @param project_id Project ID
#' @param source "published" or "unpublished"
#' @param envir Environment to assign into; defaults to the caller's own
#'   environment, matching base::load()'s default behavior.
#' @return Invisible NULL. dataExpr/picked_power/sft may not be assigned at
#'   all if the API's load_{id}.RData-equivalent doesn't exist for this
#'   project (mirrors `exists("dataExpr")` being FALSE in the local case).
api_load_wgcna <- function(project_id, source, envir = parent.frame()) {
  ensure_api_running()

  fetch_once <- function(token) {
    httr::GET(
      paste0(api_base_url, "/datasets/", source, "/", project_id, "/wgcna/bundle.rds"),
      httr::add_headers(Authorization = paste("Bearer", token))
    )
  }

  token <- api_get_token()
  resp <- fetch_once(token)

  if (httr::status_code(resp) == 401) {
    token <- api_fetch_token()
    resp <- fetch_once(token)
  }

  if (httr::status_code(resp) != 200) {
    stop(
      "Quickomics API: failed to load WGCNA results for project '", project_id,
      "' (source='", source, "'), HTTP ", httr::status_code(resp),
      ". Response: ", httr::content(resp, "text", encoding = "UTF-8")
    )
  }

  bundle <- unserialize(httr::content(resp, "raw"))

  assign_if_present <- function(var_name) {
    value <- bundle[[var_name]]
    if (!is.null(value)) {
      if (data.table::is.data.table(value)) data.table::setDT(value)
      assign(var_name, value, envir = envir)
    }
  }
  assign_if_present("netwk")
  assign_if_present("dataExpr")
  assign_if_present("picked_power")
  assign_if_present("sft")

  invisible(NULL)
}
