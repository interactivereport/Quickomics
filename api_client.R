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
# QUICKOMICS_API_AUTOSTART=true to enable it. Off by default -- and should
# stay off entirely whenever QuickOmics and the API run on different
# servers, since `docker compose up -d` only ever controls the LOCAL
# machine's Docker daemon; it has no way to reach a daemon on another host.
api_compose_dir  <- Sys.getenv("QUICKOMICS_API_COMPOSE_DIR", "")
api_autostart_on <- tolower(Sys.getenv("QUICKOMICS_API_AUTOSTART", "false")) %in% c("true", "1", "yes")

# Retry/timeout settings for data-fetching calls (api_get_project_meta,
# api_load, api_load_network, api_load_wgcna) -- NOT api_health_check or
# api_fetch_token, see api_http_get_with_retry()'s own comment for why.
#
# Plumber runs as a single R process handling one request at a time (see
# https://www.rplumber.io/articles/execution-model.html) -- with only one
# pod (minReplicas: 1) and no way to change that from this side, a request
# that arrives while the server is busy with someone else's slow call
# simply queues and waits, rather than failing outright. These settings
# make the client tolerate that wait instead of giving up too early, and
# retry a modest number of times if the wait is long enough that an
# upstream proxy/load balancer times the connection out first (a 502/503/
# 504 in that case doesn't mean the request was bad, just that it didn't
# get a chance to finish before something in the middle gave up on it).
api_max_retries        <- as.integer(Sys.getenv("QUICKOMICS_API_MAX_RETRIES", "3"))
api_timeout_sec         <- as.integer(Sys.getenv("QUICKOMICS_API_TIMEOUT_SEC", "120"))
api_retry_backoff_sec   <- as.integer(Sys.getenv("QUICKOMICS_API_RETRY_BACKOFF_SEC", "2"))

# Simple in-process token cache -- shared across all Shiny sessions in this
# R process, since the API token represents the app's own service identity,
# not an individual user. Survives for the life of the R process; a fresh
# Shiny app start re-fetches a token on first use.
.api_token_cache <- new.env(parent = emptyenv())
.api_token_cache$token <- NULL

# Also track whether we've already tried to auto-start the API this R
# session, so a genuinely broken/unreachable API doesn't get a
# `docker compose up` attempt on every single request.
.api_autostart_cache <- new.env(parent = emptyenv())
.api_autostart_cache$attempted <- FALSE

#' Perform an HTTP GET with a generous timeout and retry-with-backoff for
#' transient failures -- specifically network-level errors (connection
#' dropped, DNS hiccup, etc.) and 502/503/504 responses, which are exactly
#' what you'd see from a proxy/load balancer giving up on a request that
#' got stuck queued behind another slow one on a single-threaded plumber
#' process, rather than the API itself rejecting the request.
#'
#' Deliberately does NOT retry on 400/401/403/404 -- those are the server
#' telling you something specific and non-transient (bad parameter, bad
#' token, insufficient role, not found); retrying the identical request
#' won't change any of those outcomes. 401 already gets its own one-shot
#' refresh-token-and-retry-once handling at each call site below, which is
#' a different, complementary mechanism to this one.
#'
#' Not used by api_health_check() (a fast, cheap liveness probe -- retrying
#' it defeats its purpose of quickly telling ensure_api_running() whether
#' the API is up) or api_fetch_token() directly (POST with a different
#' body shape; add retry there too if login-time flakiness becomes an
#' issue in practice).
#'
#' @param url Full request URL
#' @param headers httr::add_headers(...) result (or similar) to pass through
#' @return The httr response object from the final attempt. Callers handle
#'   status codes exactly as before -- this only changes how many attempts
#'   happen before returning.
api_http_get_with_retry <- function(url, headers) {
  attempt <- 0L
  repeat {
    attempt <- attempt + 1L
    resp <- tryCatch(
      httr::GET(url, headers, httr::timeout(api_timeout_sec)),
      error = function(e) e
    )
    is_network_error <- inherits(resp, "error")
    is_retryable_status <- !is_network_error && httr::status_code(resp) %in% c(502L, 503L, 504L)

    if (!is_network_error && !is_retryable_status) {
      return(resp)
    }
    if (attempt >= api_max_retries) {
      if (is_network_error) {
        stop("Quickomics API: request to ", url, " failed after ", api_max_retries,
             " attempts: ", conditionMessage(resp))
      }
      return(resp)  # exhausted retries on a persistent 502/503/504 -- let the caller report it
    }

    backoff <- api_retry_backoff_sec * (2 ^ (attempt - 1L))
    reason <- if (is_network_error) conditionMessage(resp) else paste0("HTTP ", httr::status_code(resp))
    message("Quickomics API: request to ", url, " did not complete (", reason, "), attempt ",
            attempt, "/", api_max_retries, " -- retrying in ", backoff, "s...")
    Sys.sleep(backoff)
  }
}

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
#' Only actually attempts anything if QUICKOMICS_API_AUTOSTART=true --
#' see the comment on api_compose_dir above for why this must stay off
#' when QuickOmics and the API are on different servers.
.api_autostart_notification_id <- "quickomics_api_autostart"

ensure_api_running <- function(max_wait_sec = 30, poll_interval_sec = 2) {
  if (api_health_check()) return(invisible(TRUE))
  
  if (!api_autostart_on) {
    msg <- paste0(
      "Quickomics API not reachable at ", api_base_url,
      ". Please make sure the API settings (QUICKOMICS_API_URL) are correct ",
      "and the API is running."
    )
    showNotification(msg, type = "error", duration = NULL)
    stop(
      "Quickomics API is not reachable at ", api_base_url, " and auto-start ",
      "is disabled. Start it manually (`docker compose up -d` in the API ",
      "project directory, on whichever server actually runs it), or set ",
      "QUICKOMICS_API_AUTOSTART=true and QUICKOMICS_API_COMPOSE_DIR to ",
      "enable automatic startup (only meaningful when this R process runs ",
      "on the SAME server as the API's docker-compose.yml)."
    )
  }
  
  if (.api_autostart_cache$attempted) {
    msg <- paste0(
      "Quickomics API is still not reachable at ", api_base_url,
      " after a previous auto-start attempt this session. Not retrying automatically."
    )
    showNotification(msg, type = "error", duration = NULL)
    stop(
      "Quickomics API is still not reachable at ", api_base_url,
      " after a previous auto-start attempt this session. Not retrying ",
      "automatically -- check `docker compose logs` in ", api_compose_dir,
      " for why the container isn't coming up."
    )
  }
  .api_autostart_cache$attempted <- TRUE
  
  if (!nzchar(api_compose_dir)) {
    msg <- paste0(
      "Quickomics API not reachable at ", api_base_url,
      ". Auto-start is enabled but QUICKOMICS_API_COMPOSE_DIR is not set -- ",
      "please check the API settings."
    )
    showNotification(msg, type = "error", duration = NULL)
    stop(
      "QUICKOMICS_API_AUTOSTART is enabled but QUICKOMICS_API_COMPOSE_DIR ",
      "is not set -- don't know where to find docker-compose.yml."
    )
  }
  compose_file <- file.path(api_compose_dir, "docker-compose.yml")
  if (!file.exists(compose_file)) {
    msg <- paste0(
      "Quickomics API not reachable at ", api_base_url,
      ". Auto-start is enabled but no docker-compose.yml was found in '",
      api_compose_dir, "' -- please check the API settings."
    )
    showNotification(msg, type = "error", duration = NULL)
    stop("QUICKOMICS_API_COMPOSE_DIR is set to '", api_compose_dir,
         "' but no docker-compose.yml was found there.")
  }
  
  showNotification(
    paste0("Quickomics API not reachable at ", api_base_url,
           " -- attempting to start it automatically..."),
    id = .api_autostart_notification_id, type = "warning", duration = NULL
  )
  message("Quickomics API not reachable at ", api_base_url,
          " -- attempting `docker compose up -d` in ", api_compose_dir, " ...")
  
  out <- tryCatch(
    system2("docker", c("compose", "-f", shQuote(compose_file), "up", "-d"),
            stdout = TRUE, stderr = TRUE),
    error = function(e) {
      msg <- paste0("Failed to start the Quickomics API: ", conditionMessage(e),
                    ". Does this R process have permission to run docker commands?")
      showNotification(msg, id = .api_autostart_notification_id, type = "error", duration = NULL)
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
      showNotification(
        paste0("Quickomics API is now reachable at ", api_base_url,
               " (took ~", waited, "s)."),
        id = .api_autostart_notification_id, type = "message", duration = 5
      )
      message("Quickomics API is up (took ~", waited, "s).")
      return(invisible(TRUE))
    }
    showNotification(
      paste0("Quickomics API not reachable at ", api_base_url,
             " -- still attempting to start it... (", waited, "s / ",
             max_wait_sec, "s)"),
      id = .api_autostart_notification_id, type = "warning", duration = NULL
    )
  }
  
  showNotification(
    paste0("Quickomics API is still not reachable at ", api_base_url,
           " after attempting to start it (waited ", max_wait_sec, "s). ",
           "Please check the API settings."),
    id = .api_autostart_notification_id, type = "error", duration = NULL
  )
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

#' Fetch a project's registry metadata (Name, Species, ShortName,
#' ExpressionUnit, Path) from the API via GET /datasets/{source}/{id} --
#' this is the API's read_project_meta(), which reads the project's .csv
#' manifest SERVER-SIDE, where server_dir/test_dir are always valid paths,
#' regardless of what machine QuickOmics itself is running on.
#'
#' This replaces reading server_dir/test_dir's .csv manifest directly via
#' read.csv() in inputdata.R's serverfile/testfile branches -- that local
#' read only works if QuickOmics happens to share a filesystem with the
#' API's network storage, which isn't guaranteed (and defeats the point of
#' going through the API at all for everything else).
#'
#' @param project_id Project ID
#' @param source "published" or "unpublished"
#' @return A list with names matching what ProjectInfo expects directly:
#'   Name, Species, ShortName, ExpressionUnit, Path (any of the latter two
#'   may be NA if the .csv manifest doesn't have that column). Returns NULL
#'   (not an error) if the API responds 404 -- i.e. no such project -- so
#'   callers can use the same `validate(need(...))` pattern as the local
#'   file.exists() check it replaces.
api_get_project_meta <- function(project_id, source) {
  ensure_api_running()
  
  fetch_once <- function(token) {
    api_http_get_with_retry(
      paste0(api_base_url, "/datasets/", source, "/", project_id),
      httr::add_headers(Authorization = paste("Bearer", token))
    )
  }
  
  token <- api_get_token()
  resp <- fetch_once(token)
  
  if (httr::status_code(resp) == 401) {
    token <- api_fetch_token()
    resp <- fetch_once(token)
  }
  
  if (httr::status_code(resp) == 404) {
    return(NULL)
  }
  
  if (httr::status_code(resp) != 200) {
    stop(
      "Quickomics API: failed to fetch project metadata for '", project_id,
      "' (source='", source, "'), HTTP ", httr::status_code(resp),
      ". Response: ", httr::content(resp, "text", encoding = "UTF-8")
    )
  }
  
  parsed <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"))
  unwrap <- function(x) if (is.list(x) || length(x) > 1) x[[1]] else x
  
  list(
    Name            = unwrap(parsed$name),
    Species         = unwrap(parsed$species),
    ShortName       = unwrap(parsed$short_name),
    ExpressionUnit  = { v <- unwrap(parsed$expression_unit); if (is.null(v) || is.na(v)) NULL else v },
    Path            = { v <- unwrap(parsed$path);            if (is.null(v) || is.na(v)) NULL else v }
  )
}

#' Load a project's full data bundle from the Quickomics API, assigning the
#' same object names base::load(RDataFile) would have created directly
#' into `envir` -- so existing code written against `objs <- load(RDataFile)`
#' and bare object names (MetaData, ProteinGeneName, data_wide, ...) keeps
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
    api_http_get_with_retry(
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
      # serialize()/unserialize() round trip cleanly -- setDT() is
      # data.table's own official remedy; no-op for anything that isn't
      # actually a data.table.
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
#' assigns the result to a local `network` variable itself.
#'
#' @param project_id Project ID
#' @param source "published" or "unpublished"
#' @return The network data.frame, or NULL if the API responds 404 (no
#'   precomputed network for this project) -- same meaning as a missing
#'   local file, not an error; callers should fall through to computing
#'   fresh, not treat this as a failure.
api_load_network <- function(project_id, source) {
  ensure_api_running()
  
  fetch_once <- function(token) {
    api_http_get_with_retry(
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
  if (data.table::is.data.table(network)) {
    data.table::setDT(network)
  }
  network
}

#' Load a project's WGCNA bundle (netwk, dataExpr, picked_power, sft) from
#' the API, assigning the same object names that
#' `load(wgcnafile); load(load_wgcna_file)` would have created directly
#' into `envir`.
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
    api_http_get_with_retry(
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
