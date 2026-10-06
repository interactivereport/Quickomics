###########################################################################################################
## Leveled logging for the Shiny server process.
##
## Writes to stdout via cat() -- deliberately not to a dedicated file. Both
## RStudio's console and Shiny Server (which captures each app process's
## stdout/stderr into its own per-session log under /var/log/shiny-server/)
## already give every running session its own log destination; duplicating
## that with our own file management would just be redundant and, on a
## deployment where the app directory itself isn't writable by the server
## process, another thing that can fail. Each line is still tagged with a
## short session id so multiple sessions sharing one console (e.g. during
## local testing) stay distinguishable.
##
## Usage: log_message("INFO", "Project loaded: ", ProjectID)
##    or: log_info("Project loaded: ", ProjectID)
##
## ---------------------------------------------------------------------
## Which level to use -- ask: would I want to see this on every normal
## run, only when something's actually wrong, or only while actively
## digging into one specific problem?
##
## ERROR -- an operation FAILED and could not complete as intended.
##   Almost always corresponds to something the user also saw (an error
##   notification, a blank plot, a stuck spinner). If you're inside a
##   tryCatch(..., error = function(e) {...}) catching a real exception,
##   that's almost always an ERROR. Examples: a required file/RData is
##   missing and the feature can't proceed; a computation that was
##   supposed to succeed threw; Save Session restore failed to reload a
##   project.
##
## WARN -- something unexpected happened, but the app RECOVERED and kept
##   going; nothing crashed, but a human should know. Examples: a
##   validate(need(...)) guard rejecting bad input (a bad ProjectID in a
##   URL, a file that fails a format check); the Save Session
##   project-mismatch guard rejecting an upload; falling back to a
##   default because expected config was missing; data-quality issues
##   silently worked around (NAs filtered, a factor coerced, an
##   unmatched gene ID dropped).
##
## INFO -- the normal narrative of what the app is doing: expected,
##   successful events worth a record even when nothing is wrong. Think
##   "a user says something went wrong an hour ago -- what did they
##   actually do?". Examples: session start/end, a project loaded,
##   Compute/Refresh finishing, a long computation (network/WGCNA build)
##   completing, a file uploaded and accepted, Save Session save/restore.
##
## DEBUG -- fine-grained detail only useful while actively chasing one
##   specific bug; too noisy to leave on by default. Examples: function
##   entry with parameter values, intermediate values in a multi-step
##   computation, which branch of an if/else was taken, the full detail
##   behind a WARN/ERROR that's too verbose for that line itself.
##
## Rule of thumb for this app's shape (most modules are an
## observeEvent()/eventReactive() doing real work, often guarded by
## validate()/tryCatch()): INFO on the start/success boundary of
## anything that does real work (a computation, a file operation, a
## user-triggered action); WARN inside a validate(need(...)) or other
## gentle-failure branch; ERROR inside a tryCatch(error=...) catching a
## real exception; DEBUG for whatever detail you'd only want while
## staring at that one bug.
## ---------------------------------------------------------------------
###########################################################################################################

LOG_LEVELS <- c(DEBUG = 1L, INFO = 2L, WARN = 3L, ERROR = 4L)

# Global threshold: messages below this level are skipped entirely (not
# just hidden -- the cat() never runs), so leaving DEBUG-level calls in the
# code costs nothing once log_level is raised in production. Configured via
# .Renviron, same pattern as the QUICKOMICS_API_* settings.
log_level <- {
  lvl <- toupper(Sys.getenv("QUICKOMICS_LOG_LEVEL", "INFO"))
  if (!lvl %in% names(LOG_LEVELS)) {
    warning("Unknown QUICKOMICS_LOG_LEVEL '", lvl, "' -- falling back to INFO. ",
            "Valid levels: ", paste(names(LOG_LEVELS), collapse = ", "))
    lvl <- "INFO"
  }
  lvl
}

#' Write a leveled log message for the current session, if log_level allows it.
#'
#' @param level One of "DEBUG", "INFO", "WARN", "ERROR" (case-insensitive).
#' @param ... Concatenated (via paste0) to build the message text -- so
#'   call sites can pass multiple pieces directly, e.g.
#'   log_message("INFO", "Project loaded: ", ProjectID, " (", nrow(df), " rows)").
#' @param session The Shiny session this message belongs to; defaults to
#'   whichever session is currently active (works from the main server
#'   function or any moduleServer(), no need to thread `session` through
#'   call sites manually). NULL (e.g. logging from global.R at app startup,
#'   before any session exists) is handled gracefully.
log_message <- function(level, ..., session = shiny::getDefaultReactiveDomain()) {
  level <- toupper(level)
  # LOG_LEVELS is an atomic vector, not a list -- unlike a list, `[[` on an
  # atomic vector throws "subscript out of bounds" for a missing name
  # instead of returning NULL, so the name has to be checked with %in%
  # first rather than relying on an is.null() check after the fact.
  if (!level %in% names(LOG_LEVELS)) {
    stop("Unknown log level '", level, "'. Valid levels: ", paste(names(LOG_LEVELS), collapse = ", "))
  }
  lvl_num <- LOG_LEVELS[[level]]

  # A session can raise (or lower) its own verbosity at runtime via the Log
  # Level control in the Output tab, without restarting the R process --
  # useful for reproducing an intermittent problem with DEBUG on, then
  # putting it back. Falls back to the app-wide default (log_level, set
  # once from QUICKOMICS_LOG_LEVEL at startup) for any session that hasn't
  # touched that control.
  effective_level <- log_level
  if (!is.null(session) && !is.null(session$userData$log_level)) {
    effective_level <- session$userData$log_level
  }
  if (lvl_num < LOG_LEVELS[[effective_level]]) {
    return(invisible(NULL))
  }

  session_id <- if (!is.null(session)) substr(session$token, 1, 8) else "no-session"

  # Never let a logging call itself take down the app (e.g. a non-character
  # argument that paste0() can't coerce cleanly).
  tryCatch({
    cat(sprintf(
      "[%s] [%-5s] [session:%s] %s\n",
      format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"),
      level,
      session_id,
      paste0(..., collapse = "")
    ))
  }, error = function(e) {
    cat(sprintf("[%s] [WARN ] [session:%s] log_message() itself failed: %s\n",
                format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), session_id, conditionMessage(e)))
  })
  invisible(NULL)
}

#' Override the log level for one session only, at runtime.
#'
#' Lets a user (or whoever is helping them) raise verbosity to reproduce a
#' problem live and check the console/Shiny Server log afterward, without
#' restarting the R process -- which would also lose whatever app state was
#' needed to reproduce it in the first place.
#'
#' @param session The session to set this for (required -- unlike
#'   log_message()/log_info() etc., this has no sensible session-less use,
#'   so it's not defaulted via getDefaultReactiveDomain()).
#' @param level One of "DEBUG", "INFO", "WARN", "ERROR".
set_session_log_level <- function(session, level) {
  level <- toupper(level)
  if (!level %in% names(LOG_LEVELS)) {
    stop("Unknown log level '", level, "'. Valid levels: ", paste(names(LOG_LEVELS), collapse = ", "))
  }
  session$userData$log_level <- level
  # Logged unconditionally (not through the normal threshold check) so the
  # level change itself always shows up in the log, even when switching
  # to a stricter level than before.
  session_id <- substr(session$token, 1, 8)
  cat(sprintf("[%s] [INFO ] [session:%s] Log level changed to %s\n",
              format(Sys.time(), "%Y-%m-%d %H:%M:%OS3"), session_id, level))
  invisible(NULL)
}

log_debug <- function(..., session = shiny::getDefaultReactiveDomain()) log_message("DEBUG", ..., session = session)
log_info  <- function(..., session = shiny::getDefaultReactiveDomain()) log_message("INFO",  ..., session = session)
log_warn  <- function(..., session = shiny::getDefaultReactiveDomain()) log_message("WARN",  ..., session = session)
log_error <- function(..., session = shiny::getDefaultReactiveDomain()) log_message("ERROR", ..., session = session)
