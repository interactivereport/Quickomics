###########################################################################################################
## Regenerates BUILD_DATE.txt from the current git commit.
##
## Run this locally right before copying the app folder to a production
## deployment that isn't itself a git checkout. BUILD_DATE.txt travels with
## the copy (plain files, unlike .git which usually doesn't get copied), so
## global.R can still show a reliable "last updated" date in the footer even
## when the production folder has no .git directory of its own.
##
## BUILD_DATE.txt is gitignored on purpose -- it's a per-deploy stamp, not
## something that should be committed (it'd go stale the moment the next
## commit lands).
##
## Usage: Rscript update_build_date.R
###########################################################################################################

git_field <- function(args) {
  out <- suppressWarnings(system2("git", args, stdout = TRUE, stderr = FALSE))
  if (length(out) != 1 || !nzchar(out)) {
    stop("`git ", paste(args, collapse = " "), "` returned no output -- is this run from inside the repo?")
  }
  out
}

commit_hash <- git_field(c("log", "-1", "--format=%H"))
commit_date <- git_field(c("log", "-1", "--format=%cd", "--date=format:%Y-%m-%d"))
branch      <- git_field(c("rev-parse", "--abbrev-ref", "HEAD"))

writeLines(
  c(
    paste0("commit: ", commit_hash),
    paste0("date: ", commit_date),
    paste0("branch: ", branch)
  ),
  "BUILD_DATE.txt"
)

cat("Wrote BUILD_DATE.txt:\n")
cat("  commit: ", commit_hash, "\n", sep = "")
cat("  date:   ", commit_date, "\n", sep = "")
cat("  branch: ", branch, "\n", sep = "")
