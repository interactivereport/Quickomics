# Applying the api-data-source branch

This was built against a fresh clone of https://github.com/interactivereport/Quickomics
(branch `master` at the time of writing), not your actual deployed checkout —
I don't have access to wherever your real QuickOmics install lives. Pick
whichever option below fits how your checkout is set up.

## Option 1 — Apply the patch directly (preferred if your checkout is close to upstream)

```bash
cd /path/to/your/quickomics/checkout
git checkout -b api-data-source
git am /path/to/0001-api-data-source.patch
```

If `git am` complains about conflicts (likely if your checkout has diverged
from upstream `master` at all), fall back to Option 2.

## Option 2 — Apply the plain diff, resolve by hand

```bash
cd /path/to/your/quickomics/checkout
git checkout -b api-data-source
git apply --reject /path/to/full.diff
```
`--reject` will apply whatever hunks cleanly match and drop `.rej` files for
anything that doesn't — check for any `.rej` files afterward and merge them
in manually.

## Option 3 — Just copy the three files over (simplest, least "gitty")

The three files in this folder (`api_client.R`, `global.R`, `inputdata.R`)
are the complete post-change versions, not diffs. If your checkout's
`inputdata.R`/`global.R` haven't diverged much from upstream, you can just
create the branch and drop these in directly:

```bash
cd /path/to/your/quickomics/checkout
git checkout -b api-data-source
cp /path/to/api_client.R .
cp /path/to/global.R .        # only if yours hasn't diverged -- check first!
cp /path/to/inputdata.R .     # only if yours hasn't diverged -- check first!
git add api_client.R global.R inputdata.R
git commit -m "Route serverfile/testfile projects through the Quickomics API"
```

**Careful with `global.R` and `inputdata.R`** if your actual deployment has
already been customized away from upstream `master` -- copying wholesale
would silently discard any of those local changes. Check with `diff` against
your real files first, or just use Option 1/2 instead.

## Required environment variables (wherever QuickOmics actually runs)

```
QUICKOMICS_API_URL=http://<api-host>:8499
QUICKOMICS_API_USERNAME=admin
QUICKOMICS_API_PASSWORD=<the real password, not the placeholder from earlier testing>
```

### Optional: auto-start the API container

Off by default. If enabled, `api_load()` will health-check the API before
each use and, if it's down, run `docker compose up -d` and wait for it to
come up:

```
QUICKOMICS_API_AUTOSTART=true
QUICKOMICS_API_COMPOSE_DIR=/path/to/QuickomicsAPI    # dir containing docker-compose.yml
```

**This requires the R process running QuickOmics to have permission to run
`docker` commands** (the OS user needs to be in the `docker` group, or have
passwordless sudo configured for it). It's a real trade-off -- application
code gets to control infrastructure lifecycle, which is normally kept
separate. The conventional alternative is simpler: just leave the API
container always running (`docker-compose.yml` already has
`restart: unless-stopped`), managed independently of QuickOmics. Auto-start
exists purely as a dev/testing convenience (e.g. running QuickOmics from
RStudio without remembering to start the API container first) -- leave it
off in any shared/production deployment unless you've deliberately decided
the trade-off is worth it there too.

## New R package dependencies

`httr` and `jsonlite` -- confirmed not used anywhere else in this codebase
today, so make sure they're installed wherever QuickOmics' R environment
actually lives (this repo has no Dockerfile/renv.lock/DESCRIPTION, so I
can't tell what manages that for you -- check however your team currently
installs QuickOmics' R package dependencies).

## What did NOT change

Every data-loading path other than `?api_serverfile=`/`?api_testfile=`/`?datalake=` -- the local
`data/` folder, the `unlisted/` folder, `?serverfile=`/`?testfile=` (read
directly from `server_dir`/`test_dir`), and user file uploads -- is
completely untouched and still uses `base::load()` exactly as before.

## Network.R: one deliberate behavior change for API-sourced projects

For projects loaded via `?api_serverfile=`/`?api_testfile=`/`?datalake=`, `network.R` no longer
writes a `_network.RData` cache file back to disk after computing a fresh
network -- there's no local network storage for an API-sourced project to
write back to. Every other behavior (loading a cached network for the full
sample set, recomputing live for a filtered subset, falling back to live
computation if no precomputed network exists) is unchanged and matches the
original local-file behavior exactly, including for API-sourced projects
whose `_network.RData` doesn't exist on the API side (a 404 there is
treated the same as a missing local file -- not an error).

## Wgcna.R: no behavior changes needed, simpler swap than network.R

Unlike `network.R`, `wgcna.R` never writes a cache file back to disk at all
(regardless of source), so there was no write-back decision to make here --
all three call sites that used `load(wgcnafile); load(load_wgcna_file)` now
branch on `ProjectInfo$source` and call the new `api_load_wgcna()` instead,
with no other behavior change. A missing/failed WGCNA bundle on the API
side is treated the same as a missing local file at each call site.

## A real bug found and fixed along the way: data.table serialization

`data_results` (from `inputdata.R`) and `network` (from `network.R`) both
hit a known `data.table` quirk -- the internal self-reference pointer
(`.internal.selfref`) isn't guaranteed to survive a `serialize()`/
`unserialize()` round trip cleanly, which is exactly what happens when data
crosses the API boundary. This surfaced as `Error in <Anonymous>: formal
argument "dataframe" matched by multiple actual arguments` when a DT table
tried to render an affected object. Fixed in `api_client.R` by calling
`data.table::setDT()` (the package's own official remedy) on every
data.table-classed object right after `unserialize()`, in `api_load()`,
`api_load_network()`, and `api_load_wgcna()` alike.

## How to test

1. Start QuickOmics with the env vars above set.
2. Visit `?api_serverfile=TST11809` (or any real published project ID).
3. Confirm the project loads, the summary/sample table/results table all
   populate correctly, and (if you have server-side logs) that a request to
   `/datasets/published/TST11809/bundle.rds` actually hit the API rather
   than falling through to a local file read.
4. Repeat with `?api_testfile=<some unpublished project>` to confirm the
   unpublished path also works.
5. Confirm a project loaded via the normal dropdown / `?project=` (local
   data/ folder) still works exactly as before -- this path should show
   zero behavior change.
