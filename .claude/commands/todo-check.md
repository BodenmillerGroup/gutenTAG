---
description: Display all open tasks from TODO.md as a grouped, prioritised table; migrate any newly completed tasks to COMPLETED.md
allowed-tools: [Read, Edit, Bash]
---

## Step 1 — Migrate completed tasks

Before displaying the task table, scan `TODO.md` for any completed tasks that have not yet been moved to `COMPLETED.md`:
- Struck-through items: `~~…~~`
- Checked checkboxes: `- [x] …`

For each such item found:
1. Determine the current git branch with `git branch --show-current` and the package version from `DESCRIPTION` (`Version:` field).
2. Append the item (plain text, no markdown strikethrough or checkbox syntax) to `COMPLETED.md` under a heading `## Version <version> — Branch: <branch>`. If that exact heading already exists, append under it. If it does not exist, create it at the bottom of the file.
3. Remove the completed item (and its sub-bullets if any) from `TODO.md`.

If no completed tasks are found, skip this step silently.

## Step 2 — Display open tasks

Read the file `TODO.md` in the project root and display all incomplete (unchecked) tasks as a Markdown table.

Group the table by priority section in this order:
1. **General** — top-level bullet items (not under any heading) that are not struck through
2. **Priority 1** — bugs that break on clean install or produce silent wrong results
3. **Priority 2** — convention violations and R CMD CHECK failures
4. **Priority 3** — performance and style issues
5. **Priority 4** — low-priority cleanup
6. **README.md Updates** — grouped by their sub-section (Remove/correct, Add missing sections, Prose and formatting)

For each group, produce a table with columns:

| # | Task | Location |
|---|------|----------|

- **#**: sequential number within the group
- **Task**: the short task description (bold label only, no long explanation)
- **Location**: file/function reference if given, otherwise `—`

Omit any tasks that are already ticked (`[x]`) or struck through (`~~…~~`).

After the tables, print a one-line count: `X open tasks across Y groups.`

## Step 3 — Display complexity ranking

After the task table, display the following complexity ranking (simplest to most complex). Update it if tasks have been added, removed, or completed since this was last revised.

| Rank | Task | Why |
|------|------|-----|
| 1 | **Replace `== TRUE` / `== T` comparisons** | Mechanical find-and-replace across files |
| 2 | **Replace all mentions of "density" with "sparsity"** | Likely rename across a few files |
| 3 | **Remove dead and commented-out code blocks** | Straightforward deletion, low risk |
| 4 | **Fix `MulticoreParam` Windows incompatibility** | Pattern already exists in `peakDetection.R`, just copy it |
| 5 | **Add input validation to `generateSeedMz` for `density` parameter** | Small guard + roxygen entry |
| 6 | **Vectorise `cor` computation in `.finaliseCorrespondence`** | Single `cor()` call replacement |
| 7 | **Replace column loop in `.filterTIC` with index-based assignment** | Already resolved by bug fix — consider removing |
| 8 | **Replace vector-growing loops with pre-allocated vectors** | Mechanical but spread across multiple files |
| 9 | **Fix NaN propagation in `computeGearysC` for zero-variance channels** | Small guard, but need to verify downstream behaviour |
| 10 | **Replace `cbind` growing loop in `.generateFinalIntensityDF`** | Pre-allocation requires knowing dimensions upfront |
| 11 | **Add `shift = "closest"` option to `.removeDuplicates`** | Small addition but needs `choose_closest` helper and validation update |
| 12 | **Investigate `'x' should be a list object` error** | Needs reproduction and diagnosis first |
| 13 | **Determine minimum compatible Cardinal version** | Requires systematic testing across versions |
| 14 | **Add QC plot for mean-variance residuals** | New ggplot function, needs design decisions |
| 15 | **Remove/reassign untargeted metapeaks incorrectly assigned to markers** | Requires understanding of metapeak assignment logic |
| 16 | **Move discarded duplicate peaks into `Untargeted`** | Most involved — structural change to output, needs careful index tracking |
