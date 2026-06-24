#!/usr/bin/env bash
# PreToolUse(Bash) guard: refuse `git commit` while HEAD is a protected branch
# (main/master). Branch hygiene is otherwise prose-only here, and origin/main
# already carries direct non-PR commits (a2333b5, eed6f04, ...). The agent should
# branch first, then commit.
#
# Fails OPEN whenever it cannot prove the branch is protected (not a git repo,
# detached HEAD, git missing) -- a guard must never wedge unrelated Bash calls.
# Scope: only the agent's own `git commit` calls in a Claude session; an external
# auto-committer is a separate path this hook cannot see.
set -uo pipefail

payload="$(cat)"
cmd="$(printf '%s' "$payload" | python3 -c 'import json,sys; print(json.load(sys.stdin).get("tool_input",{}).get("command",""))' 2>/dev/null || true)"

# Only real commits matter. Substring match is intentional (ponytail): the branch
# check below is the actual gate, so a stray "git commit" in some other command is
# harmless unless we are genuinely on a protected branch.
case "$cmd" in
  *"git commit"*) ;;
  *) exit 0 ;;
esac

branch="$(git symbolic-ref --short HEAD 2>/dev/null || true)"
[ -z "$branch" ] && exit 0  # detached HEAD / not a repo -> fail open

case "$branch" in
  main | master)
    printf 'refuse: commit on protected branch "%s". Branch first (git switch -c <feature-branch>), then commit.\n' "$branch" >&2
    exit 2
    ;;
esac
exit 0
