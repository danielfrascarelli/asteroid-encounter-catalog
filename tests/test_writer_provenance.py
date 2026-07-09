"""Tests for provenance capture in the catalog writer.

Focus on the git-commit fallback that reads ``.git/HEAD`` directly, exercised
when the pipeline runs inside a container that ships no git binary (C16).
"""

from __future__ import annotations

import subprocess
from pathlib import Path

from src.catalog.writer import _git_commit_from_dotgit


def _make_repo(root: Path, ref: str, sha: str, *, packed: bool = False) -> None:
    """Create a minimal fake ``.git`` layout resolving HEAD -> sha."""
    git = root / ".git"
    git.mkdir()
    (git / "HEAD").write_text(f"ref: {ref}\n")
    if packed:
        (git / "packed-refs").write_text(f"# pack-refs with: peeled\n{sha} {ref}\n")
    else:
        ref_file = git / ref
        ref_file.parent.mkdir(parents=True, exist_ok=True)
        ref_file.write_text(sha + "\n")


def test_dotgit_loose_ref(tmp_path: Path, monkeypatch) -> None:
    sha = "a" * 40
    _make_repo(tmp_path, "refs/heads/main", sha)
    monkeypatch.chdir(tmp_path)
    assert _git_commit_from_dotgit() == sha


def test_dotgit_packed_ref(tmp_path: Path, monkeypatch) -> None:
    sha = "b" * 40
    _make_repo(tmp_path, "refs/heads/feature", sha, packed=True)
    monkeypatch.chdir(tmp_path)
    assert _git_commit_from_dotgit() == sha


def test_dotgit_detached_head(tmp_path: Path, monkeypatch) -> None:
    sha = "c" * 40
    git = tmp_path / ".git"
    git.mkdir()
    (git / "HEAD").write_text(sha + "\n")
    monkeypatch.chdir(tmp_path)
    assert _git_commit_from_dotgit() == sha


def test_dotgit_walks_up_to_parent(tmp_path: Path, monkeypatch) -> None:
    sha = "d" * 40
    _make_repo(tmp_path, "refs/heads/main", sha)
    sub = tmp_path / "src" / "deep"
    sub.mkdir(parents=True)
    monkeypatch.chdir(sub)
    assert _git_commit_from_dotgit() == sha


def test_dotgit_absent(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.chdir(tmp_path)
    assert _git_commit_from_dotgit() == ""


def test_matches_real_git_here() -> None:
    """In the repo (with a git binary), the fallback matches ``git rev-parse``."""
    try:
        real = subprocess.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return  # no git binary available; nothing to compare against
    assert _git_commit_from_dotgit() == real
