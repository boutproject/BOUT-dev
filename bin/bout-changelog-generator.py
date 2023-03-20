#! /usr/bin/env python3

import argparse
import os
import github


def format_pull(pull: github.PullRequest.PullRequest) -> str:
    """Format a PR as a nice changelog line"""
    url = pull.url.replace("api.", "").replace("pulls", "pull").replace("repos/", "")
    user_url = pull.user.url.replace("api.", "").replace("users/", "")

    return rf"- {pull.title} [\#{pull.number}][{url}] ([{pull.user.login}][{user_url}])"


def pulls_after_tag(repo: github.Repository.Repository, tag: str, release_branch: str):
    """Get a list of PRs merged after a given tag"""

    release = repo.get_release(tag)

    pulls = []

    for pull in repo.get_pulls(
        state="closed", sort="created", base=release_branch, direction="desc"
    ):
        # Unmerged for whatever reason
        if not pull.merged_at:
            continue

        # We're now at PRs merged before the last release
        if pull.merged_at < release.created_at:
            break

        pulls.append(pull)

    return pulls


def pulls_into_RC(repo: github.Repository.Repository, tag: str):
    """Get PRs into the release candidate branch"""
    rc_branch = f"{tag}-rc"
    try:
        repo.get_branch(rc_branch)
    except github.GithubException:
        print(f"No release candidate branch named '{rc_branch}'")
        return []

    return repo.get_pulls(
        state="closed", sort="created", base=rc_branch, direction="desc"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser("BOUT++ changelog generator")
    parser.add_argument(
        "--token",
        "-t",
        help="GitHub API token, defaults to '${GITHUB_CHANGELOG_TOKEN}' environment variable",
        default=os.environ.get("GITHUB_CHANGELOG_TOKEN", ""),
    )
    parser.add_argument("last_tag", help="Last git tag")
    parser.add_argument("next_tag", help="Next git tag")
    parser.add_argument("--owner", help="Project organisation", default="boutproject")
    parser.add_argument("--repo", help="Repository name", default="BOUT-dev")
    parser.add_argument(
        "--release-branch",
        help="Release branch name (one of 'master' or 'next')",
        default="next",
    )

    args = parser.parse_args()

    repo_name = f"{args.owner}/{args.repo}"
    repo = github.Github(args.token).get_repo(repo_name)

    print(f"## [{args.next_tag}](https://github.com/{repo_name}/tree/{args.next_tag}")
    print(
        f"\n[Full Changelog](https://github.com/{repo_name}/compare/{args.last_tag}...{args.next_tag})\n"
    )

    pulls = []
    pulls.extend(pulls_into_RC(repo, args.next_tag))
    pulls.extend(pulls_after_tag(repo, args.last_tag, args.release_branch))

    for pull in pulls:
        print(format_pull(pull))
