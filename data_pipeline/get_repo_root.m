function repoRoot = get_repo_root()
%GET_REPO_ROOT Return the absolute path to the repository root.

thisFile = mfilename('fullpath');
repoRoot = fileparts(fileparts(thisFile));
end
