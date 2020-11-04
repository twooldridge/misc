## Super basic git guide
Only really suitable as a reminder for working with myself and updating my github from the cluster.
There are probably some bad practices here, but it accomplishes what I need.

Steps:
1. Before doing anything, check the current branch with `git branch`
2. In the repo you want to push, type `git add .` -- this will add changes made to all files in the folder.
3. Next commit, and try to use an informative message: `git commit -m 'new stuff'`
4. If your local branch is `master`, then skip this step. If not, do `git merge <branch> master`
5. Finally, update the remote repo with `git push origin master`
