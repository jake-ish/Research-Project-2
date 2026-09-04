# How to create the private repo and push this

This file is for your own use — you can delete it before or after publishing.

## 0. Before you push: safety checks

- [ ] Add your actual scripts into `src/kaks/` and `src/modelling/`.
- [ ] Open every script and remove hard-coded local paths
      (e.g. `/Users/you/...`), API keys, and any credentials.
- [ ] Confirm no raw third-party data is staged (it's covered by `.gitignore`,
      but double-check `git status` before your first commit).
- [ ] Run `pip freeze > requirements.txt` in your working env to pin versions.

## 1. Create the repository on GitHub (private)

Either via the website — New repository → name it `mutation-rate-architecture`
→ set **Private** → do NOT initialise with README/.gitignore/license (you already
have them) → Create.

Or with the GitHub CLI (`gh`), from inside this folder after step 2:

```bash
gh repo create mutation-rate-architecture --private --source=. --remote=origin
```

## 2. Initialise git locally and make the first commit

From inside the unzipped `mutation-rate-architecture/` folder:

```bash
git init
git add .
git status          # <-- review carefully; make sure no raw data / secrets listed
git commit -m "Initial commit: analysis code, processed data, and docs"
```

## 3. Connect to GitHub and push (if you used the website in step 1)

```bash
git branch -M main
git remote add origin https://github.com/<your-username>/mutation-rate-architecture.git
git push -u origin main
```

(If you used `gh repo create` in step 1, it already set the remote — just
`git push -u origin main`.)

## 4. Sharing with your supervisor

While the repo is private, add them as a collaborator:
Settings → Collaborators → Add people → their GitHub username.

## 5. When the paper is accepted

- Flip the repo to public (Settings → General → Danger Zone → Change visibility).
- Consider archiving a release on Zenodo to get a citable DOI, and add that DOI
  to the paper's data-availability statement.
