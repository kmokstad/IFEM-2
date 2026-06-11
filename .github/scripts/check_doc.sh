#!/bin/bash
# Checks out previous version of the doc from the gh-pages branch, and
# compares with the newly generated doc while ignoring generated time stamps.

git checkout gh-pages
git rm -rf docs
mv doc/html docs
git add docs
if git diff -I "^Generated on .* for IFEM" HEAD --quiet; then
  echo No changes in source code documentation
  git reset --hard HEAD
  exit 0
else
  exit 1
fi
