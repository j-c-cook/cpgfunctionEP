# Contributing to cpgfunction

## Branching

To contribute code please first start an issue if one if not already started. 
Then:
- Branch from the most current `master` branch:
```
git checkout -b issuexx_ShortDescription  # checkout a branch locally from master
git push -u origin issuexx_ShortDescription  # push to your remote that you have a new branch
git branch --set-upstream-to=origin/issuexx_ShortDescription issuexx_ShortDescription
```
- Commit your changes and push them up
```
git add <filename>
git commit -m "Type 50 character message here"
git push origin issuexx_ShortDescription
```
- Create a pull request
  
- Document the changelog for your issue

The branch will then be merged into master once it has been reviewed. 

## Changelog

The cpgfunction library keeps a 
[changelog](https://github.com/j-c-cook/cpgfunction/blob/master/CHANGELOG.md)
so that changes upon each release are transparent and easily understood. Prior 
to a pull request being accepted, all changes must be marked in the changelog. 
The changes should fall under one of the following categories:

- New features - for new features
- Enhancements - for improvements made to code performance and functionality
- Maintenance - for tidying code
- Changed - for changes in functionality of the code
- Depracated - for soon-to-be removed features
- Removed - for removed features
- Fixes - for any bug fixes

## Versioning

This library makes use of [Semantic Versioning](https://semver.org/).
