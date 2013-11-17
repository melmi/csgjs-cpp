# About
This is a fork of [csgjs-cpp](https://github.com/dabroz/csgjs-cpp). The differences between this fork and the original version are:

* This fork uses `double` instead of `float`.
* This fork separates header and source files.
* The only field of `vertex` struct is `pos`.
* The functions `union`, `intersection` and `difference` renamed to `get_xxx`.
* Some geometric calculations are planned to be added.


csgjs-cpp
=========

CSG library for C++, port of https://github.com/evanw/csg.js/
