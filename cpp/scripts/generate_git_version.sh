#!/bin/bash
cd "$(dirname "$0")/.."
GIT_VERSION="$(git describe --always --tags --long 2>/dev/null || echo 'unknown')"
echo "#ifndef SEQWISH_GIT_VERSION_H"
echo "#define SEQWISH_GIT_VERSION_H"
echo "#define SEQWISH_GIT_VERSION \"${GIT_VERSION}\""
echo "#endif"
