{ pkgs ? import <nixpkgs> {} }:

pkgs.callPackage ./seqwish.nix {
  inherit (pkgs) stdenv fetchgit cmake gsl gmp makeWrapper jemalloc htslib git zlib pkg-config;
}