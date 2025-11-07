{ lib, stdenv, fetchgit, cmake, gsl, gmp, makeWrapper, jemalloc, htslib, git, zlib, pkg-config }:

stdenv.mkDerivation rec {
  pname = "seqwish";
  version = "0.7.9";

  src = fetchgit {
    url = "https://github.com/ekg/seqwish.git";
    rev = "f44b402f0c2e02988d431d9b2e5eba9727cf93a9";
    sha256 = "sha256-xEdiDyFwUYOQI7v4WH8Kv3UXQuBLIS3IjQi8Kv2ff9Q=";
    fetchSubmodules = true;
  };

  # Add -march=native to the CFLAGS
  NIX_CFLAGS_COMPILE = "-march=native";

  nativeBuildInputs = [ cmake makeWrapper ];

  buildInputs = [ 
    gsl 
    gmp
    jemalloc
    htslib
    git
    zlib
    pkg-config
  ];

  postInstall = ''
    wrapProgram $out/bin/seqwish --prefix PATH : ${lib.makeBinPath [ gsl gmp ]}
  '';

  meta = with lib; {
    description = "alignment to variation graph inducer";
    homepage = "https://github.com/ekg/seqwish";
    license = licenses.mit;
    platforms = platforms.linux;
  };
}
