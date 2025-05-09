{
  description = "Flake to get compartmap development environment";
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system: let
      pkgs = nixpkgs.legacyPackages.${system};

      Depends = with pkgs.rPackages; [
        BiocSingular
        bsseq
        csaw
        HDF5Array
        RaggedExperiment
        SummarizedExperiment
      ];

      Imports = with pkgs.rPackages; [
        DelayedArray
        DelayedMatrixStats
        GenomicRanges
        ggplot2
        impute
        Matrix
        reshape2
        RMTstat
        rtracklayer
        scales
      ];

      Suggests = with pkgs.rPackages; [
        BiocStyle
        knitr
        markdown
        minfi
        minfiData
        Rcpp
        rmarkdown
        roxygen2
        SingleCellExperiment
        testthat
      ];

      rDevDeps = with pkgs.rPackages; [
        biocthis
        covr
        devtools
        DT
        htmltools
        styler
        usethis
      ];

      sysDeps = with pkgs; [
        R
      ];

      tex = (pkgs.texlive.combine {
        inherit (pkgs.texlive) scheme-medium
        inconsolata
        xkeyval
        etoolbox;
      });
      sysDevDeps = with pkgs; [
        air-formatter
        checkbashisms
        html-tidy
        tex
      ];

      # default package
      rDeps = [ Depends Imports Suggests ];
      compartmap = pkgs.rPackages.buildRPackage {
        name = "compartmap";
        src = self;
        nativeBuildInputs = sysDeps;
        propagatedBuildInputs = rDeps;
      };
      # Create R development environment with compartmap and other useful libraries
      rvenv = pkgs.rWrapper.override {
        packages = rDeps ++ rDevDeps ++ sysDeps ++ sysDevDeps;
      };
    in {
      packages.default = compartmap;
      devShells.default = pkgs.mkShell {
          buildInputs = rDeps ++ rDevDeps ++ sysDeps ++ sysDevDeps;
          inputsFrom = pkgs.lib.singleton compartmap;
          packages = pkgs.lib.singleton rvenv;
      };
    });
}

