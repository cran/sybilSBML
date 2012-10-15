%% sybilSBML Version History

\name{NEWS}
\title{sybilSBML News}
\encoding{UTF-8}

\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}

% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.15 2012-10-15}{
  \itemize{
    \item updated package vignette
    \item slot \code{rhs} of class \code{modelorg} is not set
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.14 2012-04-20}{
  \itemize{
    \item file \file{NEWS} is now in \file{inst/NEWS.Rd}
    \item removed arguments \code{verboseMode}, \code{loglevel}, \code{logfile}
          and \code{logfileEnc} from function \code{readSBMLmod}; use function
          \code{sink} instead
    \item fixed a bug in \code{.uglyHack}: tries no longer to fix SBMLWarnings,
          which is not necessary
    \item \code{.uglyHack} uses \code{tempfile()} to generate corrected
          SBML files 
    \item fixed a bug in \code{readSBMLmod()}: slot \code{subSys} was not of
          correct data type in rare cases
    \item better SBML id formatting in \code{readSBMLmod()}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.13 2012-03-16}{
  \itemize{
    \item compatible with \CRANpkg{sybil} version 0.5.12: replaced
          \CRANpkg{SparseM} with \CRANpkg{Matrix}

    \item some additions to \code{.uglyHack()}

    \item fixed a small bug in \code{readSBMLmod()}, when parsing sub systems
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.12 2011-12-12}{
  \itemize{
    \item fixed a small bug in \code{readSBMLmod()} concerning the detection of
          gpr rules

    \item fixed a bug in \code{readSBMLmod()}: genes from reactions which
          are excluded from the model are not read
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.11 2011-12-09}{
  \itemize{
    \item very minor adjustments in \code{readSBMLmod()}

    \item \code{readSBMLmod()} will continue reading the model,
          if \pkg{rsbml} produces warnings. Check warnings in the log file;
          arguments \code{loglevel} and \code{logfile} in \code{readSBMLmod()}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.10 2011-08-30}{
  \itemize{
    \item updated documentation
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.9 2011-08-30}{
  \itemize{
    \item added boolean arguments \code{remMet}, \code{contrMet},
          \code{singletonMet} and \code{deadEndmet} to \code{readSBMLmod()},
          the two boolen arguments \code{remSingleMet} and
          \code{constrSingleMet} are removed, see \code{help("readSBMLmod")}
          for details
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.8 2011-08-15}{
  \itemize{
    \item some small changes according to \CRANpkg{codetools}

    \item slots \code{met_single} and \code{react_single} of class
          \code{modelorg} are initialized with \code{NA}, not \code{FALSE}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.7 2011-08-12}{
  \itemize{
    \item new arguments for \code{readSBMLmod()}:
          \code{mergeMet}, \code{balanceReact}, \code{remUnusedMetReact},
          \code{remSingleMet}, \code{constrSingleMet} and \code{tol}.
          See \code{help("readSBMLmod")} for details
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.6 2011-08-08}{
  \itemize{
    \item a warning is produced, if a metabolite is used only once in S
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.5 2011-08-08}{
  \itemize{
    \item \code{readSBMLmod()} removes unused metabolites and reactions
          (rows and lines without any non-zero entries)

    \item \code{readSBMLmod()} merges and balances metabolites as in
          \code{readTSVmod()} from \CRANpkg{sybil}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.4 2011-08-05}{
  \itemize{
    \item fixed a bug in \code{readSBMLmod()} concerning the number of
          metabolite names and metabolite compartments
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.3 2011-08-05}{
  \itemize{
    \item added \code{NAMESPACE}

    \item new boolean parameter \code{bndCnd} in \code{readSBMLmod()}:
          If set to \code{TRUE}, external metabolites are those,
          having the SBML tag \code{boundaryCondition} set to \code{TRUE}.
          If set to \code{FALSE}, external metabolites are recogized via
          \code{extMetFlag}

    \item usage of class \code{sybilLog} in \code{readSBMLmod()} with four new
          parameters. See \code{help("readSBMLmod")} for details.
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.2 2011-07-20}{
  \itemize{
    \item fixed a bug in \code{readSBMLmod()}:
          if an external metabolite id is removed and it occcurs in
          a reaction, that particular reaction is read correctly.
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.1 2011-07-06}{
  \itemize{
    \item new boolean parameter \code{ignoreNoAn} in \code{readSBMLmod()}:
          if set to \code{TRUE}, any notes and annotation fields in the
          \code{listOfReactions} in an SBML file are ignored
          (default: \code{TRUE}).

    \item unused metabolites having no non-zero entry in S, are removed from S
          and the metabolites list in instances of class \code{modelorg}
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.1.0 2011-05-17}{
  \itemize{
    \item updated documentation
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0.3 2011-05-11}{
  \itemize{
    \item added compatibility to \CRANpkg{sybil} version \eqn{\geq}{>=} 0.4.1:
          slot \code{subSys} is now of class \code{matrix}

    \item fixed a bug in \code{readSBMLmod()}: metablite id's differing
          only in the number of consecutive '\code{_}' are now recogized
          as different metabolite id's.

    \item \code{readSBMLmod()} can now read the Yeast model in the file
          \file{Sc_iND750_GlcMM.xml} from the Palsson group.
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0.2 2011-03-02}{
  \itemize{
    \item If the SBML file contains reactions with metabolites used as
          reactands and products at the same time, the stoichiometric
          coefficient for those metabolites is the sum of reactand and product.
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0.1 2011-03-01}{
  \itemize{
    \item fixed some regular expressions when parsing metabolite names and id's

    \item added argument \code{extMetFlag} which is used to identify external
          metabolites

    \item reads subsystem information also from annotations
  }
}


% ---------------------------------------------------------------------------- %
\section{Changes in version 1.0.0 2011-01-25}{
  \itemize{
    \item initial public release
  }
}


