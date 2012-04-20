#  readSBMLmod.R
#  FBA and friends with R.
#
#  Copyright (C) 2010-2012 Gabriel Gelius-Dietrich, Dpt. for Bioinformatics,
#  Institute for Informatics, Heinrich-Heine-University, Duesseldorf, Germany.
#  All right reserved.
#  Email: geliudie@uni-duesseldorf.de
#
#  This file is part of sybilSBML.
#
#  SyBiL is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  SyBiL is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with SyBiL.  If not, see <http://www.gnu.org/licenses/>.


################################################
# Function: readSBMLmod
#
#
# The function readSBMLmod() is inspired by the function
# readCbModel() contained in the COBRA Toolbox.
# The algorithm is basically the same.


readSBMLmod <- function(filename, description,
                        def_bnd = SYBIL_SETTINGS("MAXIMUM"),
                        checkrsbml = TRUE,
                        extMetFlag = "b",
                        bndCond = TRUE,
                        ignoreNoAn = FALSE,
                        mergeMet = TRUE,
                        balanceReact = TRUE,
                        remUnusedMetReact = TRUE,
                        singletonMet = FALSE,
                        deadEndMet = FALSE,
                        remMet = FALSE,
                        constrMet = FALSE,
                        tol = SYBIL_SETTINGS("TOLERANCE")) {


#------------------------------------------------------------------------------#
# open the model file

if ( file.exists(filename) == FALSE ) {
   stop("failed to open file ", sQuote(filename))
}

if (missing(description)) {
    mdesc <- filename
}
else {
    mdesc <- description
}


#------------------------------------------------------------------------------#
#                some functions we use to create the model                     #
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#

entryforS <- function(X) {

#------------------------------------------------------------------------------#
# X containes metabolite id's (an object of class "SpeciesReference"). The slot
# "species" contains metabolite id. The function entryforS() gets the array
# index of X from the vector "sybil::met_id(sbml)", which is the line number in
# the stoichiometric matrix S.
#------------------------------------------------------------------------------#

    n    <- length(X)
    si   <- rep(i, n)    # the current column (reaction)
    sj   <- integer(n)   # the row (metabolite)
    s_ji <- numeric(n)   # the stoichiometric coefficient

    remMet <- logical(n) # metabolite removed from the initial metabolite list?

    CURR_MET <- character(n) # metabolites in X

    t <- 0
    for (el in X) {
        t <- t + 1

        # This is possible, because the metabolite id's are unique.
        # Keep in mind!
        # The metabolite id's are removed from the metabolites list,
        # but not from the reactions list.

        CURR_MET[t] <- el@species
        if (isTRUE(mergeMet)) {
            met_indCURR <- match(CURR_MET[t], CURR_MET[-t])
        }
        else {
            met_indCURR <- NA
        }

        if (is.na(met_indCURR)) {
            sj[t]     <- match(el@species, met_id_tmp)         # the row number
            s_ji[t]   <- el@stoichiometry
            remMet[t] <- ifelse(is.na(sj[t]), FALSE, TRUE)

        }
        else {
            remMet[t] <- FALSE
            s_ji[met_indCURR] <- s_ji[met_indCURR] + el@stoichiometry
            msg <- paste("reaction no.", i, dQuote(react_id_tmp[i]),
                         "metabolite no.", t, dQuote(CURR_MET[t]),
                         "was merged")
            warning(msg, call. = FALSE)
        }

#         if (is.na(sj[t])) {                       # if the current reaction is an
#             return(FALSE)                         # exchange reaction, sj[t] will
#         }                                         # be NA. So we'll leave s_ij[t]
#         else {                                    # at zero.
#             s_ji[t] <- el@stoichiometry           # the stoichiometric coefficient
#         }
    }

    #metUnq <- unique(sj[remMet])

    return(list(sj = sj[remMet], si = si[remMet], s_ji = s_ji[remMet]))

}

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#

checkupplowbnd <- function(x) {

#------------------------------------------------------------------------------#
# Check the absolute value of a reaction bound (vmin, vmax). If it is larger
# than def_bnd, it will be replaced by def_bnd.
#------------------------------------------------------------------------------#

    if (abs(x) > def_bnd) {
        bound <- def_bnd
        if (x < 0) {
            bound <- bound * -1
        }
    }
    else {
        bound <- x
    }
    return(bound)

}

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#

fatal <- function(x, part, value) {

#------------------------------------------------------------------------------#
# A control function. For every part in the SBML file, the Id's are a mandatory
# argument. Here we check, if they are all brave and there.
#------------------------------------------------------------------------------#

  if (length(x) == 0) {
      stop("missing ", part, " ", value)
  }
  return(TRUE)

}

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#

formatSBMLid <- function(idstr) {

#------------------------------------------------------------------------------#
# beautify SBML id's
#------------------------------------------------------------------------------#

    idstr <- gsub("-DASH-",        "-",   idstr, fixed = TRUE)
    idstr <- gsub("_DASH_",        "-",   idstr, fixed = TRUE)
    #idstr <- gsub("_FSLASH_",      "/",   idstr, fixed = TRUE)
    #idstr <- gsub("_BSLASH_",      "\\",  idstr, fixed = TRUE)
    idstr <- gsub("_LPAREN_",      "(",   idstr, fixed = TRUE)
    idstr <- gsub("_RPAREN_",      ")",   idstr, fixed = TRUE)
    idstr <- gsub("_LSQBKT_",      "[",   idstr, fixed = TRUE)
    idstr <- gsub("_RSQBKT_",      "]",   idstr, fixed = TRUE)
    idstr <- gsub("_COMMA_",       ",",   idstr, fixed = TRUE)
    idstr <- gsub("_PERIOD_",      ".",   idstr, fixed = TRUE)
    idstr <- gsub("_APOS_",        "'",   idstr, fixed = TRUE)
    idstr <-  sub( "_e_?$",        "(e)", idstr)   # nicer formatting of exchange reactions
    idstr <- gsub("-",             "_",   idstr, fixed = TRUE)
    #idstr <- gsub("&amp;",         "&",   idstr, fixed = TRUE)
    #idstr <- gsub("&lt;",          "<",   idstr, fixed = TRUE)
    #idstr <- gsub("&gt;",          ">",   idstr, fixed = TRUE)
    #idstr <- gsub("&quot;",        "\"",  idstr, fixed = TRUE)

    return(idstr)
}

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#

parseNotesReact <- function(notes) {

#------------------------------------------------------------------------------#
# parse the notes field of the reactions
#------------------------------------------------------------------------------#

  if (regexpr("html:p", notes, fixed = TRUE) == -1) {
      tag <- "p"
  }
  else {
      tag <- "html:p"
  }

  split <- paste("<", tag, ">", sep = "")
  #split <- "\n"

  fields <- strsplit(notes, split, fixed = TRUE)
 # print(fields)

  start_tag  <- paste("<", tag, ">", sep = "")
  end_tag    <- paste("</", tag, ">", sep = "")
  regex      <- paste("^(?:[\\t]*\\Q", start_tag, "\\E)?", "(.*)", "\\Q", end_tag, "\\E", "(?s).*$", sep = "")
#  regex      <- paste("(.*)", end_tag, "(?s).*$", sep = "")
  #print(regex)

  fields_str <- sub(regex, "\\1", fields[[1]], perl = TRUE)
  #print(fields_str)

  subSyst   <- ""
  gpr       <- ""
  gene_rule <- NA

  for (j in 1:length(fields_str)) {
      if (grepl("GENE[_ ]?ASSOCIATION", fields_str[j])) {
      #if (charmatch("GENE", fields_str[j], nomatch = -1) != -1) {
          gpr <- sub("GENE[_ ]?ASSOCIATION: *", "", fields_str[j])
          gene_rule <- sybil:::.parseBoolean(gpr)
          #print(gene_rule)
      }
      if (charmatch("SUBSYSTEM", fields_str[j], nomatch = -1) != -1) {
          subSyst <- sub("SUBSYSTEM: *", "", fields_str[j])
          subSyst <- sub("^S_", "", subSyst, perl = TRUE)
          subSyst <- gsub("[_]+", " ", subSyst)
          if (nchar(subSyst) == 0) {
              subSyst <- "Exchange"
          }
          #print(subSyst)
      }
  }
  if (!is.list(gene_rule)) {
      gene_rule <- sybil:::.parseBoolean("")
  }

  return(list(sub_system = subSyst, genes = gene_rule$gene, rules = gene_rule$rule, gpr = gpr))

}

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#

unlistRsbml <- function(listed) {

#------------------------------------------------------------------------------#
# the following is a specific unlist function for rsbml
#------------------------------------------------------------------------------#

  return(unlist(listed, use.names = FALSE, recursive = TRUE))
}


#------------------------------------------------------------------------------#

# readTheModel <- function(name) {
#
# #------------------------------------------------------------------------------#
# # Read the SBML model from SBML file. Check rsbml version and run
# # corresponding command.
# #------------------------------------------------------------------------------#
#
#   rsbml_version = packageDescription("rsbml", fields = "Version")
#
#   if (compareVersion(rsbml_version, "1.6.0") <= 0) {
#       # for rsbml <= 1.6.0 (bioconductor 2.1)
#       goodModel <- rsbml::rsbml_read(name, validate = FALSE, quiet = TRUE)
#       return(goodModel)
#   }
#   if (compareVersion(rsbml_version, "1.6.0") > 0) {
#       # for rsbml >= 2.0.0 (bioconductor 2.2)
#       goodModel <- rsbml::rsbml_read(name, dom = FALSE)
#       return(goodModel)
#   }
#   else {
#       # this is for rsbml < 1.5.2 (maybe I should leave that out)
#       goodModel <- rsbml::rsbml_read(name, check = FALSE, validate = FALSE)
#       return(goodModel)
#   }
#
#
# }

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                            reading the model                                 #
#------------------------------------------------------------------------------#

message("reading the model ... ", appendLF = FALSE)

ModelTmp <- try(rsbml::rsbml_read(filename, dom = FALSE), silent = TRUE)

if (is(ModelTmp, "try-error")) {
    message("seems to be bad formated, trying to fix ... ", appendLF = FALSE)
    hackedModel <- sybilSBML:::.uglyHack(filename)
    ModelTmp <- rsbml::rsbml_read(hackedModel, dom = FALSE)
    unlink(hackedModel)
    remove(hackedModel)
}

message("OK")


#------------------------------------------------------------------------------#
#                              check the model                                 #
#------------------------------------------------------------------------------#

if (checkrsbml == TRUE) {
    message("checking the model ... ", appendLF = FALSE)
    check <- rsbml::rsbml_check(ModelTmp)                   # check for errors

    if (check == FALSE) {

        rsbmlProb <- rsbml::rsbml_problems(ModelTmp)
        if ((length(rsbml::errors(rsbmlProb)) != 0) ||
            (length(rsbml::fatals(rsbmlProb)) != 0)) {
            msg <- paste("FAILED review file", dQuote(filename),
                         "carefully, returning rsbml_problems object")
            warning(msg)
            return(rsbmlProb)
        }
        else {
            msg <- "found warnings concerning SBML, check them carefully ... "
            message(msg, appendLF = FALSE)
            lapply(rsbml::warns(rsbmlProb),
                   function(x) {
                       wmsg <- gsub("\\s", " ", x@msg)
                       tmsg <- paste("[",
                                     x@line,
                                     " ",
                                     x@column, "] ",
                                     wmsg, sep = "")
                       warning(tmsg, call. = FALSE)
                    }
            )
        }

    }
    message("OK")
}


#------------------------------------------------------------------------------#
#                           translate the model                                #
#------------------------------------------------------------------------------#

message("translating the model ... ", appendLF = FALSE)
Mod <- try(rsbml::rsbml_dom(ModelTmp), silent = TRUE)

if (is(Mod, "try-error")) {
    message("rsbml_dom failed, trying to fix ... ", appendLF = FALSE)
    hackedModel <- sybilSBML:::.uglyHack(filename, remapply = TRUE)
    ModelTmp <- rsbml::rsbml_read(hackedModel, dom = FALSE)
    unlink(hackedModel)
    remove(hackedModel)
    Mod2 <- try(rsbml::rsbml_dom(ModelTmp))
    if (is(Mod2, "try-error")) {
        rsbmlProb <- rsbml::rsbml_problems(ModelTmp)
        msg <- paste("FAILED review file", dQuote(filename),
                     "carefully, returning rsbml_problems object")
        warning(msg)
        return(rsbmlProb)
    }
    else {
        Mod <- Mod2
    }
}

mid   <- ifelse(length(rsbml::id(rsbml::model(Mod))) == 0,   filename, rsbml::id(rsbml::model(Mod)))
mname <- ifelse(length(rsbml::name(rsbml::model(Mod))) == 0, filename, rsbml::name(rsbml::model(Mod)))

sbml <- sybil::modelorg(mid, mname)              # S4 object of class modelorg

message("OK")


#------------------------------------------------------------------------------#
#                             model description                                #
#------------------------------------------------------------------------------#

if (mdesc == filename) {
    mdesc  <- sub("^(.+)\\.xml", "\\1", basename(filename))
}

sybil::mod_desc(sbml) <- mdesc


#------------------------------------------------------------------------------#
#                                   units                                      #
#------------------------------------------------------------------------------#

# I have to do this, but later (2007-08-14)


#------------------------------------------------------------------------------#
#                               compartments                                   #
#------------------------------------------------------------------------------#

mod_compart_tmp          <- lapply(rsbml::compartments(rsbml::model(Mod)), rsbml::id)
sybil::mod_compart(sbml) <- unlistRsbml(mod_compart_tmp)

# check weather all Id's are there.
sapply(sybil::mod_compart(sbml), function(x) fatal(x, part = " compartment ", value = "Id"))



#------------------------------------------------------------------------------#
#                           initial reactions list                             #
#------------------------------------------------------------------------------#
# check for the reactions slot

react_id_tmp <- lapply(rsbml::reactions(rsbml::model(Mod)), rsbml::id)  # all reaction id's
react_id_tmp <- unlistRsbml(react_id_tmp)

# check weather all Id's are there.
sapply(react_id_tmp, function(x) fatal(x, part = " reaction ", value = "Id"))

# number of reactions (columns)
numreact <- length(react_id_tmp)


#------------------------------------------------------------------------------#
#                         initial metabolites list                             #
#------------------------------------------------------------------------------#
# check for the metabolites slot

# position of all metabolites except external ones
#met_id_pos <- grep("_?[^b]$", sapply(rsbml::species(rsbml::model(Mod)), rsbml::id))

metSpIds <- lapply(rsbml::species(rsbml::model(Mod)), rsbml::id)
metSpIds <- unlistRsbml(metSpIds)

if (isTRUE(bndCond)) {
    metSpBnd <- sapply(rsbml::species(rsbml::model(Mod)), rsbml::boundaryCondition)
    met_id_pos <- !metSpBnd
}
else {
    # regular expression to identify external metabolites
    extMetRegEx <- paste("_", extMetFlag, "$", sep = "")

    met_id_pos  <- grep(extMetRegEx, metSpIds, invert = TRUE)
    #met_id_pos <- grep("_b$", sapply(rsbml::species(rsbml::model(Mod)), rsbml::id), invert = TRUE)
    #met_id_pos <- grep("(?!_b)$", sapply(rsbml::species(rsbml::model(Mod)), rsbml::id), perl = TRUE)
}

met_id_tmp <- metSpIds[met_id_pos]

sapply(met_id_tmp, function(x) fatal(x, part = " metabolite ", value = "Id"))

# number of metabolites
nummet <- length(met_id_tmp)


#------------------------------------------------------------------------------#
#                            reversibilities                                   #
#------------------------------------------------------------------------------#

# boolean vector with reversibilities, default = TRUE
react_rev_tmp <- lapply(rsbml::reactions(rsbml::model(Mod)), rsbml::reversible)
react_rev_tmp <- unlistRsbml(react_rev_tmp)          # all reaction id's
react_rev_tmp <- react_rev_tmp != 0


#------------------------------------------------------------------------------#
#                           data structures                                    #
#------------------------------------------------------------------------------#


#S <- matrix(0, nummet, numreact)
St <- Matrix::Matrix(0, nrow = nummet, ncol = numreact, sparse = TRUE)

lbnd <- numeric(numreact)        # v min
ubnd <- numeric(numreact)        # v max
ocof <- numeric(numreact)        # objective coefficients


#------------------------------------------------------------------------------#
#                       S matrix and constraints, gpr                          #
#------------------------------------------------------------------------------#

message("creating S and parsing constraints ... ", appendLF = FALSE)

# for the gpr stuff
subSys   <- character(numreact)
genes    <- list(numreact)
rules    <- character(numreact)
gpr      <- character(numreact)
#allGenes <- character(0)

# Only one entry, because if one reaction has a notes
# field, all others are supposed to have one. Otherwise
# the gpr stuff does not make sense.
hasNotes <- FALSE
hasAnnot <- FALSE

for (i in 1 : numreact) {

    # the notes field
    notes <- rsbml::reactions(rsbml::model(Mod))[[i]]@notes
    annot <- rsbml::reactions(rsbml::model(Mod))[[i]]@annotation

    if (length(notes) != 0) {

        hasNotes    <- TRUE
        notes_field <- parseNotesReact(notes)
        #print(notes_field)
        subSys[i]   <- notes_field$sub_system   # the reaction's sub system: glykolysis, TCA, ...
        genes[[i]]  <- notes_field$genes        # list of genes
        rules[i]    <- notes_field$rules        # rules
        gpr[i]      <- notes_field$gpr          # original gpr association
        #allGenes    <- c(allGenes, genes[[i]])

    }
    else {

        if (length(annot) != 0) {
            hasAnnot    <- TRUE
            pn <- regexpr("Pathway Name: [^<]+", annot, perl = TRUE)
            subSys[i] <- substr(annot, (pn+14), pn + ((attr(pn, "match.length"))-1))
        }

    }

    # Check here if reactants and products lists exist, same for the stoichiometry slot

    # Entries for S -- the reactants
    S_tmp <- entryforS(rsbml::reactions(rsbml::model(Mod))[[i]]@reactants)
    #print(S_tmp)
    if (is.list(S_tmp) == TRUE) {
        St[S_tmp$sj, i] <- (S_tmp$s_ji * -1)
        #St[S_tmp$sj, S_tmp$si] <- (S_tmp$s_ji * -1)
    }
# Check here if S_tmp is FALSE. Should only be the case in
# the products slot due to the exclusion of external metabolites.
# In that case, the current reaction must be an exchange reaction.

#    else {
#        print(rsbml::reactions(rsbml::model(Mod))[[i]]@id)
#        print(S_tmp)
#        stop("something is wrong here")
#    }

    # Entries for S -- the products
    S_tmp <- entryforS(rsbml::reactions(rsbml::model(Mod))[[i]]@products)

    if (is.list(S_tmp) == TRUE) {

        if (isTRUE(balanceReact)) {
            nnull <- St[S_tmp$sj, i] %in% 0
            St[S_tmp$sj, i] <- St[S_tmp$sj, i] + S_tmp$s_ji

            if ( any(nnull == FALSE) ) {
                msg <- paste("reaction no.", i,
                             dQuote(react_id_tmp[i]), sum(!nnull),
                             ngettext(sum(!nnull),
                                      "metabolite was balanced",
                                      "metabolites were balanced:\n\t"),
                             paste(dQuote(met_id_tmp[S_tmp$sj[!nnull]]),
                                   collapse = "\n\t "))
                warning(msg, call. = FALSE)
            }

        }
        else {
            St[S_tmp$sj, i] <- S_tmp$s_ji
        }
    }
#    else {
#        print(rsbml::reactions(rsbml::model(Mod))[[i]]@id)
#        print(S_tmp)
#        stop("something is wrong here")
#    }

    # the constraints
    parm <- try((rsbml::model(Mod)@reactions[[i]]@kineticLaw)@parameters, silent = TRUE)
    if (!is(parm, "try-error")) {
#     parm <- (rsbml::model(Mod)@reactions[[i]]@kineticLaw)@parameters
#     if (length(parm) != 0) {
        for (p in parm) {

            if (p@id == "LOWER_BOUND") {
                #sbml@lowbnd[i] <- -1 * checkupplowbnd(p@value)
                #sbml@lowbnd[i] <- checkupplowbnd(p@value)
                lbnd[i] <- checkupplowbnd(p@value)
            }
            if (p@id == "UPPER_BOUND") {
                #sbml@uppbnd[i] <- checkupplowbnd(p@value)
                ubnd[i] <- checkupplowbnd(p@value)
            }
            if (p@id == "OBJECTIVE_COEFFICIENT") {
                #sbml@obj_coef[i] <- p@value
                ocof[i] <- p@value
            }
            # flux value?   (sbml file)
            # reduced cost?  (sbml file)
        }

    }
    else {
        ubnd[i] <- def_bnd
        if (isTRUE(react_rev_tmp[i])) {
            lbnd[i] <- -1 * def_bnd
        }
        else {
            lbnd[i] <- 0
        }
        ocof[i] <- 0
    }

}


# ---------------------------------------------------------------------------- #
# search for unused metabolites and unused reactions

# binary version of stoichiometric matrix
#Stb <- St != 0
Stb <- abs(St) > tol

SKIP_METABOLITE   <- rowSums(Stb) != 0
SKIP_REACTION     <- colSums(Stb) != 0


if (isTRUE(remUnusedMetReact)) {
    did <- "and therefore removed from S:"
}
else {
    did <- "in S:"
}


# ---------------------------------------------------------------------------- #
# empty rows

if ( any(SKIP_METABOLITE == FALSE) ) {
    met_list  <- paste(dQuote(met_id_tmp[!SKIP_METABOLITE]),
                       collapse = "\n\t")
    nmet_list <- sum(!SKIP_METABOLITE)
    msg_part  <- paste("not used in any reaction", did)
    msg <- sprintf(ngettext(nmet_list,
                            "%d metabolite is %s %s",
                            "%d metabolites are %s\n\t%s"),
                   nmet_list, msg_part, met_list)
    warning(msg, call. = FALSE)
}


# ---------------------------------------------------------------------------- #
# empty columns

if ( any(SKIP_REACTION == FALSE) ) {
    react_list  <- paste(dQuote(react_id_tmp[!SKIP_REACTION]),
                       collapse = "\n\t")
    nreact_list <- sum(!SKIP_REACTION)
    msg_part    <- paste("not used", did)
    msg <- sprintf(ngettext(nreact_list,
                            "%d reaction is %s %s",
                            "%d reactions are %s\n\t%s"),
                   nreact_list, msg_part, react_list)
    warning(msg, call. = FALSE)
}

if (!isTRUE(remUnusedMetReact)) {
    SKIP_METABOLITE[!SKIP_METABOLITE] <- TRUE
    SKIP_REACTION[!SKIP_REACTION]     <- TRUE
}


# ---------------------------------------------------------------------------- #
# single metabolites

sing_met   <- rep(NA, nrow(St))
sing_react <- rep(NA, ncol(St))

if (isTRUE(singletonMet)) {

    message("identifying reactions containing single metabolites ... ", appendLF = FALSE)

	singleton <- sybil:::.singletonMetabolite(mat = Stb)

	sing_met[!singleton$smet]     <- FALSE
	sing_react[!singleton$sreact] <- FALSE
	sing_met[singleton$smet]      <- TRUE
	sing_react[singleton$sreact]  <- TRUE

    # singleton metabolites found?
	if (sum(singleton$smet) > 0) {

        if ( xor(isTRUE(constrMet), isTRUE(remMet)) ) {

            if (isTRUE(constrMet)) {
                # set to zero
                did_watm <- "identified"
                did_watr <- "constrained"
            }
            else {
                # remove
                SKIP_METABOLITE[singleton$smet] <- FALSE
                SKIP_REACTION[singleton$sreact] <- FALSE
                did_watm <- "removed"
                did_watr <- "removed"
            }

            met_list  <- paste(dQuote(met_id_tmp[singleton$smet]),
                               collapse = "\n\t")
            nmet_list <- sum(singleton$smet)
            react_list  <- paste(dQuote(react_id_tmp[singleton$sreact]),
                               collapse = "\n\t")
            nreact_list <- sum(singleton$sreact)

            msgm <- sprintf(ngettext(nmet_list,
                                    "%s %d singleton metabolite: %s",
                                    "%s %d singleton metabolites:\n\t%s"),
                            did_watm, nmet_list, met_list)

            msgr <- sprintf(ngettext(nreact_list,
                                    "%s %d reaction containing singleton metabolites: %s",
                                    "%s %d reactions containing singleton metabolites:\n\t%s"),
                            did_watr, nreact_list, react_list)

            #warning(paste(msgm, msgr, sep = "\n\t "))
            warning(msgm, call. = FALSE)
            warning(msgr, call. = FALSE)

        }
        else {

            met_list  <- paste(dQuote(met_id_tmp[singleton$smet]),
                               collapse = "\n\t")
            nmet_list <- sum(singleton$smet)
            msg <- sprintf(ngettext(nmet_list,
                                   "%d metabolite is singleton in S: %s",
                                   "%d metabolites are singletons in S:\n\t%s"),
                           nmet_list, met_list)
            warning(msg, call. = FALSE)
        }

    }
    else {
        message("nothing found ... ", appendLF = FALSE)
        sing_met   <- logical(nrow(St))
        sing_react <- logical(ncol(St))
    }
}


# ---------------------------------------------------------------------------- #
# dead end metabolites

de_met   <- rep(NA, nrow(St))
de_react <- rep(NA, ncol(St))

if (isTRUE(deadEndMet)) {

    message("identifying reactions containing dead end metabolites ... ", appendLF = FALSE)

	demr <- sybil:::.deadEndMetabolite(mat   = St,
					         		   lb    = lbnd,
							           exclM = sing_met,
							           exclR = sing_react,
							           tol   = tol)

	de_met[!demr$dem]   <- FALSE
	de_react[!demr$der] <- FALSE
	de_met[demr$dem]    <- TRUE
	de_react[demr$der]  <- TRUE

    # dead end metabolites found?
	if (sum(demr$dem) > 0) {

		if ( xor(isTRUE(constrMet), isTRUE(remMet)) ) {

            if (isTRUE(constrMet)) {
                # set to zero
                did_watm <- "identified"
                did_watr <- "constrained"
            }
            else {
                # remove
                SKIP_METABOLITE[demr$dem] <- FALSE
                SKIP_REACTION[demr$der]   <- FALSE
                did_watm <- "removed"
                did_watr <- "removed"
            }

            met_list  <- paste(dQuote(met_id_tmp[demr$dem]),
                               collapse = "\n\t")
            nmet_list <- sum(demr$dem)
            react_list  <- paste(dQuote(react_id_tmp[demr$der]),
                               collapse = "\n\t")
            nreact_list <- sum(demr$der)

            msgm <- sprintf(ngettext(nmet_list,
                                    "%s %d dead end metabolite: %s",
                                    "%s %d dead end metabolites:\n\t%s"),
                            did_watm, nmet_list, met_list)

            msgr <- sprintf(ngettext(nreact_list,
                                    "%s %d reaction containing dead end metabolites: %s",
                                    "%s %d reactions containing dead end metabolites:\n\t%s"),
                            did_watr, nreact_list, react_list)

            warning(msgm, call. = FALSE)
            warning(msgr, call. = FALSE)

		}
		else {

			met_list  <- paste(dQuote(met_id_tmp[demr$dem]),
							   collapse = "\n\t")
			nmet_list <- sum(demr$dem)
			msg <- sprintf(ngettext(nmet_list,
								   "%d dead end metabolite in S: %s",
								   "%d dead end metabolites in S:\n\t%s"),
						   nmet_list, met_list)
			warning(msg, call. = FALSE)
		}

    }
    else {
        message("nothing found ... ", appendLF = FALSE)
        de_met   <- logical(nrow(St))
        de_react <- logical(ncol(St))
    }
}

# ---------------------------------------------------------------------------- #
# S

St <- St[SKIP_METABOLITE, , drop = FALSE]
St <- St[ , SKIP_REACTION, drop = FALSE]

sybil::S(sbml) <- St
#remove(St)

#sbml@S <- S
#remove(S)

# NNZ <- nonZeroElements(S(sbml))
#
# Sne(sbml) <- NNZ$ne
# Sia(sbml) <- NNZ$ia
# Sja(sbml) <- NNZ$ja
# Sar(sbml) <- NNZ$ar


numreact <- sum(SKIP_REACTION)

sybil::met_num(sbml)   <- sum(SKIP_METABOLITE)
sybil::react_num(sbml) <- numreact

sybil::met_single(sbml)   <- sing_met[SKIP_METABOLITE]
sybil::react_single(sbml) <- sing_react[SKIP_REACTION]

sybil::met_de(sbml)   <- de_met[SKIP_METABOLITE]
sybil::react_de(sbml) <- de_react[SKIP_REACTION]

if (isTRUE(constrMet)) {
    lbnd[sing_react] <- 0
    ubnd[sing_react] <- 0
    lbnd[de_react]   <- 0
    ubnd[de_react]   <- 0
}
else {}


sybil::lowbnd(sbml)    <- lbnd[SKIP_REACTION]
sybil::uppbnd(sbml)    <- ubnd[SKIP_REACTION]
sybil::obj_coef(sbml)  <- ocof[SKIP_REACTION]


message("OK")


#------------------------------------------------------------------------------#
#                        gene to reaction mapping                              #
#------------------------------------------------------------------------------#

if (isTRUE(ignoreNoAn)) {
    sybil::gprRules(sbml)   <- character(numreact)
    sybil::genes(sbml)      <- vector(mode = "list", length = numreact)
    sybil::gpr(sbml)        <- character(numreact)
    sybil::allGenes(sbml)   <- character(numreact)
    sybil::rxnGeneMat(sbml) <- Matrix::Matrix(FALSE, nrow = numreact, ncol = numreact, sparse = TRUE)
    sybil::subSys(sbml)     <- Matrix::Matrix(FALSE, nrow = numreact, ncol = 1, sparse = TRUE)
}
else {

    subSys <- subSys[SKIP_REACTION]
    genes  <- genes[SKIP_REACTION]
    rules  <- rules[SKIP_REACTION]
    gpr    <- gpr[SKIP_REACTION]

    if (isTRUE(hasNotes)) {
        message("GPR mapping ... ", appendLF = FALSE)

        #allGenes <- unique(allGenes)
        #allGenesTMP <- unique(allGenes)
        allGenesTMP <- unique(unlist(genes))
        temp <- nchar(allGenesTMP)
        allGenes <- allGenesTMP[which(temp != 0)]


        rxnGeneMat <- Matrix::Matrix(FALSE,
                                     nrow = numreact,
                                     ncol = length(allGenes),
                                     sparse = TRUE)

        for (i in 1 : numreact) {

            if ( (length(genes[[i]] == 1)) && (genes[[i]] != "") ) {
                geneInd <- match(genes[[i]], allGenes)
                rxnGeneMat[i, geneInd] <- TRUE
    
                for (j in 1 : length(geneInd)) {
                    pat  <- paste("x(", j, ")", sep = "")
                    repl <- paste("x[", geneInd[j], "]", sep = "")
    
                    rules[i] <- gsub(pat, repl, rules[i], fixed = TRUE)
                }
            }
        }

        sybil::genes(sbml)      <- genes
        sybil::gpr(sbml)        <- gpr
        sybil::allGenes(sbml)   <- allGenes
        sybil::gprRules(sbml)   <- rules
        sybil::rxnGeneMat(sbml) <- rxnGeneMat
        #sybil::subSys(sbml)     <- subSys
        sybil::subSys(sbml)     <- sybil:::.prepareSubSysMatrix(subSys, numreact)

        #sbml@gprRules <- rules
        #sbml@genes <- genes
        #sbml@gpr <- gpr
        #sbml@allGenes <- allGenes
        #sbml@subSys <- subSys

        message("OK")
    }
    else {
        sybil::rxnGeneMat(sbml) <- Matrix::Matrix(NA, nrow = 0, ncol = 0)
        if (isTRUE(hasAnnot)) {
            #subSys(sbml)     <- subSys
            sybil::subSys(sbml) <- sybil:::.prepareSubSysMatrix(subSys, numreact)
        }
        else {
            sybil::subSys(sbml) <- Matrix::Matrix(FALSE,
                                                  nrow = numreact,
                                                  ncol = 1,
                                                  sparse = TRUE)
        }
    }

}


#------------------------------------------------------------------------------#
#                              reaction id's                                   #
#------------------------------------------------------------------------------#
message("cleaning up ... ", appendLF = FALSE)

react_id_tmp   <- sub( "^R[_]+",        "", react_id_tmp[SKIP_REACTION])   # remove the leading R_
#react_id_tmp   <- gsub("_LPAREN_",      "(",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_RPAREN_",      ")",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_LSQBKT_",      "[",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_RSQBKT_",      "]",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_COMMA_",       ",",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_APOS_",        "'",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- gsub("_DASH_",        "-",   react_id_tmp, fixed = TRUE)
#react_id_tmp   <- sub( "_e_?$",         "(e)", react_id_tmp)   # nicer formatting of exchange reactions
#sybil::react_id(sbml) <- gsub("-",      "_",   react_id_tmp, fixed = TRUE)
sybil::react_id(sbml) <- formatSBMLid(react_id_tmp)


#------------------------------------------------------------------------------#
#                             reaction names                                   #
#------------------------------------------------------------------------------#

react_name_tmp   <- lapply(rsbml::reactions(rsbml::model(Mod)), rsbml::name)[SKIP_REACTION]
react_name_tmp[lapply(react_name_tmp, length) == 0] <- ""
react_name_tmp   <- unlistRsbml(react_name_tmp)
react_name_tmp   <- sub( "^R[_]+", "",  react_name_tmp)
react_name_tmp   <- gsub("[_]+",   " ", react_name_tmp)
react_name_tmp   <- sub( "\\s+$",  "",  react_name_tmp, perl = TRUE)
sybil::react_name(sbml) <- react_name_tmp


#------------------------------------------------------------------------------#
#                             metabolite id's                                  #
#------------------------------------------------------------------------------#

met_id_tmp     <- sub( "^[MSms]_+", "", met_id_tmp[SKIP_METABOLITE])   # remove the leading M_ or S_
#met_id_tmp     <- gsub( "[_]+",           "_",    met_id_tmp)
met_id_tmp     <- sub( "_([A-Za-z0-9]+)$",     "[\\1]", met_id_tmp)   # put the compartment id into square brackets
sybil::met_id(sbml) <- gsub("-",      "_",    met_id_tmp, fixed = TRUE)
#sybil::met_id(sbml) <- gsub("[-_]+",  "_",     met_id_tmp)


#------------------------------------------------------------------------------#
#                        metabolite compartments                               #
#------------------------------------------------------------------------------#

met_comp_tmp   <- lapply(rsbml::species(rsbml::model(Mod))[met_id_pos], rsbml::compartment)[SKIP_METABOLITE]

# check weather all Id's are there.
sapply(met_comp_tmp, function(x) fatal(x, part = " metabolite ", value = "compartment"))

sybil::met_comp(sbml) <- match(met_comp_tmp, sybil::mod_compart(sbml))


#------------------------------------------------------------------------------#
#                            metabolite names                                  #
#------------------------------------------------------------------------------#

met_name_tmp   <- lapply(rsbml::species(rsbml::model(Mod))[met_id_pos], rsbml::name)[SKIP_METABOLITE]
met_name_tmp[lapply(met_name_tmp, length) == 0] <- ""
met_name_tmp   <- unlistRsbml(met_name_tmp)
met_name_tmp   <- sub( "^[MS]?[_]+", "",  met_name_tmp)
met_name_tmp   <- gsub("[-_]+",   "-", met_name_tmp)
met_name_tmp   <- sub("-$",       "",  met_name_tmp)
met_name_tmp   <- sub( "\\s+$",   "",  met_name_tmp, perl = TRUE)
sybil::met_name(sbml) <- met_name_tmp


#------------------------------------------------------------------------------#
#                             right hand site                                  #
#------------------------------------------------------------------------------#

# add a right hand site to the model, all zeros - the null space: Sv = 0

sybil::rhs(sbml) <- numeric(sybil::met_num(sbml))

#------------------------------------------------------------------------------#
#                             check reversibilities                            #
#------------------------------------------------------------------------------#

# check up with the matlab version
# check the reversibilities
react_rev_tmp <- react_rev_tmp[SKIP_REACTION]
isrev <- which(sybil::lowbnd(sbml) < 0 & sybil::uppbnd(sbml) > 0)
#print(isrev)
react_rev_tmp[isrev]   <- TRUE
sybil::react_rev(sbml) <- react_rev_tmp

message("OK")


#------------------------------------------------------------------------------#
#                               validate the model                             #
#------------------------------------------------------------------------------#

message("validating object ... ", appendLF = FALSE)

check <- validObject(sbml, test = TRUE)

if (check != TRUE) {
    msg <- paste("Validity check failed:", check, sep = "\n    ")
    warning(msg)
}

message("OK")


#------------------------------------------------------------------------------#
#                                return the model                              #
#------------------------------------------------------------------------------#

# Returns sbml, an object of the class modelorg
return(sbml)

}
