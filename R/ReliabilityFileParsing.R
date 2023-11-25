#' Estimates system reliability for a specified system.
#' @import stringr
#' @param file A string with the file path containing the reliability diagram (see details)
#' @param priorList A list with named elements containing bsp object representing the prior for each component in the reliability diagram
#' @param dataList A list with named elements for each component or subsystem for which reliability was measured
#'
#' @return A named list of posterior bsp objects, one for each named component in the system
#' (if no data was provided the posterior will be equivalent to the provided prior)
#' @export
#'
#' @section file Details:
#' The file is a text file specifying your reliability block diagram. Parallel and series subsystems can be specified as in the following example:
#'
#' S(Engine, GasDelivery):GasPropulsion
#'
#' S(Motor, Batteries, Controller, Belt):ElectricPropulsion
#'
#' P(GasPropulsion, ElectricPropulsion):Propulsion
#'
#' S and P represent series and parallel relationships respectively.
#' The first line of this diagram can be read as
#' "Engine and GasDelivery are components in series
#' of a subsystem named GasPropulsion." For each component on the right hand
#' side of an equation, the prior of that component will be determined
#' by the posterior of its subcomponents. Data, if available, can be used to update its
#' posterior by adding the data to dataList (e.g. assigning a matrix to dataList$GasPropulsion).
#' For each component on the left hand side of the colon, priors can be defined
#' by adding a bsp object to priorList. If no prior is provided, a non-informative
#' prior will be used, with a warning.
#'
#' The following syntax conventions for the file will be enforced:
#'\itemize{
#' \item Each line begins with S( or P( indicating series or parallel, one relationship per line
#'
#'\item Components can be uppercase, lowercase, and include numbers, but
#'no special characters (names are case-sensitive).
#'
#'\item All relationships must be named by adding a colon and a label after the end of the relationship
#'
#' \item Two or more elements can be specified in a relationship
#'
#' \item A check will be done to ensure there are no circular dependencies among the named components.
#' For example, Propulsion cannot be used as an input to GasPropulsion because GasPropulsion
#' is an input to Propulsion
#'
#' \item Currently, all component names on the LHS of the diagram must be unique
#'}
#'
#'@section priorList Details:
#'The priorList contains the bsp prior objects used as objects for the components.
#'For example, if "valve" is a component appearing only on the left hand side of a relationship
#'in the reliability diagram, priorList$valve must contain a betaStaceyProcess object
#'representing the prior. These objects can be specified with the bsp() function.
#'Objects appearing on the right hand side of a relationship do not need a prior, as
#'their prior will be computed from the posteriors of the component pieces.
#'If a prior is not provided for a LHS component, a non-informative prior will be
#'created and used, with a warning.
#'
#'@section dataList Details:
#'A dataList element should be provided for any component on which data has been gathered.
#'Each dataList element (eg dataList$valve) should be a nx2 matrix with the first column
#'being the observed failurer time (or end or test) and the second column contains the censoring variable
#'0 - if right censored 1 - if fully observed.
#'
#'Any object for which data is provided will have a new bsp object computed as the
#'posterior of the prior and the provided data.
#'
#'
#' @examples
#' file <- textConnection("S(Engine, GasDelivery):GasPropulsion")
#' priorList <- list(Engine = bsp(1, .5, 1), GasDelivery=bsp(3,.6, 1))
#' dataList <- list(Engine=matrix(c(2, 1, 3, 0), byrow=TRUE, nrow=2),
#'     GasDelivery=matrix(c(5, 1, 6, 1), byrow=TRUE, nrow=2))
#' estimateSystemReliability(file, priorList, dataList)
#' 
estimateSystemReliability<-function(file, priorList, dataList){

  text=suppressWarnings(readLines(file))
  if(length(text)==0)stop("Text file contained 0 elements")
  text=stringr::str_replace_all(text, " ", "")

  #Checks to see if valid
  if(!all(stringr::str_sub(text, 1, 2) %in% c('p(', 'P(', 's(', 'S(')))stop("All expressions must start with P( or S(")


  #Makes sure the pattern after the opening parenthesis is #,#,#,....):#
  if(!(all(stringr::str_detect(stringr::str_sub(text, 3),
            pattern = "^([a-zA-Z0-9]+,)+[a-zA-Z0-9]+\\):[a-zA-Z0-9]+$"))))
    stop("One or more strings did not meet grammar requirements (eg. P(V1,V2,V3):P1)")

  #http://stla.github.io/stlapblog/posts/Numextract.html
  wordExtract <- function(string){
    middle=stringr::str_sub(stringr::str_extract(string,
                          '\\([a-zA-Z0-9,]+\\)'),2,-2)
    middleWords<-unlist(strsplit(middle, ','))
    endWord=stringr::str_sub(stringr::str_extract(string, ":[a-zA-Z0-9]+"), 2, -1)
    return(c(middleWords, endWord))
  }

  words<-wordExtract(text)

  ##The first time through we will see if it is logically defined
  #The 2nd time through we will perform the computations
  needsBuilt<-stringr::str_sub(stringr::str_extract(text, ":[a-zA-Z0-9]+"), 2, -1)

  if (length(needsBuilt)!=length(unique(needsBuilt)))stop("Non-unique definitions given to component (check right side of colon for non-unique values)")
  #Anything is ready if it is not on the RHS of a definition
  ready<-setdiff(words, needsBuilt)
  #Require list of priors
  #Require list of data

  ####Check to make sure all components only appear once on the LHS of the relations (enforce that data not used twice)
  counts<-table(words)
  countIndices<-names(counts)%in% names(table(needsBuilt))
  if(any(counts[countIndices]>2) | !all(counts[!countIndices]==1))
    stop("Component names must be unique on LHS of equations")

  done<-rep(F, length(text))
  i=0
  oneChange=FALSE
  while(any(!done)){
    i=i+1
    if (i > length(text) &!oneChange)stop("The diagram specified contains circular relationships that makes the model impossible")
    else if (i>length(text) &oneChange){
      i =1
      oneChange=FALSE
    }
    if (done[i]) next
    compNeeded <- wordExtract(text[i])
    Name<-compNeeded[length(compNeeded)]
    compNeeded<-compNeeded[-length(compNeeded)]
    if (all(compNeeded%in% ready)){
      ##We are ready to perform the computations needed in order to provide the prior for Name
      needsBuilt<-setdiff(needsBuilt, Name)
      ready<-c(ready, Name)
      oneChange=TRUE
      done[i]<-TRUE
    }
  }


  #2nd time through
  needsBuilt<-stringr::str_sub(stringr::str_extract(text, ":[a-zA-Z0-9]+"), 2, -1)
  #Anything is ready if it is not on the RHS of a definition
  ready<-setdiff(wordExtract(text), needsBuilt)
  posteriorList<-list()

  for (comp in ready){
    if(is.null(priorList[[comp]])){
      warning(paste0("No prior provided for ", comp, ". A non-informative prior was created"))
      if (is.null(dataList[[comp]]))stop(paste("No prior or data provided for component", comp))
      priorList[[comp]]<-createUninformativePrior(dataList[[comp]])
    }
    if (!is.null(dataList[[comp]])){
      posteriorList[[comp]]<-bspPosterior(priorList[[comp]], dataList[[comp]])
    }else{
      posteriorList[[comp]]<-priorList[[comp]]
    }
  }
  #Require list of priors
  #Require list of data
  done<-rep(F, length(text)) #Done matches up with the lines in variable text
  i=0
  oneChange=FALSE
  while(any(!done)){
    i=i+1
    if (i>length(text) &oneChange){
      i =1
      oneChange=FALSE
    }
    if (done[i]) next
    compNeeded <- wordExtract(text[i])
    merging_function=ifelse(substr(text[i],1,1)%in%c("S", "s"), E1E2_series, E1E2_parallel)
    Name<-compNeeded[length(compNeeded)]
    compNeeded<-compNeeded[-length(compNeeded)]
    if (all(compNeeded%in% ready)){
      ##We are ready to perform the computations needed in order to provide the prior for Name

      #This handles the case where there are more than two components in the list
      temp = posteriorList[[compNeeded[[1]]]]
      for (j in 2:length(compNeeded)){
        temp<-bspFromMoments(merging_function(temp, posteriorList[[compNeeded[j]]]))
        #d<<-temp
        #if(compNeeded[j]=="generator7")skldf
      }

      priorList[[Name]]<-temp

      ###As a placeholder, just create an uninformative prior for each subsystem
      #Remove these lines later
      #if(!is.null(dataList[[Name]])) data<-dataList[[Name]] else data<-matrix(c(1,1), byrow=T, nrow=T)
      #priorList[[Name]]<-createUninformativePrior(data)
      ####################


      needsBuilt<-setdiff(needsBuilt, Name)
      ready<-c(ready, Name)
      oneChange=TRUE
      done[i]<-TRUE
      if (!is.null(dataList[[Name]])){
        posteriorList[[Name]]<-bspPosterior(priorList[[Name]], dataList[[Name]])
      }else{ #if there's no data to augment it, just pass the prior through as the posterior
        posteriorList[[Name]]<-priorList[[Name]]
      }
    }
  }

  return(structure(posteriorList, class="bspPosteriorList"))

}


#' Creates a non-informative prior
#'
#' @param data The data matrix, which is used to ensure only the points necessary appear on the support
#'
#' @return A bsp object with 0 precision.
#' @export
#'
createUninformativePrior<-function(data){
  #mindata<-min(data[,1])
  maxdata<-max(data[,1])
  bsp(c(maxdata), centeringMeasure = c(.999), precision = 0)

}
