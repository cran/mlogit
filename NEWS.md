# mlogit 0.4-1

* the Cracker, Catsup and Car data set are back in mlogit since AER,
  flexmix and mlogitBMA run examples based on them.

* the alt vector in the index is now carrefully checked in case of
  alternative subseting or reference level change.

# mlogit 0.4-0

* the main vignette is improved, writen in markdown and now and split
  by sections

* the Exercises vignette is splited and is now writen in markdown

* importantly, the Cholesky matrix is now coerced to a vector by rows
  and note by columns.
  
* the mlogit function was checked and improved.

* implementation of the computation of the standard deviations of the
  covariance matrix of the random parameters, using the delta method.
  
* some data sets are removed 

# mlogit 0.3-0

## New features

* zbu and zbt distributions are added : these are one-parameter
  distributions for which the lower bond is 0,

* a logsum function is provided to compute the log-sum or the
  inclusive utility of a random utility model,

* group-hetheroscedastic model can be estimated by setting the
  relevant covariates in the 4th part of the formula,

* the linear predictor is now returned by mlogit,

* correlation can still be a boolean, but can also be a character
  vector if one wants that a subset of the random parameters being
  correlated.

## data sets

* the RiskyTransport data set (used in the vignette to illustrate the
  estimation of the mixed logit model

* the NOx data set (used in the vignette to illustrate the estimation
  of the multinomial and group-heteroscedastic logit model),

* the JapaneseFDI data set (used in the vignette to illustrate the
  estimation of the nested logit model)
  
## vignette

A new vignette called mlogit2 is added ; this is the draft version of
an article submitted to the Journal of Statistical Software ; it is
less exhaustive, but better writen thant the original mlogit vignette.

## miscellanous

* the id series (one observation per choice situation) was badly
  constructed, it is now fixed

* the levels of the choice variable are now equalized to the those of
  the alt variable, allowing the case were some alternatives are never
  chosen

* mlogit is now able to estimate models with singular matrix of
  covariates. At the end of model.matrix.mformula, the linear
  dependent columns of X are removed

* there was a bug in the triangular distribution which is now fixed

* bug in the effects method fixed


# mlogit  0.2-4

* the list of primes used to generate halton sequences was too
	short, its length has been increased

* halton sequences where used to estimate mixed logit even for
	the default value of halton (NULL), this has been fixed

* the contribution of each observation to the gradient is not
	returned as the 'gradient' element of mlogit objects

* the distributions are now checked for rpar and an error is
	returned in case of unknown distribution

# mlogit 0.2-3
	
* some sys.frame() changed to parent.frame() 

# mlogit 0.2-2

* ranked-order models can be now estimated ; a new argument called
	'ranked' is introduced in mlogit.data which performs the relevant
	transformation of the data.frame. The estimated model is then a
	standard multinomial logit model

* multinomial probit model is now estimated by setting the new
	probit arguments to TRUE

* for the mixed logit model, different draws are now used for each
	observation

* a predict method is now available for mlogit objects

* a coef method is added which removes the fixed argument

* constPar can now be a named numeric vector. In this case, default
	starting values are changed according to constPar

* the vcov method for mlogit objects is greatly enhanced.

* mlogit objects now have two elements which indicate the fitted
	probabilities : fitted is the estimated probability for the
	outcome and probabilities contains the fitted probabilities for
	all the alternatives

* mentions to 'alt' in the names of the effects is canceled ;
	moreover, the intercepts are now called altname:('intercept')

* a 'choice' attribute is added to mlogit.data objects

* an effects method is provided, which computes the marginal
	effects of a covariate

# mlogit 0.2-1

* all the rda files are now compressed

# mlogit 0.2-0

* all the models could normally be estimated on unbalanced data

* the three tests are added, i.e. a new scoretest function and
	specific methods for waldtest and lrtest from the lmtest package

* the model.matrix method for mlogit objects is now exported

# mlogit 0.1-8

* mFormula modified so that models can be updated

* likelihood has been rewriten for the heteroscedastic logit
	model, the computation is now much faster

* nested logit models with overlapping nests are now supported;
	nests = "pcl" enables the estimation of the pair combinatorial
	logit model

* the norm argument is added to rpar

* the logLik argument is now of class logLik

* mlogit.data is modified so that an id argument can be used with
	data in long shape

* the argument of mlogit.data used to define longitudinal data is
	now called id.var

* mlogit.lnls is corrected so that the estimation of multinomial
	models can handle unbalanced data (pb with Reduce)
 
* the three tests are temporary removed
	
# mlogit 0.1-7

* a bug in mFormula (effects vs variable) is fixed

# mlogit 0.1-6

* a third part of the formula is added : it concerns alternative
	specific variables with alternative specific coefficients

* improved presentation for the Fishing dataset.

* a bug (forgotten drop = FALSE) corrected in
	model.matrix.mFormula

* Electricity and ModeCanada datasets are added

# mlogit 0.1-5

* if the choice variable is not an ordered factor, use as.factor()
	instead of class() <- "factor"

* cov.mlogit, cor.mlogit, rpar , med, rg, stdev, mean functions
	are added to extract and analyse random coefficients.
	
* a panel argument is added to mlogit so that mixed models with
	repeated observation can be estimated using panel methods or not

* a problem with the weights argument is fixed

* the estimation of nested logit models with a unique elasticity
	is now possible using un.nest.el = TRUE

* the estimation of nested logit models can now be done with or
	without normalization depending on the value of the argument
	unscaled

# mlogit 0.1-4

* mlogit didn't work when the dependent variable was an ordered
	factor in a "wide-shaped" data.frame.

* the reflevel argument didn't work any more in version 0.1-3.


# mlogit 0.1-3

* major change, most of the package has been rewriten

* it is now possible to estimate heteroscedastic, nested and mixed
	effects logit model

* the package doesn't depend any more on maxLik but a specific
	optimization function is provided for efficiency reason

# mlogit 0.1-2

* robust inference is provided with meat and estfunc methods
	defined for mlogit models.

* subset argument is added to mlogit so that the model may be
	estimated on a subset of alternatives.

* reflevel argument is added to mlogit which defines the base
	alternative.

* hmftest implements the Hausman McFadden test for the IIA hypothesis.

* mlogit.data function has been rewriten. It now use the reshape
  function.

* logitform class is provided to describe a logit model: update,
	model.matrix and model.frame methods are available.

# mlogit 0.1-1

# mlogit 0.1-0



	

	
	
	




	
