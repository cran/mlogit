#' Stated preference data for the choice of electricity suppliers
#' 
#' A sample of 2308 households in the United States
#' 
#' @name Electricity
#' @docType data
#' @format A dataframe containing :
#'
#' - choice: the choice of the individual, one of 1, 2, 3, 4,
#' - id: the individual index,
#' - pfi: fixed price at a stated cents per kWh, with the price
#' varying over suppliers and experiments, for scenario i=(1, 2, 3,
#' 4),
#' - cli: the length of contract that the supplier offered, in years
#' (such as 1 year or 5 years.) During this contract period, the
#' supplier guaranteed the prices and the buyer would have to pay a
#' penalty if he/she switched to another supplier. The supplier could
#' offer no contract in which case either side could stop the
#' agreement at any time. This is recorded as a contract length of 0,
#' - loci: is the supplier a local company,
#' - wki: is the supplier a well-known company,
#' - todi: a time-of-day rate under which the price is 11 cents per
#' kWh from 8am to 8pm and 5 cents per kWh from 8pm to 8am. These TOD
#' prices did not vary over suppliers or experiments: whenever the
#' supplier was said to offer TOD, the prices were stated as above.
#' - seasi: a seasonal rate under which the price is 10 cents per kWh
#' in the summer, 8 cents per kWh in the winter, and 6 cents per kWh
#' in the spring and fall. Like TOD rates, these prices did not vary.
#' Note that the price is for the electricity only, not transmission
#' and distribution, which is supplied by the local regulated utility.
#' 
#' @source
#' [Kenneth Train's home page](http://elsa.berkeley.edu/~train/).
#' @references
#' \insertRef{HUBE:TRAI:00}{mlogit}
#' 
#' \insertRef{REVE:TRAI:01}{mlogit}
#' @keywords datasets
#' @importFrom Rdpack reprompt

NULL


#' Choice of Fishing Mode
#' 
#' A sample of 1182 individuals in the United-States for the choice of
#' 4 alternative fishing modes.
#' 
#' @name Fishing
#' @docType data
#' @format A dataframe containing :
#' 
#' - mode: recreation mode choice, one of : beach, pier, boat and
#' charter,
#' - price.beach: price for beach mode
#' - price.pier: price for pier mode,
#' - price.boat: price for private boat mode,
#' - price.charter: price for charter boat mode,
#' - catch.beach: catch rate for beach mode,
#' - catch.pier: catch rate for pier mode,
#' - catch.boat: catch rate for private boat mode,
#' - catch.charter: catch rate for charter boat mode,
#' - income: monthly income,
#' 
#' @source
#' \insertRef{CAME:TRIV:05}{mlogit}
#' @references
#' \insertRef{HERR:KLIN:99}{mlogit}
#' @keywords datasets
NULL

#' Ranked data for gaming platforms
#' 
#' A sample of 91 Dutch individuals
#' 
#' The data are also provided in long format (use in this case
#' `data(Game2)`. In this case, the alternative and the choice
#' situation are respectively indicated in the `platform` and
#' `chid` variables.
#' 
#' @name Game
#' @aliases Game Game2
#' @docType data
#' @format A dataframe containing :
#' 
#' - ch.Platform: where `platform` is one of `Xbox`,
#' `PlayStation`, `PSPortable`, `GameCube`,
#' `GameBoy` and `PC`. This variables contain the ranking of
#' the platforms from 1 to 6,
#' - own.Platform: these 6 variables are dummies which indicate
#' whether the given plaform is already owned by the respondent,
#' - age: the age of the respondent,
#' - hours: hours per week spent on gaming.,
#' 
#' @source
#' [Journal of Applied Econometrics data archive](http://jae.wiley.com/jae/).
#' @references
#' \insertRef{FOK:PAAP:VAND:12}{mlogit}
#' @keywords datasets
NULL

#' Heating and Cooling System Choice in Newly Built Houses in California
#'
#' A sample of 250 Californian households
#' 
#' @name HC
#' @docType data
#' @format A dataframe containing :
#' 
#' - depvar: heating system, one of `gcc` (gas central heat with
#' cooling), `ecc` (electric central resistence heat with cooling), `erc`
#' (electric room resistence heat with cooling), `hpc` (electric heat
#' pump which provides cooling also), `gc` (gas central heat without
#' cooling), `ec` (electric central resistence heat without cooling), `er`
#' (electric room resistence heat without cooling),
#' - ich.z: installation cost of the heating portion of the
#' system,
#' - icca: installation cost for cooling,
#' - och.z: operating cost for the heating portion of the system,
#' - occa: operating cost for cooling,
#' - income: annual income of the household.
#' 
#' @source
#' [Kenneth Train's home page](http://elsa.berkeley.edu/~train/).
#' @keywords datasets
NULL

#' Heating System Choice in California Houses
#' 
#' A sample of 900 Californian households#'
#' 
#' @name Heating
#' @docType data
#' @format A dataframe containing:
#'
#' - idcase: id,
#' - depvar: heating system, one of gc (gas central), gr (gas room), ec
#' (electric central), er (electric room), hp (heat pump),
#' - ic.z: installation cost for heating system z (defined for the
#' 5 heating systems),
#' - oc.z: annual operating cost for heating system z (defined for
#' the 5 heating systems),
#' - pb.z: ratio oc.z/ic.z ,
#' - income: annual income of the household,
#' - agehed: age of the household head
#' - rooms: numbers of rooms in the house,
#'
#' @source
#' [Kenneth Train's home page](http://elsa.berkeley.edu/~train/).
#' @keywords datasets
NULL

#' Japanese Foreign Direct Investment in European Regions
#' 
#' A sample of 452 Japanese production units in Europe
#' #' 
#' @name JapaneseFDI
#' @docType data
#' @format A dataframe containing :
#'
#' - firm: the investment id,
#' - country: the country,
#' - region: the region (nuts1 nomenclature),
#' - choice: a dummy indicating the chosen region ,
#' - choice.c: the chosen country,
#' - wage: wage rate in the region,
#' - unemp: unemployment rate in the region,
#' - elig: is the country eligible to european subsidies,
#' - area: the area of the region,
#' - scrate:  social charge rate (country level),
#' - ctaxrate: corporate tax rate (country level),
#' - gdp: regional gdp,
#' - harris: harris' market potential,
#' - krugman: krugman's market potential,
#' - domind: domestic industry count,
#' - japind: japan industry count,
#' - network: network count.
#' 
#' @references
#' \insertRef{HEAD:MAYE:04}{mlogit}
#' @source kindly provided by Thierry Mayer
#' @keywords datasets
NULL

#' Mode Choice
#' 
#' A sample of 453 individuals for 4 transport modes.
#' 
#' @name Mode
#' @docType data
#' @format A dataframe containing :
#' 
#' - choice: one of car, carpool, bus or rail,
#' - cost.z: cost of mode z,
#' - time.z: time of mode z.
#' 
#' @source
#' [Kenneth Train's home page](http://elsa.berkeley.edu/~train/).
#' @keywords datasets
NULL

#' Mode Choice for the Montreal-Toronto Corridor
#'
#' A sample of 3880 travellers for the Montreal-Toronto corridor
#' 
#' @name ModeCanada
#' @docType data
#' @format A dataframe containing
#' 
#' - case: the individual index,
#' - alt: the alternative, one of train, car, bus and air,
#' - choice: one if the mode is chosen, zero otherwise,
#' - cost: monetary cost,
#' - ivt: in vehicule time,
#' - ovt: out vehicule time,
#' - frequency: frequency,
#' - income: income,
#' - urban: urban,
#' - noalt: the number of alternatives available.
#' 
#' @source
#' kindly provided by S. Koppelman
#' @references
#' \insertRef{BHAT:95}{mlogit}
#' 
#' \insertRef{KOPP:WEN:00}{mlogit}
#' 
#' \insertRef{WEN:KOPP:01}{mlogit}
#' @keywords datasets
#' @examples
#' data("ModeCanada", package = "mlogit")
#' bususers <- with(ModeCanada, case[choice == 1 & alt == "bus"])
#' ModeCanada <- subset(ModeCanada, ! case %in% bususers)
#' ModeCanada <- subset(ModeCanada, noalt == 4)
#' ModeCanada <- subset(ModeCanada, alt != "bus")
#' ModeCanada$alt <- ModeCanada$alt[drop = TRUE]
#' KoppWen00 <- mlogit.data(ModeCanada, shape='long', chid.var = 'case',
#'                          alt.var = 'alt', choice='choice',
#'                          drop.index=TRUE)
#' pcl <- mlogit(choice~freq+cost+ivt+ovt, KoppWen00, reflevel='car',
#'               nests='pcl', constPar=c('iv:train.air'))
#' 
NULL

#' Technologies to reduce NOx emissions
#' 
#' A sample of 632 American production units
#' 
#' @name NOx
#' @docType data
#' @format A dataframe containing:
#' 
#' - chid: the plant id,
#' - alt: the alternative,
#' - id: the owner id,
#' - choice: the
#' chosen alternative,
#' - available: a dummy indicating that the alternative is
#' available,
#' - env: the regulatory environment, one of `'regulated'`,
#' `'deregulated'` and `'public'`,
#' - post: dummy for post-combustion polution control technology,
#' - cm: dummy for combustion modification technology,
#' - lnb: dummy for low NOx burners technology,
#' - age: age of the plant (in deviation from the mean age).,
#' - vcost: variable cost,
#' - kcost: capital cost.
#' 
#' @references
#' \insertRef{FOWL:10}{mlogit}
#' @source
#' [American Economic Association data archive](http://aeaweb.org/aer/).
#' @keywords datasets
NULL

#' Risky Transportation Choices
#' 
#' 1793 choices by 561 individuals of a transport mode at Freetwon
#' airport
#' 
#' 
#' @name RiskyTransport
#' @docType data
#' @format A dataframe containing:
#'
#' - id: individual id,
#' - choice: 1 for the chosen mode,
#' - mode: one of `Helicopter`,`WaterTaxi`, `Ferry,
#' and `Hovercraft`,
#' - cost: the generalised cost of the transport mode,
#' - risk: the fatality rate, numbers of death per 100,000 trips,
#' - weight: weights,
#' - seats: ,
#' - noise: ,
#' - crowdness: ,
#' - convloc: ,
#' - clientele: ,
#' - chid: choice situation id,
#' - african: `yes` if born in Africa, `no` otherwise,
#' - lifeExp: declared life expectancy,
#' - dwage: declared hourly wage,
#' - iwage: imputed hourly wage,
#' - educ: level of education, one of `low` and `high`,
#' - fatalism: self-ranking of the degree of fatalism,
#' - gender: gender, one of `female` and `male`,
#' - age: age,
#' - haveChildren: `yes` if the traveler has children,
#' `no` otherwise,
#' - swim: `yes` if the traveler knows how to swim, `no,
#' otherwise.
#'
#' @references
#' \insertRef{LEON:MIGU:17}{mlogit} 
#' @source
#' [American Economic Association data archive](http://aeaweb.org/aer/).
#' @keywords datasets
NULL

# Gianmarco, Leon and Edward Miguel (2017)
#     "Risky Transportation Choices and the Value of a Statistical Life",
#     *American Economic Journal: Applied Economics*, **9(1)**,
#     202-28.


#' Stated Preferences for Train Traveling
#' 
#' A sample of 235 Dutch individuals facing 2929 choice situations
#' 
#' @name Train
#' @docType data
#' @format A dataframe containing:
#' 
#' - id: individual identifient,
#' - choiceid: choice identifient,
#' - choice: one of 'A' or 'B',
#' - price_z: price of proposition z (z = 'A', 'B') in cents of
#' guilders,
#' - time_z: travel time of proposition z (z = 'A', 'B') in
#' minutes,
#' - comfort_z: comfort of proposition z (z = 'A', 'B'), 0, 1 or
#' 2 in decreasing comfort order,
#' - change_z: number of changes for proposition z (z = 'A', 'B').
#' @source
#' [Journal of Applied Econometrics data archive](http://jae.wiley.com/jae/).
#' 
#' @references
#' \insertRef{BENA:BOLD:BRAD:93}{mlogit}
#' 
#' \insertRef{MEIJ:ROUW:06}{mlogit}
#' @keywords datasets
NULL

#' Choice of Brand for Crakers
#' 
#' a sample of 3292 individualscross-section
#' 
#' @name Cracker
#' @docType data
#' @format A dataframe containing :
#'
#' - id: individuals identifiers,
#' - choice: one of sunshine, keebler, nabisco, private,
#' - disp.z: is there a display for brand z?
#' - feat.z: is there a newspaper feature advertisement for brand z?
#' - price.z: price of brand z.
#' 
#' @source
#' [Journal of Business Economics and Statistics web site](https://www.amstat.org).
#' 
#' @references
#' \insertRef{DIPA:NAUF:CHIN:94}{mlogit}
#'
#' \insertRef{PAAP:FRAN:00}{mlogit}
#' @keywords datasets
NULL

#' mlogit package: estimation of random utility discrete choice models
#' by maximum likelihood
#'
#' mlogit provides a model description interface (enhanced
#' formula-data), a very versatile estimation function and a testing
#' infrastructure to deal with random utility models.
#'
#' @name mlogit-package
#' @docType package
#' @details For a gentle and comprehensive introduction to the
#'     package, see the package's vignettes.
NULL


#' Choice of Brand for Catsup
#' 
#' a sample of 2798 individuals
#' 
#' @name Catsup
#' @docType data
#' @format A dataframe containing :
#'
#' - id: individuals identifiers,
#' - choice: one of heinz41, heinz32, heinz28, hunts32,
#' - disp.z: is there a display for brand z ?
#' - feat.z: is there a newspaper feature advertisement for brand z?
#' - price.z: price of brand z.
#' @source
#' [Journal of Business Economics and Statistics web site](https://www.amstat.org).
#' 
#' @references
#' \insertRef{DIPA:NAUF:CHIN:94}{mlogit}
#' @keywords datasets
NULL

#' Stated Preferences for Car Choice
#' 
#' a sample of 4654 individuals

#' @name Car
#' @docType data
#' @format A dataframe containing :
#'
#' - choice: choice of a vehicule amoung 6 propositions, 
#' - college: college education?,
#' - hsg2: size of household greater than 2? 
#' - coml5: commulte lower than 5 miles a day?, 
#' - typez: body type, one of regcar (regular car), sportuv (sport utility vehicule), sportcar, stwagon (station wagon), truck, van, for each proposition z from 1 to 6, 
#' - fuelz: fuel for proposition z, one of gasoline, methanol, cng (compressed natural gas), electric.,
#' - pricez: price of vehicule divided by the logarithme of income,
#' - rangez: hundreds of miles vehicule can travel between refuelings/rechargings, 
#' - accz: acceleration, tens of seconds required to reach 30 mph from stop, 
#' - speedz: highest attainable speed in hundreds of mph, 
#' - pollutionz: tailpipe emissions as fraction of those for new gas vehicule, 
#' - sizez: 0 for a mini, 1 for a subcompact, 2 for a compact and 3 for a mid--size or large vehicule, 
#' - spacez: fraction of luggage space in comparable new gas vehicule, 
#' - costz: cost per mile of travel (tens of cents) : home recharging for electric vehicule, station refueling otherwise, 
#' - stationz: fraction of stations that can refuel/recharge vehicule.
#'
#' @source
#' [Journal of Applied Econometrics data archive](http://jae.wiley.com/jae/).
#' 
#' @references
#' \insertRef{MCFA:TRAI:00}{mlogit}
#' 
#' @keywords datasets
NULL

